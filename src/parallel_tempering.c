#include <signal.h>
#include <gsl/gsl_sf.h>
#include <omp.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include "parallel_tempering.h"
#include "parallel_tempering_interaction.h"
#include "debug.h"

int run;
int dumpflag;

void set_beta(mcmc * m, double newbeta) {
	((parallel_tempering_mcmc *) m->additional_data)->beta = newbeta;
	((parallel_tempering_mcmc *) m->additional_data)->swapcount = 0;
}
double get_beta(const mcmc * m) {
	return ((parallel_tempering_mcmc *) m->additional_data)->beta;
}
void inc_swapcount(mcmc * m) {
	((parallel_tempering_mcmc *) m->additional_data)->swapcount++;
}
unsigned long get_swapcount(const mcmc * m) {
	return ((parallel_tempering_mcmc *) m->additional_data)->swapcount;
}

static void ctrl_c_handler(int signalnr) {
	printf("\nreceived Ctrl-C (%d). Stopping ... (please be patient)\n\n",
			signalnr);
	run = 0;
}
static void sigusr_handler(int signalnr) {
	printf("\nreceived SIGUSR (%d). Will dump at next opportunity.\n\n",
			signalnr);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);
	dumpflag = 1;
}

#include <time.h>

static int get_duration() {
	static clock_t stored = 0;
	clock_t new = stored;
	stored = clock();
	return stored - new;
}

static void print_current_positions(const mcmc ** sinmod, const int n_beta) {
	int i;
	printf("printing chain parameters: \n");
	for (i = 0; i < n_beta; i++) {
		printf("\tchain %d: swapped %lu times: ", i, get_swapcount(sinmod[i]));
		printf("\tchain %d: current %f: ", i, get_prob(sinmod[i]));
		dump_vectorln(get_params(sinmod[i]));
		printf("\tchain %d: best %f: ", i, get_prob_best(sinmod[i]));
		dump_vectorln(get_params_best(sinmod[i]));

	}
	fflush(stdout);
}

static void report(const mcmc ** sinmod, const int n_beta) {
	int i = 0;
	char buf[100];
	print_current_positions(sinmod, n_beta);
	printf("\nwriting out visited parameters ");
	while (1) {
		printf(".");
		sprintf(buf, "-chain%d", i);
		mcmc_dump_probabilities(sinmod[i], -1, buf);
		fflush(stdout);
#ifndef DUMP_ALL_CHAINS
		break;
#endif
		if (i > n_beta)
			break;
		i++;
	}
	printf("done.\n");
}

double equidistant_beta(const unsigned int i, const unsigned int n_beta,
		const double beta_0) {
	return beta_0 + i * (1 - beta_0) / (n_beta - 1);
}
double equidistant_temperature(const unsigned int i,
		const unsigned int n_beta, const double beta_0) {
	return 1 / (1 / beta_0 + i * (1 - 1 / beta_0) / (n_beta - 1));
}
double chebyshev_temperature(const unsigned int i,
		const unsigned int n_beta, const double beta_0) {
	return 1 / (1 / beta_0 + (1 - 1 / beta_0) / 2 * (1 - cos(i * M_PI / (n_beta
			- 1))));
}
double chebyshev_beta(const unsigned int i, const unsigned int n_beta,
		const double beta_0) {
	return beta_0 + (1 - beta_0) / 2 * (1 - cos(i * M_PI / (n_beta - 1)));
}

double equidistant_stepwidth(const unsigned int i,
		const unsigned int n_beta, const double beta_0) {
	return beta_0 + pow(i / (n_beta - 1), 2) * (1 - beta_0);
}

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Defines how the beta value should be assigned/distributed between the chains.
 *
 * beta = 1 / temperature.
 *
 * You can choose equidistant (linear) distribution of the temperature or beta.
 * Or, what often proves to be a good choice, you can use Chebyshev nodes
 * for the distribution of the temperature or beta.
 *
 * e.g.: BETA_DISTRIBUTION=equidistant_beta <br>
 * also available: equidistant_temperature, chebyshev_beta,
 * chebyshev_temperature, equidistant_stepwidth
 *
 */
#define BETA_DISTRIBUTION
#endif

static double get_chain_beta(unsigned int i, unsigned int n_beta, double beta_0) {
#ifndef BETA_DISTRIBUTION
#define BETA_DISTRIBUTION chebyshev_temperature
#endif
	if (n_beta == 1)
		return 1.0;
	/* this reverts the order so that beta(0) = 1.0. */
	return BETA_DISTRIBUTION(n_beta - i - 1, n_beta, beta_0);
}

static void analyse(mcmc ** sinmod, int n_beta, unsigned int n_swap);

void parallel_tempering(const char * params_filename,
		const char * data_filename, const int n_beta, const double beta_0,
		const unsigned long burn_in_iterations, const double rat_limit,
		const unsigned long iter_limit, const double mul, const int n_swap) {
	int n_par;
	int i;
	const char ** params_descr;
	mcmc ** sinmod;

	sinmod = (mcmc**) mem_calloc(n_beta, sizeof(mcmc*));
	assert(sinmod != NULL);

	printf("Initializing parallel tempering for %d chains\n", n_beta);
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, get_chain_beta(i, n_beta, beta_0));
		/* That is kind of stupid (duplicate execution) and could be optimized.
		 * not critical though. */
		sinmod[i] = mcmc_load_params(params_filename);
		if (i == 0)
			mcmc_load_data(sinmod[i], data_filename);
		else
			mcmc_reuse_data(sinmod[i], sinmod[0]);
		mcmc_check(sinmod[i]);
		sinmod[i]->additional_data
				= mem_malloc(sizeof(parallel_tempering_mcmc));
		set_beta(sinmod[i], get_chain_beta(i, n_beta, beta_0));
		printf("\tsteps: ");
		dump_vectorln(get_steps(sinmod[i]));
		mcmc_check(sinmod[i]);
	}
	params_descr = get_params_descr(sinmod[0]);
	n_par = get_n_par(sinmod[0]);

	/* here was code for removing *_results.dat files.
	 * this should be done externally before calling the program. */

	printf("Initializing models\n");
	fflush(stdout);
	for (i = 0; i < n_beta; i++) {
		calc_model(sinmod[i], NULL);
		mcmc_check(sinmod[i]);
	}
	printf("Starting markov chain calibration\n");
	fflush(stdout);
	markov_chain_calibrate(sinmod[0], burn_in_iterations, rat_limit,
			iter_limit, mul, DEFAULT_ADJUST_STEP);
	printf("You can update your parameters file to the initial steps:\n");
	for (i = 0; i < n_par; i++) {
		printf("\t%s\t%f\n", get_params_descr(sinmod[0])[i], gsl_vector_get(
				get_steps(sinmod[0]), i) / (gsl_vector_get(
				sinmod[0]->params_max, i) - gsl_vector_get(
				sinmod[0]->params_min, i)));
	}

	printf("Setting startingpoint for the calibration of all hotter "
		"distribution to \n  the best parameter values of the (beta=1)"
		"-distribution\n");
	fflush(stdout);

#pragma omp parallel for
	for (i = 0 + 1; i < n_beta; i++) {
		printf("\tCalibrating chain %d\n", i);
		fflush(stdout);
		gsl_vector_free(sinmod[i]->params_step);
		sinmod[i]->params_step = dup_vector(get_steps(sinmod[0]));
		gsl_vector_scale(sinmod[i]->params_step, pow(get_beta(sinmod[i]), -0.5));
		set_params(sinmod[i], dup_vector(get_params_best(sinmod[0])));
		calc_model(sinmod[i], NULL);

		/*markov_chain_calibrate(sinmod[i], burn_in_iterations, rat_limit,
		 iter_limit, mul, DEFAULT_ADJUST_STEP);*/
	}
	fflush(stdout);
	printf("all chains calibrated.\n");
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, get_beta(sinmod[i]));
		printf("\tsteps: ");
		dump_vectorln(get_steps(sinmod[i]));
	}

	signal(SIGINT, ctrl_c_handler);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);

	analyse(sinmod, n_beta, n_swap);
	for (i = 0; i < n_beta; i++) {
		mem_free(sinmod[i]->additional_data);
		if (i != 0) {
			set_x(sinmod[i], NULL);
			set_y(sinmod[i], NULL);
		}
		sinmod[i] = mcmc_free(sinmod[i]);
		mem_free(sinmod[i]);
	}
	mem_free(sinmod);
}

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Enable Random Walk Metropolis (adaptive MCMC method)
 * Also important: #MINIMAL_STEPWIDTH, #MAXIMAL_STEPWIDTH
 */
#define RWM
#endif

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Enable a constant but small rescaling of the step width to keep the
 * acceptance rate up.
 */
#define ADAPT
#endif

static void adapt(mcmc ** sinmod, const unsigned int n_beta,
		const unsigned int n_swap) {
	unsigned int i;

#ifdef RWM
	static double prob_old[100];
#endif

#ifdef RWM
	for (i = 0; i < n_beta; i++) {
		prob_old[i] = get_prob(sinmod[i]);
		markov_chain_step(sinmod[i], 0);
		rmw_adapt_stepwidth(sinmod[i], prob_old[i]);
	}
#endif
#ifdef ADAPT
	for (i = 0; i < n_beta; i++) {
		if (get_params_accepts_sum(sinmod[i]) + get_params_rejects_sum(
						sinmod[i]) < 20000) {
			continue;
		}
		if (get_params_accepts_sum(sinmod[i]) * 1.0
				/ get_params_rejects_sum(sinmod[i]) < TARGET_ACCEPTANCE_RATE - 0.05) {
			dump_i("too few accepts, scaling down", i);
			gsl_vector_scale(get_steps(sinmod[i]), 0.99);
		} else if (get_params_accepts_sum(sinmod[i]) * 1.0
				/ get_params_rejects_sum(sinmod[i]) > TARGET_ACCEPTANCE_RATE + 0.05) {
			dump_i("too many accepts, scaling up", i);
			gsl_vector_scale(get_steps(sinmod[i]), 1 / 0.99);
		}
		if (get_params_accepts_sum(sinmod[i]) + get_params_rejects_sum(
						sinmod[i]) > 100000) {
			reset_accept_rejects(sinmod[i]);
		}
	}
#endif
	(void) sinmod;
	i = n_beta + n_swap;
}

void dump(const mcmc ** sinmod, const unsigned int n_beta,
		const unsigned long iter, FILE * acceptance_file) {
	unsigned int i;
	if (iter % PRINT_PROB_INTERVAL == 0) {
		if (dumpflag) {
			report(sinmod, n_beta);
			dumpflag = 0;
		}
		fprintf(acceptance_file, "%lu\t", iter);
		for (i = 0; i < n_beta; i++) {
			fprintf(acceptance_file, "%lu\t", get_params_accepts_sum(sinmod[i]));
			fprintf(acceptance_file, "%lu\t", get_params_rejects_sum(sinmod[i]));
		}
		fprintf(acceptance_file, "\n");
		fflush(acceptance_file);
		IFDEBUG {
			debug("dumping distribution");
			dump_ul("iteration", iter);
			dump_ul("acceptance rate: accepts", get_params_accepts_sum(sinmod[0]));
			dump_ul("acceptance rate: rejects", get_params_rejects_sum(sinmod[0]));
			dump_mcmc(sinmod[0]);
		} else {
			printf("iteration: %lu, a/r: %lu/%lu v:", iter,
					get_params_accepts_sum(sinmod[0]), get_params_rejects_sum(
							sinmod[0]));
			dump_vector(get_params(sinmod[0]));
			printf(" [%d/%lu ticks]\r", get_duration(), CLOCKS_PER_SEC);
			fflush(stdout);
		}
	}
}

static void analyse(mcmc ** sinmod, const int n_beta, const unsigned int n_swap) {
	int i;
	unsigned long iter = sinmod[0]->n_iter;
	unsigned int subiter;
	FILE * acceptance_file;

#ifdef DUMP_PROBABILITIES
	FILE ** probabilities_file = (FILE**) mem_calloc(n_beta, sizeof(FILE*));
	char buf[100];
	assert(probabilities_file != NULL);
	for (i = 0; i < n_beta; i++) {
		sprintf(buf, "prob-chain%d.dump", i);
		probabilities_file[i] = fopen(buf, "w");
		if (probabilities_file[i] == NULL) {
			fprintf(stderr, "opening file %s failed\n", buf);
			perror("opening file failed");
			exit(1);
		}
	}
#endif
	assert(n_beta < 100);

	acceptance_file = fopen("acceptance_rate.dump", "w");
	assert(acceptance_file != NULL);
	get_duration();
	run = 1;
	dumpflag = 0;
	printf("starting the analysis\n");
	fflush(stdout);

	while (run
#ifdef MAX_ITERATIONS
	&& iter < MAX_ITERATIONS
#endif
	) {
#pragma omp parallel for
		for (i = 0; i < n_beta; i++) {
			for (subiter = 0; subiter < n_swap; subiter++) {
				markov_chain_step(sinmod[i]);
				mcmc_check_best(sinmod[i]);
				mcmc_append_current_parameters(sinmod[i]);
#ifdef DUMP_PROBABILITIES
				fprintf(probabilities_file[i], "%6e\n", get_prob(sinmod[i]));
#endif
			}
		}
		adapt(sinmod, n_beta, iter);
		iter += n_swap;
		tempering_interaction(sinmod, n_beta, iter);
		dump((const mcmc **)sinmod, n_beta, iter, acceptance_file);
	}
	if (fclose(acceptance_file) != 0) {
		assert(0);
	}
#ifdef DUMP_PROBABILITIES
	for (i = 0; i < n_beta; i++) {
		if (fclose(probabilities_file[i]) != 0) {
			assert(0);
		}
	}
#endif
	report((const mcmc **)sinmod, n_beta);
}

