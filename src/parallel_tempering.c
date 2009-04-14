#include <signal.h>
#include <gsl/gsl_sf.h>
#include <omp.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include "parallel_tempering.h"
#include "tempering_interaction.h"
#include "debug.h"

int run;
int dumpflag;

void set_beta(mcmc * m, double newbeta) {
	((parallel_tempering_mcmc *) m->additional_data)->beta = newbeta;
	((parallel_tempering_mcmc *) m->additional_data)->swapcount = 0;
}
double get_beta(mcmc * m) {
	return ((parallel_tempering_mcmc *) m->additional_data)->beta;
}
void inc_swapcount(mcmc * m) {
	((parallel_tempering_mcmc *) m->additional_data)->swapcount++;
}
unsigned long get_swapcount(mcmc * m) {
	return ((parallel_tempering_mcmc *) m->additional_data)->swapcount;
}

void ctrl_c_handler(int signalnr) {
	printf("\nreceived Ctrl-C (%d). Stopping ... (please be patient)\n\n",
			signalnr);
	run = 0;
}
void sigusr_handler(int signalnr) {
	printf("\nreceived SIGUSR (%d). Will dump at next opportunity.\n\n",
			signalnr);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);
	dumpflag = 1;
}

#include <time.h>

int get_duration() {
	static clock_t stored = 0;
	clock_t new = stored;
	stored = clock();
	return stored - new;
}

void print_current_positions(mcmc ** sinmod, int n_beta) {
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

void report(mcmc ** sinmod, int n_beta) {
	print_current_positions(sinmod, n_beta);
	mcmc_dump_probabilities(sinmod[0], -1);
}

int run;
int dumpflag;

void analyse(mcmc ** sinmod, int n_beta, int n_swap);

double equidistant_beta(unsigned int i, unsigned int n_beta, double beta_0) {
	return 1 - i * (1 - beta_0) / (n_beta - 1);
}
double equidistant_temperature(unsigned int i, unsigned int n_beta,
		double beta_0) {
	return 1 / (1 / beta_0 + i * (1 - 1 / beta_0) / (n_beta - 1));
}
double chebyshev_temperature(unsigned int i, unsigned int n_beta, double beta_0) {
	return 1 / (1 / beta_0 + (1 - 1 / beta_0) / 2 * (1 - cos(i * M_PI / (n_beta
			- 1))));
}
double chebyshev_beta(unsigned int i, unsigned int n_beta, double beta_0) {
	return beta_0 + (1 - beta_0) / 2 * (1 - cos(i * M_PI / (n_beta - 1)));
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
 * chebyshev_temperature
 *
 */
#define BETA_DISTRIBUTION
#endif

double get_chain_beta(unsigned int i, unsigned int n_beta, double beta_0) {
#ifndef BETA_DISTRIBUTION
#error "choose a beta distribution, e.g.: BETA_DISTRIBUTION=equidistant_beta"
#error "also available: equidistant_temperature, chebyshev_beta, chebyshev_temperature"
#endif
	if (n_beta == 1)
		return 1.0;
	return BETA_DISTRIBUTION(i, n_beta, beta_0);
}

void parallel_tempering(const char * params_filename,
		const char * data_filename, int n_beta, double beta_0,
		unsigned long burn_in_iterations, double rat_limit,
		unsigned long iter_limit, double mul, int n_swap) {
	int n_par;
	int i;
	const char ** params_descr;
	mcmc ** sinmod;

	sinmod = (mcmc**) calloc(n_beta, sizeof(mcmc*));
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
		sinmod[i]->additional_data = malloc(sizeof(parallel_tempering_mcmc));
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
	for (i = 0; i < n_beta; i++) {
		calc_model(sinmod[i], NULL);
		mcmc_check(sinmod[i]);
	}
	printf("Starting markov chain calibration\n");
	wait();
	markov_chain_calibrate(sinmod[0], burn_in_iterations, rat_limit,
			iter_limit, mul, DEFAULT_ADJUST_STEP);

	printf("Setting startingpoint for the calibration of all hotter "
		"distribution to \n  the best parameter values of the (beta=1)"
		"-distribution\n");
	wait();

#pragma omp parallel for
	for (i = 0 + 1; i < n_beta; i++) {
		printf("\tCalibrating chain %d\n", i);
		set_params(sinmod[i], dup_vector(get_params_best(sinmod[0])));
		calc_model(sinmod[i], NULL);

		markov_chain_calibrate(sinmod[i], burn_in_iterations, rat_limit,
				iter_limit, mul, DEFAULT_ADJUST_STEP);
	}
	printf("all chains calibrated.\n");
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, get_beta(sinmod[i]));
		printf("\tsteps: ");
		dump_vectorln(get_steps(sinmod[i]));
	}
	wait();

	signal(SIGINT, ctrl_c_handler);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);

	analyse(sinmod, n_beta, n_swap);
	for (i = 0; i < n_beta; i++) {
		sinmod[i] = mcmc_free(sinmod[i]);
		free(sinmod[i]);
	}
	free(sinmod);
}

void analyse(mcmc ** sinmod, int n_beta, int n_swap) {
	int i;
	unsigned long iter = sinmod[0]->n_iter;
	int subiter;
	get_duration();
	run = 1;
	dumpflag = 0;
	printf("starting the analysis\n");
	wait();

	while (run
#ifdef MAX_ITERATIONS
	&& iter < MAX_ITERATIONS
#endif
	) {
#pragma omp parallel for
		for (i = 0; i < n_beta; i++) {
			for (subiter = 0; subiter < n_swap; subiter++) {
				markov_chain_step(sinmod[i], 0);
				mcmc_check_best(sinmod[i]);
				if (i == 0) {
					mcmc_append_current_parameters(sinmod[0]);
				}
			}
		}
		iter += n_swap;
		tempering_interaction(sinmod, n_beta, iter);
		if (iter % PRINT_PROB_INTERVAL == 0) {
			if (dumpflag) {
				mcmc_dump_probabilities(sinmod[0], DUMP_PROB_LENGTH);
				print_current_positions(sinmod, n_beta);
				dumpflag = 0;
			}
			IFDEBUG {
				debug("dumping distribution");
				dump_ul("iteration", iter);
				dump_ul("acceptance rate: accepts", get_params_accepts_sum(sinmod[0]));
				dump_ul("acceptance rate: rejects", get_params_rejects_sum(sinmod[0]));
				dump(sinmod[0]);
			} else {
				printf("iteration: %lu, a/r: %lu/%lu v:", iter,
						get_params_accepts_sum(sinmod[0]),
						get_params_rejects_sum(sinmod[0]));
				dump_vector(get_params(sinmod[0]));
				printf(" [%d/%lu ticks]\r", get_duration(), CLOCKS_PER_SEC);
				fflush(stdout);
			}
		}
	}
	report(sinmod, n_beta);
}

