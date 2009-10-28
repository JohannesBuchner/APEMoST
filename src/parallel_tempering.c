#include <omp.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_beta.h"
#include "parallel_tempering_interaction.h"
#include "debug.h"
#include "define_defaults.h"
#include "gsl_helper.h"
#include "parallel_tempering_run.h"

void register_signal_handlers();

void report(const mcmc ** sinmod, const int n_beta) {
	int i = 0;
	print_current_positions(sinmod, n_beta);
	printf("\nwriting out visited parameters ");
	while (1) {
		printf(".");
		mcmc_dump_flush(sinmod[i]);
		fflush(stdout);
#ifndef DUMP_ALL_CHAINS
		break;
#endif
		if (i >= n_beta)
			break;
		i++;
	}
	printf("done.\n");
}
void analyse(mcmc ** sinmod, int n_beta, unsigned int n_swap);

void write_params_file(mcmc * m) {
	unsigned int i;
	FILE * f = fopen(PARAMS_FILENAME "_suggested", "w");
	if (f != NULL) {
		for (i = 0; i < get_n_par(m); i++) {
			fprintf(f, "%f\t%f\t%f\t%s\t%f\n", 
					gsl_vector_get(get_params_best(m), i),
					gsl_vector_get(get_params_min(m), i),
					gsl_vector_get(get_params_max(m), i),
					get_params_descr(m)[i], 
					gsl_vector_get(get_steps(m), i));
		}
		fclose(f);
		printf("new suggested parameters file has been written\n");
	} else {
		fprintf(stderr, "Could not write to file " PARAMS_FILENAME "_suggested\n");
	}
}

void write_calibration_summary(mcmc ** sinmod, unsigned int n_chains) {
	unsigned int i;
	unsigned int j;
	double beta_0 = get_beta(sinmod[0]);
	unsigned int n_pars = get_n_par(sinmod[0]);
	FILE * f = fopen("calibration_summary", "w");
	if (f != NULL) {
		fprintf(f, "Summary of calibrations\n");
		fprintf(f, "\nBETA TABLE\n");
		fprintf(f, "Chain # | Calculated | Calibrated\n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "Chain %d | %f | %f\n", i, 
					get_chain_beta(i, n_chains, beta_0), get_beta(sinmod[i]));
		}
		fprintf(f, "\nSTEPWIDTH TABLE\n");
		fprintf(f, "Chain # | Calibrated stepwidths... \n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "%d", i);
			for (j = 0; j < n_pars; j++) {
				fprintf(f, "\t%f", get_steps_for(sinmod[i], j));
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\nSTEPWIDTH ESTIMATE TABLE\n");
		fprintf(f, "If you find that the estimate deviates much or "
			"systematically from the");
		fprintf(f, "calibrated stepwidths, please contact the authors.\n");
		fprintf(f, "Chain # | Calculated stepwidths... \n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "%d", i);
			for (j = 0; j < n_pars; j++) {
				fprintf(f, "\t%f", get_steps_for(sinmod[i], j));
			}
			fprintf(f, "\n");
		}
		fclose(f);
		printf("calibration summary has been written\n");
	} else {
		fprintf(stderr, "Could not write to file calibration_summary\n");
	}
}
/*
	parallel_tempering(PARAMS_FILENAME, DATA_FILENAME, N_BETA, BETA_0,
			BURN_IN_ITERATIONS, RAT_LIMIT, ITER_LIMIT, MUL, N_SWAP);
*/
void parallel_tempering(const char * params_filename,
		const char * data_filename, const int n_beta, double beta_0,
		const unsigned long burn_in_iterations, const double rat_limit,
		const unsigned long iter_limit, const double mul, int n_swap) {
	int n_par;
	int i;
	const char ** params_descr;
	mcmc ** sinmod;
	gsl_vector * stepwidth_factors;

	sinmod = (mcmc**) mem_calloc(n_beta, sizeof(mcmc*));
	assert(sinmod != NULL);

	if (n_swap < 0) {
		n_swap = 2000 / n_beta;
		printf("automatic n_swap: %d\n", n_swap);
	}

	printf("Initializing parallel tempering for %d chains\n", n_beta);
	for (i = 0; i < n_beta; i++) {
		sinmod[i] = mcmc_load_params(params_filename);
		if (i == 0) {
			mcmc_load_data(sinmod[i], data_filename);
			sinmod[i]->additional_data
					= mem_malloc(sizeof(parallel_tempering_mcmc));
			set_beta(sinmod[i], 1.0);
		} else {
			mcmc_reuse_data(sinmod[i], sinmod[0]);
		}
		mcmc_check(sinmod[i]);
	}
	params_descr = get_params_descr(sinmod[0]);
	n_par = get_n_par(sinmod[0]);
	stepwidth_factors = gsl_vector_alloc(n_par);
	gsl_vector_set_all(stepwidth_factors, 1);

	printf("Starting markov chain calibration\n");
	fflush(stdout);
	calc_model(sinmod[0], NULL);
	mcmc_check(sinmod[0]);
	markov_chain_calibrate(sinmod[0], burn_in_iterations, rat_limit,
			iter_limit, mul, DEFAULT_ADJUST_STEP);
	mcmc_open_dump_files(sinmod[0], "-chain", 0);
	
	write_params_file(sinmod[0]);

	printf("Calibrating chains\n");
	fflush(stdout);

	i = 1;
	if (n_beta >= 2) {
		sinmod[i]->additional_data = mem_malloc(
				sizeof(parallel_tempering_mcmc));
		if (beta_0 < 0)
			set_beta(sinmod[i], get_chain_beta(i, n_beta, calc_beta_0(
					sinmod[0], stepwidth_factors)));
		else
			set_beta(sinmod[i], get_chain_beta(i, n_beta, beta_0));
		gsl_vector_free(get_steps(sinmod[i]));
		sinmod[i]->params_step = dup_vector(get_steps(sinmod[0]));
		gsl_vector_scale(get_steps(sinmod[i]), pow(get_beta(sinmod[i]), -0.5));
		set_params(sinmod[i], dup_vector(get_params_best(sinmod[0])));
		calc_model(sinmod[i], NULL);
		printf("Calibrating second chain to infer stepwidth factor\n");
		printf("\tChain %2d - ", i);
		printf("beta = %f\tsteps: ", get_beta(sinmod[i]));
		dump_vectorln(get_steps(sinmod[i]));
		fflush(stdout);
		markov_chain_calibrate(sinmod[i], burn_in_iterations, rat_limit,
				iter_limit, mul, DEFAULT_ADJUST_STEP);
		gsl_vector_scale(stepwidth_factors, pow(get_beta(sinmod[i]), -0.5));
		gsl_vector_mul(stepwidth_factors, get_steps(sinmod[0]));
		gsl_vector_div(stepwidth_factors, get_steps(sinmod[i]));
		mem_free(sinmod[i]->additional_data);
	}

	printf("stepwidth factors: ");
	dump_vectorln(stepwidth_factors);

	if (beta_0 < 0) {
		beta_0 = calc_beta_0(sinmod[0], stepwidth_factors);
		printf("automatic beta_0: %f\n", beta_0);
	}

	fflush(stdout);

#pragma omp parallel for
	for (i = 1; i < n_beta; i++) {
		printf("\tChain %2d - ", i);
		fflush(stdout);
		sinmod[i]->additional_data
				= mem_malloc(sizeof(parallel_tempering_mcmc));
		set_beta(sinmod[i], get_chain_beta(i, n_beta, beta_0));
		gsl_vector_free(get_steps(sinmod[i]));
		sinmod[i]->params_step = dup_vector(get_steps(sinmod[0]));
		gsl_vector_scale(get_steps(sinmod[i]), pow(get_beta(sinmod[i]), -0.5));
		gsl_vector_mul(get_steps(sinmod[i]), stepwidth_factors);
		set_params(sinmod[i], dup_vector(get_params_best(sinmod[0])));
		calc_model(sinmod[i], NULL);
		printf("beta = %f\tsteps: ", get_beta(sinmod[i]));
		dump_vectorln(get_steps(sinmod[i]));
		fflush(stdout);
#ifndef SKIP_CALIBRATE_ALLCHAINS
		markov_chain_calibrate(sinmod[i], burn_in_iterations, rat_limit,
				iter_limit, mul, DEFAULT_ADJUST_STEP);
#endif
#ifdef DUMP_ALL_CHAINS
		mcmc_open_dump_files(sinmod[i], "-chain", i);
#endif
	}
	gsl_vector_free(stepwidth_factors);
	fflush(stdout);
	printf("all chains calibrated.\n");
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f \tsteps: ", i, get_beta(sinmod[i]));
		dump_vectorln(get_steps(sinmod[i]));
	}
	write_calibration_summary(sinmod, n_beta);

	register_signal_handlers();

	analyse(sinmod, n_beta, n_swap);
	for (i = 0; i < n_beta; i++) {
		mem_free(sinmod[i]->additional_data);
		if (i != 0) {
			/* this was reused, thus avoid double free */
			set_data(sinmod[i], NULL);
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

void adapt(mcmc ** sinmod, const unsigned int n_beta,
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
		   printf("iteration: %lu, a/r: %lu/%lu = %f, v:", iter,
			get_params_accepts_sum(sinmod[0]), get_params_rejects_sum(
			sinmod[0]), (double) get_params_accepts_sum(sinmod[0])/
			(double) (get_params_accepts_sum(sinmod[0]) + get_params_rejects_sum(sinmod[0])));
			dump_vector(get_params(sinmod[0]));
			printf(" [%d/%lu ticks]\r", get_duration(), get_ticks_per_second());
			fflush(stdout);
		}
	}
}

void analyse(mcmc ** sinmod, const int n_beta, const unsigned int n_swap) {
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
#ifdef DUMP_ALL_CHAINS
				if (i == 0)
#endif
				mcmc_append_current_parameters(sinmod[i]);
#ifdef DUMP_PROBABILITIES
				fprintf(probabilities_file[i], "%6e\n", get_prob(sinmod[i]));
#endif
			}
		}
		adapt(sinmod, n_beta, iter);
		iter += n_swap;
		tempering_interaction(sinmod, n_beta, iter);
		dump((const mcmc **) sinmod, n_beta, iter, acceptance_file);
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
	report((const mcmc **) sinmod, n_beta);
}

