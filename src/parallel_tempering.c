/*
    APEMoST - Automated Parameter Estimation and Model Selection Toolkit
    Copyright (C) 2009  Johannes Buchner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <omp.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_beta.h"
#include "parallel_tempering_interaction.h"
#include "parallel_tempering_config.h"
#include "debug.h"
#include "define_defaults.h"
#include "gsl_helper.h"
#include "parallel_tempering_run.h"
#include "utils.h"

void register_signal_handlers();

void run_sampler(mcmc ** chains, int n_beta, unsigned int n_swap, char * mode);

void report(const mcmc ** chains, const int n_beta) {
	int i = 0;
	print_current_positions(chains, n_beta);
	printf("\nwriting out visited parameters ");
	while (1) {
		printf(".");
		mcmc_dump_flush(chains[i]);
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


#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * how many iterations should be spent on burn-in
 */
#define BURN_IN_ITERATIONS
#endif

/**
 * needs:
 * params file
 * BURN_IN_ITERATIONS
 * start values (in params file)
 **
 * does:
 * calibrate first chain (beta = 1)
 * writes beta, stepwidths and start values as first line in file
 *   calibration_result
 **
 * provides:
 * stepwidths of first chain (calibration_result)
 * new params file (params_suggest)
 * new start values (calibration_result)
 **/
void calibrate_first() {
	mcmc ** chains = setup_chains();
	const double desired_acceptance_rate = TARGET_ACCEPTANCE_RATE;
	const double max_ar_deviation = MAX_AR_DEVIATION;
	const unsigned long burn_in_iterations = BURN_IN_ITERATIONS;
	const unsigned long iter_limit = ITER_LIMIT;
	const double mul = MUL;

	printf("Starting markov chain calibration\n");
	fflush(stdout);
	calc_model(chains[0], NULL);
	mcmc_check(chains[0]);
	markov_chain_calibrate(chains[0], burn_in_iterations,
			desired_acceptance_rate, max_ar_deviation, iter_limit, mul,
			DEFAULT_ADJUST_STEP);
	write_calibrations_file(chains, 1);
	write_params_file(chains[0]);
}

/**
 * needs:
 * params file
 * BURN_IN_ITERATIONS
 * first line in calibration_result
 * BETA_ALIGNMENT
 * BETA_0
 * SKIP_CALIBRATE_ALLCHAINS
 **
 * does:
 * calibrate remaining chains (beta < 1)
 * writes all betas, stepwidths and start values in file calibration_result
 **
 * provides:
 * stepwidths of first chain (calibration_result)
 * new params file (params_suggest)
 * new start values (calibration_result)
 **/
void calibrate_rest() {
	int n_beta = N_BETA;
	const double desired_acceptance_rate = TARGET_ACCEPTANCE_RATE;
	const double max_ar_deviation = MAX_AR_DEVIATION;
	double beta_0 = BETA_0;
	const unsigned long burn_in_iterations = BURN_IN_ITERATIONS;
	const unsigned long iter_limit = ITER_LIMIT;
	const double mul = MUL;
	unsigned int n_par;
	int i;
	gsl_vector * stepwidth_factors;
	mcmc ** chains = setup_chains();

	read_calibration_file(chains, 1);

	printf("Calibrating chains\n");
	fflush(stdout);
	n_par = get_n_par(chains[0]);
	stepwidth_factors = gsl_vector_alloc(n_par);
	gsl_vector_set_all(stepwidth_factors, 1);

	i = 1;
	if (n_beta > 1) {
		if (beta_0 < 0)
			set_beta(chains[i], get_chain_beta(i, n_beta, calc_beta_0(
					chains[0], stepwidth_factors)));
		else
			set_beta(chains[i], get_chain_beta(i, n_beta, beta_0));
		gsl_vector_free(get_steps(chains[i]));
		chains[i]->params_step = dup_vector(get_steps(chains[0]));
		gsl_vector_scale(get_steps(chains[i]), pow(get_beta(chains[i]), -0.5));
		set_params(chains[i], dup_vector(get_params_best(chains[0])));
		calc_model(chains[i], NULL);
		mcmc_check(chains[i]);
		printf("Calibrating second chain to infer stepwidth factor\n");
		printf("\tChain %2d - ", i);
		printf("beta = %f\tsteps: ", get_beta(chains[i]));
		dump_vectorln(get_steps(chains[i]));
		fflush(stdout);
		markov_chain_calibrate(chains[i], burn_in_iterations,
				desired_acceptance_rate, max_ar_deviation, iter_limit, mul,
				DEFAULT_ADJUST_STEP);
		gsl_vector_scale(stepwidth_factors, pow(get_beta(chains[i]), -0.5));
		gsl_vector_mul(stepwidth_factors, get_steps(chains[0]));
		gsl_vector_div(stepwidth_factors, get_steps(chains[i]));
		mem_free(chains[i]->additional_data);
	}

	printf("stepwidth factors: ");
	dump_vectorln(stepwidth_factors);

	if (beta_0 < 0) {
		beta_0 = calc_beta_0(chains[0], stepwidth_factors);
		printf("automatic beta_0: %f\n", beta_0);
	}

	fflush(stdout);

#pragma omp parallel for
	for (i = 1; i < n_beta; i++) {
		printf("\tChain %2d - ", i);
		fflush(stdout);
		chains[i]->additional_data
				= mem_malloc(sizeof(parallel_tempering_mcmc));
		set_beta(chains[i], get_chain_beta(i, n_beta, beta_0));
		gsl_vector_free(get_steps(chains[i]));
		chains[i]->params_step = dup_vector(get_steps(chains[0]));
		gsl_vector_scale(get_steps(chains[i]), pow(get_beta(chains[i]), -0.5));
		gsl_vector_mul(get_steps(chains[i]), stepwidth_factors);
		set_params(chains[i], dup_vector(get_params_best(chains[0])));
		calc_model(chains[i], NULL);
		mcmc_check(chains[i]);
		printf("beta = %f\tsteps: ", get_beta(chains[i]));
		dump_vectorln(get_steps(chains[i]));
		fflush(stdout);
#ifndef SKIP_CALIBRATE_ALLCHAINS
		markov_chain_calibrate(chains[i], burn_in_iterations,
				desired_acceptance_rate, max_ar_deviation, iter_limit, mul,
				DEFAULT_ADJUST_STEP);
#else
		burn_in(chains[i], burn_in_iterations);
#endif
	}
	gsl_vector_free(stepwidth_factors);
	fflush(stdout);
	printf("all chains calibrated.\n");
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f \tsteps: ", i, get_beta(chains[i]));
		dump_vectorln(get_steps(chains[i]));
	}
	write_calibration_summary(chains, n_beta);
	write_calibrations_file(chains, n_beta);
}

void prepare_and_run_sampler(int append) {
	unsigned int n_beta = N_BETA;
	unsigned int i = 0;
	int n_swap = N_SWAP;

	mcmc ** chains = setup_chains();

	read_calibration_file(chains, n_beta);

	mcmc_open_dump_files(chains[i], "-chain", i, (append == 1 ? "a" : "w"));

#ifdef DUMP_ALL_CHAINS
	for (i = 1; i < n_beta; i++) {
		mcmc_open_dump_files(chains[i], "-chain", i, (append == 1 ? "a" : "w"));
	}
#endif

	if (n_swap < 0) {
		n_swap = 2000 / n_beta;
		printf("automatic n_swap: %d\n", n_swap);
	}

	register_signal_handlers();
	run_sampler(chains, n_beta, n_swap, (append == 1 ? "a" : "w"));

	report((const mcmc **) chains, n_beta);

	for (i = 0; i < n_beta; i++) {
		mem_free(chains[i]->additional_data);
		if (i != 0) {
			/* this was reused, thus avoid double free */
			set_data(chains[i], NULL);
		}
		chains[i] = mcmc_free(chains[i]);
		mem_free(chains[i]);
	}
	mem_free(chains);
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

void adapt(mcmc ** chains, const unsigned int n_beta, const unsigned int n_swap) {
	unsigned int i;

#ifdef RWM
	static double prob_old[100];
#endif

#ifdef RWM
	for (i = 0; i < n_beta; i++) {
		prob_old[i] = get_prob(chains[i]);
		markov_chain_step(chains[i], 0);
		rmw_adapt_stepwidth(chains[i], prob_old[i]);
	}
#endif
#ifdef ADAPT
	for (i = 0; i < n_beta; i++) {
		if (get_params_accepts_sum(chains[i]) + get_params_rejects_sum(
						chains[i]) < 20000) {
			continue;
		}
		if (get_params_accepts_sum(chains[i]) * 1.0
				/ get_params_rejects_sum(chains[i]) < TARGET_ACCEPTANCE_RATE - 0.05) {
			dump_i("too few accepts, scaling down", i);
			gsl_vector_scale(get_steps(chains[i]), 0.99);
		} else if (get_params_accepts_sum(chains[i]) * 1.0
				/ get_params_rejects_sum(chains[i]) > TARGET_ACCEPTANCE_RATE + 0.05) {
			dump_i("too many accepts, scaling up", i);
			gsl_vector_scale(get_steps(chains[i]), 1 / 0.99);
		}
		if (get_params_accepts_sum(chains[i]) + get_params_rejects_sum(
						chains[i]) > 100000) {
			reset_accept_rejects(chains[i]);
		}
	}
#endif
	/* avoid unused warnings if adapt is disabled */
	(void) chains;
	i = n_beta + n_swap;
}

void dump(const mcmc ** chains, const unsigned int n_beta,
		const unsigned long iter, FILE * acceptance_file,
		FILE ** probabilities_file) {
	unsigned int i;
	if (iter % PRINT_PROB_INTERVAL == 0) {
		if (dumpflag) {
			report(chains, n_beta);
			dumpflag = 0;
			for (i = 0; i < n_beta; i++) {
				fflush(probabilities_file[i]);
			}
		}
		fprintf(acceptance_file, "%lu", iter);
		for (i = 0; i < n_beta; i++) {
			fprintf(acceptance_file, "\t%lu", get_params_accepts_global(
					chains[i]));
		}
		fprintf(acceptance_file, "\n");
		fflush(acceptance_file);
		IFDEBUG {
			debug("dumping distribution");
			dump_ul("iteration", iter);
			dump_ul("acceptance rate: accepts", get_params_accepts_global(chains[0]));
			dump_ul("acceptance rate: rejects", get_params_rejects_global(chains[0]));
			dump_mcmc(chains[0]);
		} else {
			printf("iteration: %lu, a/r: %.3f(%lu/%lu), v:", iter,
					(double) get_params_accepts_global(chains[0])
							/ (double) (get_params_accepts_global(chains[0])
									+ get_params_rejects_global(chains[0])),
					get_params_accepts_global(chains[0]),
					get_params_rejects_global(chains[0]));
			dump_vector(get_params(chains[0]));
			printf(" [%d/%lu ticks]\r", get_duration(), get_ticks_per_second());
			fflush(stdout);
		}
	}
}

void run_sampler(mcmc ** chains, const int n_beta, const unsigned int n_swap,
		char * mode) {
	int i;
	unsigned long iter = chains[0]->n_iter;
	unsigned int subiter;
	FILE * acceptance_file;

	FILE ** probabilities_file = (FILE**) mem_calloc(n_beta, sizeof(FILE*));
	char buf[100];
	assert(probabilities_file != NULL);
	for (i = 0; i < n_beta; i++) {
		sprintf(buf, "prob-chain%d.dump", i);
		probabilities_file[i] = fopen(buf, mode);
		if (probabilities_file[i] == NULL) {
			fprintf(stderr, "opening file %s failed\n", buf);
			perror("opening file failed");
			exit(1);
		}
	}
	assert(n_beta < 100);

	acceptance_file = fopen("acceptance_rate.dump.gnuplot", "w");
	if (acceptance_file != NULL) {
		fprintf(acceptance_file,
				"# format: iteration | number of accepts for each chain\n");
		fprintf(acceptance_file, "plot ");
		for (i = 0; i < n_beta; i++) {
			fprintf(
					acceptance_file,
					"\"acceptance_rate.dump\" u 1:%d title \"chain %d, beta = %f\"",
					i + 2, i, get_beta(chains[i]));
			if (i != n_beta - 1)
				fprintf(acceptance_file, ", ");
		}
		fprintf(acceptance_file, "\n");
		fclose(acceptance_file);
	}
	acceptance_file = fopen("acceptance_rate.dump", mode);
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
				markov_chain_step(chains[i]);
				mcmc_check_best(chains[i]);
				mcmc_append_current_parameters(chains[i]);
				fprintf(probabilities_file[i], "%6e\t%6e\n",
						get_prob(chains[i]), get_prob(chains[i]) - get_prior(
								chains[i]));
			}
		}
		adapt(chains, n_beta, iter);
		iter += n_swap;
		tempering_interaction(chains, n_beta, iter);
		dump((const mcmc **) chains, n_beta, iter, acceptance_file,
				probabilities_file);
	}
	if (fclose(acceptance_file) != 0) {
		assert(0);
	}
	for (i = 0; i < n_beta; i++) {
		if (fclose(probabilities_file[i]) != 0) {
			assert(0);
		}
	}
	printf("handled %lu iterations on %d chains\n", iter, n_beta);
}


