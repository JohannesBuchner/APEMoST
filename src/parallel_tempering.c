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

void run_sampler(mcmc ** sinmod, int n_beta, unsigned int n_swap);

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

void write_params_file(mcmc * m) {
	unsigned int i;
	FILE * f = fopen(PARAMS_FILENAME "_suggested", "w");
	if (f != NULL) {
		for (i = 0; i < get_n_par(m); i++) {
			fprintf(f, "%f\t%f\t%f\t%s\t%f\n", gsl_vector_get(
					get_params_best(m), i),
					gsl_vector_get(get_params_min(m), i), gsl_vector_get(
							get_params_max(m), i), get_params_descr(m)[i],
					gsl_vector_get(get_steps(m), i));
		}
		fclose(f);
		printf("new suggested parameters file has been written\n");
	} else {
		fprintf(stderr,
				"Could not write to file " PARAMS_FILENAME "_suggested\n");
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
			fprintf(f, "Chain %d | %f | %f\n", i, get_chain_beta(i, n_chains,
					beta_0), get_beta(sinmod[i]));
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
mcmc ** setup_chains() {
	unsigned int i;
	mcmc ** sinmod;
	const char * params_filename = PARAMS_FILENAME;
	const char * data_filename = DATA_FILENAME;
	const unsigned int n_beta = N_BETA;
	sinmod = (mcmc**) mem_calloc(n_beta, sizeof(mcmc*));
	assert(sinmod != NULL);

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
		sinmod[i]->additional_data = mem_malloc(
				sizeof(parallel_tempering_mcmc));
		set_beta(sinmod[i], 1);
	}
	return sinmod;
}

#define CALIBRATION_FILE "calibration_results"

/**
 * read betas, stepwidths and start values of all chains
 *
 * @return lines read
 **/
void read_calibration_file(mcmc ** sinmod, unsigned int n_chains) {
	unsigned int i;
	unsigned int j;
	unsigned int err = 0;
	unsigned int n_par = get_n_par(sinmod[0]);
	double v;
	FILE * f;

	f = fopen(CALIBRATION_FILE, "r");
	if (f == NULL) {
		fprintf(f, "could not read calibration file %s\n", CALIBRATION_FILE);
		exit(1);
	}

	for (i = 0; i < n_chains && !feof(f); i++) {
		if (fscanf(f, "%lf", &v) != 1)
			err++;
		set_beta(sinmod[i], v);
		for (j = 0; j < n_par && !feof(f) && err == 0; j++) {
			if (fscanf(f, "%lf", &v) != 1)
				err++;
			set_steps_for(sinmod[i], v, j);
		}
		for (j = 0; j < n_par && !feof(f) && err == 0; j++) {
			if (fscanf(f, "%lf", &v) != 1)
				err++;
			set_params_for(sinmod[i], v, j);
		}
		if (feof(f) && j < n_par) {
			fprintf(f, "could not read %d chain calibrations. \nError with "
				"line %d.\n", n_chains, i + 1);
		}
	}

	fclose(f);
}
void write_calibrations_file(mcmc ** sinmod, const unsigned int n_chains) {
	FILE * f;
	unsigned int i;
	unsigned int j;
	unsigned int n_par = get_n_par(sinmod[0]);

	f = fopen(CALIBRATION_FILE, "w");
	if (f == NULL) {
		perror("error writing to calibration results file");
		exit(1);
	}
	for (j = 0; j < n_chains; j++) {
		fprintf(f, "%f", get_beta(sinmod[j]));
		for (i = 0; i < n_par; i++) {
			/* FIXME: precision */
			fprintf(f, "\t%f", get_steps_for(sinmod[j], i));
		}
		for (i = 0; i < n_par; i++) {
			/* FIXME: precision */
			fprintf(f, "\t%f", get_params_best_for(sinmod[j], i));
		}
		fprintf(f, "\n");
	}
	fclose(f);
	printf("wrote calibration results for %d chains to %s\n", n_chains,
			CALIBRATION_FILE);

}

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
	mcmc ** sinmod = setup_chains();
	const double rat_limit = RAT_LIMIT;
	const unsigned long burn_in_iterations = BURN_IN_ITERATIONS;
	const unsigned long iter_limit = ITER_LIMIT;
	const double mul = MUL;

	printf("Starting markov chain calibration\n");
	fflush(stdout);
	calc_model(sinmod[0], NULL);
	mcmc_check(sinmod[0]);
	markov_chain_calibrate(sinmod[0], burn_in_iterations, rat_limit,
			iter_limit, mul, DEFAULT_ADJUST_STEP);
	write_calibrations_file(sinmod, 1);
	write_params_file(sinmod[0]);
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
	const double rat_limit = RAT_LIMIT;
	double beta_0 = BETA_0;
	const unsigned long burn_in_iterations = BURN_IN_ITERATIONS;
	const unsigned long iter_limit = ITER_LIMIT;
	const double mul = MUL;
	unsigned int n_par;
	int i;
	gsl_vector * stepwidth_factors;
	mcmc ** sinmod = setup_chains();

	read_calibration_file(sinmod, 1);

	printf("Calibrating chains\n");
	fflush(stdout);
	n_par = get_n_par(sinmod[0]);
	stepwidth_factors = gsl_vector_alloc(n_par);
	gsl_vector_set_all(stepwidth_factors, 1);

	i = 1;
	if (n_beta > 1) {
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
	}
	gsl_vector_free(stepwidth_factors);
	fflush(stdout);
	printf("all chains calibrated.\n");
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f \tsteps: ", i, get_beta(sinmod[i]));
		dump_vectorln(get_steps(sinmod[i]));
	}
	write_calibration_summary(sinmod, n_beta);
	write_calibrations_file(sinmod, n_beta);

	register_signal_handlers();
}

void prepare_and_run_sampler() {
	unsigned int n_beta = N_BETA;
	unsigned int i;
	int n_swap = N_SWAP;

	mcmc ** sinmod = setup_chains();

	read_calibration_file(sinmod, n_beta);

	for (i = 0; i < n_beta; i++) {
		mcmc_open_dump_files(sinmod[i], "-chain", i);
	}

	if (n_swap < 0) {
		n_swap = 2000 / n_beta;
		printf("automatic n_swap: %d\n", n_swap);
	}

	run_sampler(sinmod, n_beta, n_swap);

	report((const mcmc **) sinmod, n_beta);

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

void adapt(mcmc ** sinmod, const unsigned int n_beta, const unsigned int n_swap) {
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
	/* avoid unused warnings if adapt is disabled */
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
							sinmod[0]), (double) get_params_accepts_sum(
							sinmod[0]) / (double) (get_params_accepts_sum(
							sinmod[0]) + get_params_rejects_sum(sinmod[0])));
			dump_vector(get_params(sinmod[0]));
			printf(" [%d/%lu ticks]\r", get_duration(), get_ticks_per_second());
			fflush(stdout);
		}
	}
}

void run_sampler(mcmc ** sinmod, const int n_beta, const unsigned int n_swap) {
	int i;
	unsigned long iter = sinmod[0]->n_iter;
	unsigned int subiter;
	FILE * acceptance_file;

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
				fprintf(probabilities_file[i], "%6e\n", get_prob(sinmod[i]));
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
	for (i = 0; i < n_beta; i++) {
		if (fclose(probabilities_file[i]) != 0) {
			assert(0);
		}
	}
}

