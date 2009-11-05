#include <omp.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_beta.h"
#include "parallel_tempering_interaction.h"
#include "debug.h"
#include "define_defaults.h"
#include "gsl_helper.h"
#include "parallel_tempering_run.h"
#include "histogram.h"
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

void write_params_file(mcmc * m) {
	unsigned int i;
	FILE * f = fopen(PARAMS_FILENAME "_suggested", "w");
	if (f != NULL) {
		for (i = 0; i < get_n_par(m); i++) {
			fprintf(
					f,
					DUMP_FORMAT "\t" DUMP_FORMAT "\t" DUMP_FORMAT "\t%s\t" DUMP_FORMAT "\n",
					gsl_vector_get(get_params_best(m), i), gsl_vector_get(
							get_params_min(m), i), gsl_vector_get(
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

void write_calibration_summary(mcmc ** chains, unsigned int n_chains) {
	unsigned int i;
	unsigned int j;
	double beta_0 = get_beta(chains[n_chains - 1]);
	unsigned int n_pars = get_n_par(chains[0]);
	FILE * f = fopen("calibration_summary", "w");
	if (f != NULL) {
		fprintf(f, "Summary of calibrations\n");
		fprintf(f, "\nBETA TABLE\n");
		fprintf(f, "Chain # | Calculated | Calibrated\n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "Chain %d | " DUMP_FORMAT " | %f\n", i, get_chain_beta(
					i, n_chains, beta_0), get_beta(chains[i]));
		}
		fprintf(f, "\nSTEPWIDTH TABLE\n");
		fprintf(f, "Chain # | Calibrated stepwidths... \n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "%d", i);
			for (j = 0; j < n_pars; j++) {
				fprintf(f, "\t" DUMP_FORMAT, get_steps_for(chains[i], j));
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\nSTEPWIDTH ESTIMATE TABLE\n");
		fprintf(f, "If you find that the estimate deviates much or "
			"systematically from the");
		fprintf(f, "calibrated stepwidths, please notify the authors.\n");
		fprintf(f, "Chain # | Calculated stepwidths... \n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "%d", i);
			for (j = 0; j < n_pars; j++) {
				fprintf(f, "\t" DUMP_FORMAT, get_steps_for(chains[i], j));
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
	mcmc ** chains;
	const char * params_filename = PARAMS_FILENAME;
	const char * data_filename = DATA_FILENAME;
	const unsigned int n_beta = N_BETA;
	chains = (mcmc**) mem_calloc(n_beta, sizeof(mcmc*));
	assert(chains != NULL);

	printf("Initializing %d chains\n", n_beta);
	for (i = 0; i < n_beta; i++) {
		chains[i] = mcmc_load_params(params_filename);
		if (i == 0) {
			mcmc_load_data(chains[i], data_filename);
			chains[i]->additional_data
					= mem_malloc(sizeof(parallel_tempering_mcmc));
			set_beta(chains[i], 1.0);
		} else {
			mcmc_reuse_data(chains[i], chains[0]);
		}
		mcmc_check(chains[i]);
		chains[i]->additional_data = mem_malloc(
				sizeof(parallel_tempering_mcmc));
		set_beta(chains[i], 1);
	}
	return chains;
}

/**
 * read betas, stepwidths and start values of all chains
 *
 * @return lines read
 **/
void read_calibration_file(mcmc ** chains, unsigned int n_chains) {
	unsigned int i;
	unsigned int j;
	unsigned int err = 0;
	unsigned int n_par = get_n_par(chains[0]);
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
		set_beta(chains[i], v);
		for (j = 0; j < n_par && !feof(f) && err == 0; j++) {
			if (fscanf(f, "%lf", &v) != 1)
				err++;
			set_steps_for(chains[i], v, j);
		}
		for (j = 0; j < n_par && !feof(f) && err == 0; j++) {
			if (fscanf(f, "%lf", &v) != 1)
				err++;
			set_params_for(chains[i], v, j);
		}
		if (feof(f) && j < n_par) {
			fprintf(f, "could not read %d chain calibrations. \nError with "
				"line %d.\n", n_chains, i + 1);
			exit(1);
		}
		set_params_best(chains[i], get_params(chains[i]));
	}

	fclose(f);
}
void write_calibrations_file(mcmc ** chains, const unsigned int n_chains) {
	FILE * f;
	unsigned int i;
	unsigned int j;
	unsigned int n_par = get_n_par(chains[0]);

	f = fopen(CALIBRATION_FILE, "w");
	if (f == NULL) {
		perror("error writing to calibration results file");
		exit(1);
	}
	for (j = 0; j < n_chains; j++) {
		fprintf(f, DUMP_FORMAT, get_beta(chains[j]));
		for (i = 0; i < n_par; i++) {
			fprintf(f, "\t" DUMP_FORMAT, get_steps_for(chains[j], i));
		}
		for (i = 0; i < n_par; i++) {
			fprintf(f, "\t" DUMP_FORMAT, get_params_for(chains[j], i));
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

#ifdef DUMP_ALLCHAINS
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
			printf("iteration: %lu, a/r: %lu/%lu = %f, v:", iter,
					get_params_accepts_global(chains[0]),
					get_params_rejects_global(chains[0]),
					(double) get_params_accepts_global(chains[0])
							/ (double) (get_params_accepts_global(chains[0])
									+ get_params_rejects_global(chains[0])));
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

/*
 * calculate data probability
 */
void analyse_data_probability() {
	unsigned int i;
	unsigned int j;
	unsigned long n = 0;
	double v;
	double w;
	double sums[100];
	double previous_beta;
	double data_logprob;
	char buf[100];
	FILE * f;
	unsigned int n_beta = N_BETA;
	mcmc ** chains = setup_chains();

	read_calibration_file(chains, n_beta);

	assert(n_beta < 100);
	for (i = 0; i < n_beta; i++) {
		sprintf(buf, "prob-chain%d.dump", i);
		dump_s("summing up probability file", buf);
		printf("reading probabilities of chain %d\r", i);
		fflush(stdout);
		f = fopen(buf, "r");
		if (f == NULL) {
			fprintf(stderr,
					"calculating data probability failed: file %s not found\n",
					buf);
			return;
		}
		n = 0;
		sums[i] = 0;
		while (!feof(f)) {
			if (fscanf(f, "%le\t%le", &w, &v) == 2) {
				/*
				 * note: rounding errors could occur here
				 * since the values are of the same magnitude, they hopefully won't
				 */
				sums[i] += v;
				n++;
			}
		}
		if (n == 0) {
			fprintf(stderr, "calculating data probability failed: "
				"no data points found in %s\n", buf);
			return;
		}
		sums[i] = sums[i] / get_beta(chains[i]) / n;
	}

	data_logprob = 0;
	previous_beta = 0;
	/* calculate the integral by an estimate */
	for (j = n_beta - 1;; j--) {
		assert((get_beta(chains[j]) > previous_beta));

		data_logprob += sums[j] * (get_beta(chains[j]) - previous_beta);

		if (j == 0)
			break;
		previous_beta = get_beta(chains[j]);
	}

	printf("Model probability ln(p(D|M, I)): [about 10^%.0f] %.5f"
		"\n"
		"\nTable to compare support against other models (Jeffrey):\n"
		" other model ln(p(D|M,I)) | supporting evidence for this model\n"
		" --------------------------------- \n"
		"        >  %04.1f \tnegative (supports other model)\n"
		"  %04.1f .. %04.1f \tBarely worth mentioning\n"
		"  %04.1f .. %04.1f \tSubstantial\n"
		"  %04.1f .. %04.1f \tStrong\n"
		"  %04.1f .. %04.1f \tVery strong\n"
		"        <  %04.1f \tDecisive\n", data_logprob / gsl_sf_log(10),
			data_logprob, data_logprob, data_logprob, data_logprob - 10,
			data_logprob - 10, data_logprob - 23, data_logprob - 23,
			data_logprob - 34, data_logprob - 34, data_logprob - 46,
			data_logprob - 46);
	printf("\nbe careful.\n");
}

#ifndef NBINS
#define NBINS 200
#endif

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * If not set, the marginal distribution will be calculated for the whole
 * parameter range. Pros: faster, comparable. Cons: Not as detailed.
 *
 * If set, the maximum and minimum values found are used for
 * the histogram. Pros: more detailed in the area of interest
 */
#define HISTOGRAMS_MINMAX
#endif

void calc_marginal_distribution(mcmc ** chains, unsigned int n_beta,
		unsigned int param, int find_minmax) {
	const unsigned int nbins = NBINS;
	char ** filenames;
	unsigned int filecount = n_beta;

	unsigned int i;
	gsl_vector * min;
	gsl_vector * max;
	gsl_histogram * h;
	FILE * outfile;
	char outfilename[100];
	const char * paramname = get_params_descr(chains[0])[param];

#ifdef HISTOGRAMS_ALLCHAINS
	filecount = n_beta;
#else
	filecount = 1;
#endif

	filenames = (char **) calloc(filecount + 1, sizeof(char*));
	for (i = 0; i < filecount; i++) {
		filenames[i] = (char *) malloc(100 * sizeof(char));
		sprintf(filenames[i], "%s-chain-%d.prob.dump", paramname, i);
	}
	sprintf(outfilename, "%s.histogram", paramname);

	min = gsl_vector_alloc(1);
	max = gsl_vector_alloc(1);

	gsl_vector_set(min, 0, get_params_min_for(chains[0], param));
	gsl_vector_set(max, 0, get_params_max_for(chains[0], param));

	if (find_minmax != 0) {
		for (i = 0; i < filecount; i++) {
			printf("minmax search : chain %3d parameter %s   \r", i, paramname);
			fflush(stdout);
			if (1 != get_column_count(filenames[i])) {
				fprintf(
						stderr,
						"number of columns different in file %s: %i vs %i in %s\n",
						filenames[i], 1, get_column_count(filenames[i]),
						filenames[0]);
				exit(1);
			}
			if (i == 0)
				find_min_max(filenames[0], min, max);
			else
				update_min_max(filenames[i], min, max);
		}
		dump_v("minima", min);
		dump_v("maxima", max);
	}
	debug("creating histogram ...");
	h = create_hist(nbins, gsl_vector_get(min, 0), gsl_vector_get(max, 0));
	debug("filling histogram... ");
	for (i = 0; i < filecount; i++) {
		printf("reading values: chain %3d parameter %s   \r", i, paramname);
		fflush(stdout);
		dump_s("with file", filenames[i]);
		append_to_hists(&h, 1, filenames[i]);
		free(filenames[i]);
	}
	free(filenames);
	gsl_histogram_scale(h, (gsl_vector_get(max, 0) - gsl_vector_get(min, 0))
			/ nbins / gsl_histogram_sum(h));

	outfile = fopen(outfilename, "w");
	debug("writing histogram... ");
	assert(outfile != NULL);
	gsl_histogram_fprintf(outfile, h, DUMP_FORMAT, DUMP_FORMAT);
	/*for (i = 0; i < gsl_histogram_bins(h); i++) {
	 gsl_histogram_get_range(h, i, &lower, &upper);

	 fprintf(outfile, DUMP_FORMAT "\t" DUMP_FORMAT "\t" DUMP_FORMAT "\t",
	 lower, upper, gsl_histogram_get(h, i));
	 }*/
	dump_s("histogram file done", outfilename);
	fclose(outfile);

	gsl_histogram_free(h);
}

#ifndef GNUPLOT_STYLE
#define GNUPLOT_STYLE "with histeps"
#endif

void analyse_marginal_distributions() {
	unsigned int i;
	int find_minmax = 0;
	unsigned int n_beta = N_BETA;
	mcmc ** chains = setup_chains();
	FILE * plotplate;

	read_calibration_file(chains, n_beta);

#ifdef HISTOGRAMS_MINMAX
	find_minmax = 1;
#endif

	for (i = 0; i < get_n_par(chains[0]); i++) {
		calc_marginal_distribution(chains, n_beta, i, find_minmax);
	}

	plotplate = fopen("marginal_distributions.gnuplot", "w");
	assert(plotplate != NULL);
	fprintf(plotplate, "# set terminal png size %d,%d; set output "
		"\"marginal_distributions.png\"\n", 600, 300 * get_n_par(chains[0]));
	fprintf(plotplate, "set multiplot\n");
	fprintf(plotplate, "set size 1,%f\n", 1. / get_n_par(chains[0]));
	for (i = 0; i < get_n_par(chains[0]); i++) {
		fprintf(plotplate, "set origin 0,%f\n", (get_n_par(chains[0]) - i - 1)
				* 1. / get_n_par(chains[0]));
		fprintf(plotplate, "plot \"%s.histogram\" u 1:3 title \"%s\" "
		GNUPLOT_STYLE
		"\n", get_params_descr(chains[0])[i], get_params_descr(chains[0])[i]);
	}
	fprintf(plotplate, "unset multiplot\n");
	fclose(plotplate);
}

