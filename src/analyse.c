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
#include "parallel_tempering_config.h"
#include "debug.h"
#include "define_defaults.h"
#include "gsl_helper.h"
#include "parallel_tempering_run.h"
#include "histogram.h"
#include "utils.h"

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
			data_logprob, data_logprob, data_logprob, 
			data_logprob - gsl_sf_log(3), data_logprob - gsl_sf_log(3), 
			data_logprob - gsl_sf_log(10), data_logprob - gsl_sf_log(10),
			data_logprob - gsl_sf_log(30), data_logprob - gsl_sf_log(30),
			data_logprob - gsl_sf_log(100), data_logprob - gsl_sf_log(100)
       );
       printf("\nbe careful.\n");
}

double calc_mcmc_error(const double mean, const char * filename,
		unsigned long batchsize) {
	double v;
	unsigned long n = 0;
	int nbatches = 0;
	FILE * f = openfile(filename);
	double batchsum = 0;
	double batchmean;
	double batcherror;
	double errorsum = 0;

	while (!feof(f)) {
		if (fscanf(f, "%lf", &v) == 1) {
			n++;
			batchsum += v;
			if (n % batchsize == batchsize - 1) {
				batchmean = batchsum / batchsize;
				/*printf("batchmean: %f (batchsize = %lu, batchnr %d)\n",
				 batchmean, batchsize, nbatches);*/
				batcherror = /*batchsize * */pow(batchmean - mean, 2);
				errorsum += batcherror;
				batchsum = 0;
				nbatches++;
			}
		}
	}
	return sqrt(errorsum / nbatches);
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
	double iter;
	double mean;
	double sigma;
	double mcmcerror;
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
	h = create_hist(nbins, gsl_vector_get(min, 0), gsl_vector_get(max, 0));
	debug("filling histogram... ");
	for (i = 0; i < filecount; i++) {
		printf("reading values: chain %3d parameter %s   \r", i, paramname);
		fflush(stdout);
		dump_s("with file", filenames[i]);
		append_to_hists(&h, 1, filenames[i]);
	}
	iter = gsl_histogram_sum(h);
	gsl_histogram_scale(h, (gsl_vector_get(max, 0) - gsl_vector_get(min, 0))
			/ nbins / iter);

	outfile = fopen(outfilename, "w");
	debug("writing histogram... ");
	assert(outfile != NULL);
	gsl_histogram_fprintf(outfile, h, DUMP_FORMAT, DUMP_FORMAT);
	dump_s("histogram file done", outfilename);
	fclose(outfile);

	mean = gsl_histogram_mean(h);
	sigma = gsl_histogram_sigma(h);
	for (i = 0; i < filecount; i++) {
		mcmcerror = calc_mcmc_error(mean, filenames[i], sqrt(iter));
		printf("mcmc error estimate of %s: %f %s\n", paramname, mcmcerror,
				(mcmcerror > sigma * 0.01 ? "** high!" : " (ok)"));
		free(filenames[i]);
	}
	printf("Note: Include a error estimate in your publication!\n");
	free(filenames);

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

