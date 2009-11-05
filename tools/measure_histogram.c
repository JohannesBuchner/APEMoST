#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "gsl_helper.h"
#include "debug.h"
#include "utils.h"
#include "histogram.h"

void usage(char * progname) {
	fprintf(stderr, "%s: SYNOPSIS: min max file\n"
		"\n"
		"\tfile\tfile containing one row of data\n"
		"\tmin\tminimal value\n"
		"\tmax\tmaximum value\n"
		"\n"
		"This program interpretes a histogram from datafiles.\n"
		"It prints out the median, deviation and probability of "
		"the features.\n", progname);
}

void output_hist(gsl_histogram * h) {
	double lower;
	double upper;
	unsigned int count;
	unsigned int i;

	for (i = 0; i < gsl_histogram_bins(h); i++) {
		gsl_histogram_get_range(h, i, &lower, &upper);
		count = gsl_histogram_get(h, i);

		printf("%f..%f\t:\t%10d\n", lower, upper, count);
	}
}

void analyse_part(gsl_histogram * h, unsigned int left, unsigned int right,
		unsigned long nvalues, double * median, double * leftquartile,
		double * rightquartile, double * percentage) {
	unsigned long n = 0;
	unsigned long c = 0;
	unsigned int i;
	double min;
	double max;
	for (i = left; i <= right; i++) {
		n += gsl_histogram_get(h, i);
	}
	dump_ul("number of values in this part", n);
	*percentage = 1.0 * n / nvalues;
	dump_d("equal to percentage of", *percentage);
	for (i = left; i <= right; i++) {
		c += gsl_histogram_get(h, i);
		gsl_histogram_get_range(h, i, &min, &max);
		if (c < (right - left) * 1 / 4) {
			*leftquartile = min;
		}
		if (c < (right - left) * 2 / 4) {
			*median = min;
		}
		if (c < (right - left) * 3 / 4) {
			*rightquartile = min;
		}
	}
	dump_d("leftquartile", *leftquartile);
	dump_d("median", *median);
	dump_d("rightquartile", *rightquartile);
}


void run(char * filename, unsigned long nmaxvalues, double min, double max) {
	unsigned int i;
	gsl_histogram * h;
	unsigned int nbins;
	unsigned long nvalues;
	unsigned int empty_in_a_row = 0;
	unsigned int empty_in_a_row_needed;
	unsigned int left;
	unsigned int right;
	double median, leftquartile, rightquartile, percentage;
	gsl_vector * medians = gsl_vector_alloc(100);
	gsl_vector * leftquartiles = gsl_vector_alloc(100);
	gsl_vector * rightquartiles = gsl_vector_alloc(100);
	gsl_vector * percentages = gsl_vector_alloc(100);
	gsl_vector * vectors[4];
	unsigned int npeaks = 0;
	vectors[0] = percentages;
	vectors[1] = medians;
	vectors[2] = leftquartiles;
	vectors[3] = rightquartiles;

	dump_ul("nmaxvalues", nmaxvalues);
	h = create_hist(nmaxvalues, min, max);
	append_to_hists(&h, 1, filename);
	nvalues = gsl_histogram_sum(h);
	debug("histogram ready!");

	nbins = gsl_histogram_bins(h);
	dump_ui("nbins", nbins);

	/* 1% of parameter space */
	empty_in_a_row_needed = nbins / 100;

	dump_ui("nbins", nbins);
	right = 0;
	while (1) {
		left = right;
		while (left < nbins && (gsl_histogram_get(h, left) == 0))
			left++;
		left--;
		if (left == nbins - 1)
			break;
		dump_ui("found left side", left);
		right = left + 1;
		for (; right < nbins; right++) {
			if (gsl_histogram_get(h, right) == 0) {
				empty_in_a_row++;
				if (empty_in_a_row > empty_in_a_row_needed)
					break;
			} else {
				empty_in_a_row = 0;
			}
		}
		right--;
		dump_ui("found right side", right);
		analyse_part(h, left, right, nvalues, &median, &leftquartile,
				&rightquartile, &percentage);
		gsl_vector_set(medians, npeaks, median);
		gsl_vector_set(leftquartiles, npeaks, leftquartile);
		gsl_vector_set(rightquartiles, npeaks, rightquartile);
		gsl_vector_set(percentages, npeaks, percentage);
		npeaks++;

		if (right == nbins - 1)
			break;
	}

#ifndef NOSORT
	sort(vectors, 4, npeaks);
#endif
	printf("median\t-\t+\tpercent\n");
	fflush(NULL);
	for (i = 0; i < npeaks; i++) {
		printf("%f\t%f\t%f\t%f\n", gsl_vector_get(medians, i), gsl_vector_get(
				medians, i) - gsl_vector_get(leftquartiles, i), gsl_vector_get(
				rightquartiles, i) - gsl_vector_get(medians, i),
				gsl_vector_get(percentages, i));
	}

	free(h);
}

/* see usage. nvalues min max file */
int main(int argc, char ** argv) {
	if (argc <= 3) {
		usage(argv[0]);
	} else {
		run(argv[3], countlines(argv[3]), atof(argv[1]), atof(argv[2]));
	}
	return 0;
}
