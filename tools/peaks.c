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
 
/**
 * peaks (this program) measures the median and quartiles of distinctive 
 * features of a "histogram".
 * In reality, no histogram is created from the input file, but rather the
 * file is sorted and the quartiles directly read out, which makes the 
 * calculations exact.
 * "Distinctive features" are peaks that have 1% of parameter space in between.
 * The output also includes the share each peak accounts for.
 */

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

void analyse_part(gsl_vector * v, unsigned int left, unsigned int right,
		unsigned long nvalues, double * median, double * leftquartile,
		double * rightquartile, double * percentage) {
	unsigned long n = right - left + 1;
	unsigned long c = 0;
	unsigned int i;
	dump_ul("number of values in this part", n);
	*percentage = 1.0 * n / nvalues;
	dump_d("equal to percentage of", *percentage);
	for (i = left; i <= right; i++) {
		c++;
		if (c <= n * 1 / 4) {
			*leftquartile = gsl_vector_get(v, i);
		}
		if (c <= n * 2 / 4) {
			*median = gsl_vector_get(v, i);
		}
		if (c <= n * 3 / 4) {
			*rightquartile = gsl_vector_get(v, i);
		}
	}
	dump_d("leftquartile", *leftquartile);
	dump_d("median", *median);
	dump_d("rightquartile", *rightquartile);
}

unsigned int fill_vector(gsl_vector * values, char * filename, double min,
		double max) {
	FILE * input;
	int col;
	double v;
	int line = 0;
	unsigned int i = 0;

	input = openfile(filename);

	while (!feof(input)) {
		col = fscanf(input, "%lf", &v);
		if (col != 1) {
			if (feof(input))
				break;
			fprintf(stderr, "field could not be read: %d, line %d in %s\n", i
					+ 1, line + 1, filename);
			exit(1);
		}
		if (v >= min && v <= max) {
			gsl_vector_set(values, i, v);
			i++;
		}
		line++;
	}
	dump_i("read lines", line);
	dump_i("accepted lines", i);
	fclose(input);
	return i;
}

int double_comp(const void * av, const void * bv) {
	double a = *(double*) av;
	double b = *(double*) bv;
	if (a < b)
		return -1;
	else if (a == b)
		return 0;
	else
		return 1;
}

void run(char * filename, unsigned long nmaxvalues, double min, double max) {
	unsigned int i;
	unsigned int nvalues;
	double empty_space_needed;
	unsigned int left;
	unsigned int right;
	double median, leftquartile, rightquartile, percentage;
	gsl_vector * values = gsl_vector_alloc(nmaxvalues);
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
	nvalues = fill_vector(values, filename, min, max);
	dump_ul("nvalues", nmaxvalues);
	debug("data ready!");
	qsort(values->data, nvalues, sizeof(double), double_comp);
	debug("data sorted!");

	/* 0.1% of parameter space */
	empty_space_needed = (max - min) / 100;
	dump_d("empty_space_needed", empty_space_needed);
	dump_d("min", min);
	dump_d("max", max);

	assert(gsl_vector_get(values, 0) >= min);
	assert(gsl_vector_get(values, 0) <= max);
	assert(gsl_vector_get(values, nvalues - 1) >= min);
	assert(gsl_vector_get(values, nvalues - 1) <= max);

	left = 0;
	while (1) {
		dump_ui("found left side", left);
		dump_d("left side", gsl_vector_get(values, left));
		right = left + 1;
		for (; right < nvalues; right++) {
			if (gsl_vector_get(values, right) - gsl_vector_get(values, right
					- 1) > empty_space_needed) {
				dump_d("gap till next", gsl_vector_get(values, right) - gsl_vector_get(values, right
								- 1));
				break;
			}
		}
		right--;
		dump_ui("found right side", right);
		dump_d("right side", gsl_vector_get(values, right));
		analyse_part(values, left, right, nvalues, &median, &leftquartile,
				&rightquartile, &percentage);
		gsl_vector_set(medians, npeaks, median);
		gsl_vector_set(leftquartiles, npeaks, leftquartile);
		gsl_vector_set(rightquartiles, npeaks, rightquartile);
		gsl_vector_set(percentages, npeaks, percentage);
		npeaks++;
		assert(npeaks < 100);

		if (right == nvalues - 1)
			break;
		left = right + 1;
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

	gsl_vector_free(values);
	gsl_vector_free(medians);
	gsl_vector_free(leftquartiles);
	gsl_vector_free(rightquartiles);
	gsl_vector_free(percentages);
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

