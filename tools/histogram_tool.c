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
#include "histogram.c"

gsl_histogram * create_hist(int nbins, double min, double max) {
	gsl_histogram * h;

	h = gsl_histogram_alloc(nbins);
	gsl_histogram_set_ranges_uniform(h, min, max);

	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += 1;
	return h;
}

void usage(char * progname) {
	fprintf(
			stderr,
			"%s: SYNOPSIS: nbins file1 file2 ...\n"
				"\n"
				"\tnbins\tNumber of bins to use\n"
				"\tfiles\tfiles to include. these should contain float values in one or more rows\n"
				"\n"
				"This program calculates a histogram from datafiles.\n",
			progname);
}

/**
 * s = rep + sep + rep + sep + ... + rep + end
 *                (n times)
 */
char * str_rep_join(const char * rep, const char * sep, const char * end, int n) {
	int i;
	char * s = (char*) calloc(strlen(sep) * (n - 1) + strlen(rep) * n + strlen(
			end) + 1, sizeof(char));
	sprintf(s, "%s", rep);
	for (i = 1; i < n; i++) {
		sprintf(s + strlen(s), "%s%s", sep, rep);
	}
	sprintf(s + strlen(s), "%s", end);
	return s;
}

void append_to_hists(gsl_histogram ** hists, unsigned int n,
		const char * filename) {
	FILE * input;
	unsigned int i;
	int col;
	double v;
	int line = 0;

	input = fopen(filename, "r");
	if (input == NULL) {
		fprintf(stderr, "error opening file %s\n", filename);
		perror("file could not be opened");
		exit(1);
	}
	while (!feof(input)) {
		for (i = 0; i < n; i++) {
			col = fscanf(input, "%lf", &v);
			if (col != 1) {
				if (feof(input))
					break;
				fprintf(stderr, "field could not be read: %d, line %d in %s\n",
						i + 1, line + 1, filename);
				exit(1);
			}
			gsl_histogram_increment(hists[i], v);
		}
		line++;
	}
	dump_i("read lines", line);
	fclose(input);
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

void run(unsigned int nbins, char ** filenames, unsigned int filecount) {
	unsigned int i;
	gsl_vector * min;
	gsl_vector * max;
	gsl_histogram ** hists;
	unsigned int ncolumns;

	debug("looking for number of columns ...");
	ncolumns = get_column_count(filenames[0]);
	dump_i("number of columns", ncolumns);
	min = gsl_vector_alloc(ncolumns);
	max = gsl_vector_alloc(ncolumns);
	hists = (gsl_histogram **) calloc(ncolumns, sizeof(gsl_histogram *));

	debug("finding minima/maxima ...");
	dump_s("in file", filenames[0]);
	find_min_max(filenames[0], min, max);
	for (i = 1; i < filecount; i++) {
		if (ncolumns != get_column_count(filenames[i])) {
			fprintf(stderr, "number of columns different in file %s: %i vs %i in %s\n",
					filenames[i], ncolumns, get_column_count(filenames[i]),
					filenames[0]);
			exit(1);
		}
		dump_s("in file", filenames[i]);
		update_min_max(filenames[i], min, max);
	}
	dump_v("minima", min);
	dump_v("maxima", max);
	debug("creating histograms ...");
	for (i = 0; i < ncolumns; i++) {
		hists[i] = create_hist(nbins, gsl_vector_get(min, i), gsl_vector_get(
				max, i));
	}
	debug("filling histograms ... ");
	for (i = 0; i < filecount; i++) {
		dump_s("with file", filenames[i]);
		append_to_hists(hists, ncolumns, filenames[i]);
	}

	for (i = 0; i < ncolumns; i++) {
		printf("column %d:\n", i + 1);
		output_hist(hists[i]);
		gsl_histogram_free(hists[i]);
		printf("\n");
	}
	free(hists);
}

int main(int argc, char ** argv) {
	if (argc <= 2) {
		usage(argv[0]);
	} else {
		run(atoi(argv[1]), &argv[2], argc - 2);
	}
	return 0;
}
