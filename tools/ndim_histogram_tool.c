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

typedef struct {
	long * values;
	/* minima */
	gsl_vector * min;
	/* maxima */
	gsl_vector * max;
	unsigned int nbins;
	unsigned int size;
	unsigned int dimensions;
} ndim_histogram;

int get_index_for(ndim_histogram * h, gsl_vector * value) {
	unsigned int i;
	unsigned int index = 0;
	unsigned int subindex;

	assert(value->size == h->dimensions);
	for (i = 0; i < h->dimensions; i++) {
		subindex = (h->nbins - 1) * ((gsl_vector_get(value, i)
				- gsl_vector_get(h->min, i)) / (gsl_vector_get(h->max, i)
				- gsl_vector_get(h->min, i)));
		index = index * h->nbins + subindex;
	}
	return index;
}

gsl_vector * get_lower_corner_for_index(ndim_histogram * h, unsigned int index) {
	int subindex;
	int i;
	gsl_vector * lower = gsl_vector_alloc(h->dimensions);

	for (i = h->dimensions - 1; i >= 0; i--) {
		subindex = index % h->nbins;
		gsl_vector_set(lower, i, subindex * (gsl_vector_get(h->max, i)
				- gsl_vector_get(h->min, i)) / (h->nbins - 1) + gsl_vector_get(
				h->min, i));
		index /= h->nbins;
	}
	return lower;
}

gsl_vector * get_upper_corner_for_index(ndim_histogram * h, unsigned int index) {
	int subindex;
	int i;
	gsl_vector * upper = gsl_vector_alloc(h->dimensions);

	for (i = h->dimensions - 1; i >= 0; i--) {
		subindex = index % h->nbins;
		gsl_vector_set(upper, i, (subindex + 1) * (gsl_vector_get(h->max, i)
				- gsl_vector_get(h->min, i)) / (h->nbins - 1) + gsl_vector_get(
				h->min, i));
		index /= h->nbins;
	}
	return upper;
}

ndim_histogram * ndim_histogram_alloc(unsigned int nbins, gsl_vector * min,
		gsl_vector * max) {
	ndim_histogram * h;

	h = (ndim_histogram*) malloc(sizeof(ndim_histogram));
	assert(h!=NULL);
	h->min = min;
	h->max = max;
	h->nbins = nbins;
	assert(min->size == max->size);
	h->dimensions = min->size;
	if (pow(h->nbins, h->dimensions) > INT_MAX) {
		printf("I can not make %e cubes. choose less bins.\n", pow(h->nbins,
				h->dimensions));
		printf("Maximum: %d^%e, %d\n", INT_MAX, (1 / (double) h->dimensions),
				(int) pow((double) INT_MAX, 1 / (double) h->dimensions));
		exit(1);
	}
	h->size = pow(h->nbins, h->dimensions);

	dump_i("used size", h->size);

	h->values = (long*) calloc(h->size, sizeof(long));
	if (h->values == NULL) {
		perror("could not allocate that much");
		exit(1);
	} else {
		debug("allocating succeeded.");
	}

	return h;
}
void ndim_histogram_free(ndim_histogram * h) {
	free(h->values);
	free(h);
}

void ndim_histogram_increment(ndim_histogram * h, gsl_vector * current) {
	h->values[get_index_for(h, current)]++;
}

long ndim_histogram_get(ndim_histogram * h, int index) {
	return h->values[index];
}

void append_to_hist(ndim_histogram * h, const char * filename) {
	FILE * input;
	unsigned int i;
	int col;
	double v;
	int line = 0;
	gsl_vector * current = gsl_vector_alloc(h->dimensions);

	input = fopen(filename, "r");
	if (input == NULL) {
		fprintf(stderr, "error opening file %s\n", filename);
		perror("file could not be opened");
		exit(1);
	}
	while (!feof(input)) {
		for (i = 0; i < h->dimensions; i++) {
			col = fscanf(input, "%lf", &v);
			if (col != 1) {
				if (feof(input))
					break;
				fprintf(stderr, "field could not be read: %d, line %d in %s\n",
						i + 1, line + 1, filename);
				exit(1);
			}
			gsl_vector_set(current, i, v);
		}
		ndim_histogram_increment(h, current);
		line++;
	}
	dump_i("read lines", line);
	fclose(input);
}

void output_hist(ndim_histogram * h) {
	unsigned int i;
	long value;

	for (i = 0; i < h->size; i++) {
		value = ndim_histogram_get(h, i);
		if (value == 0)
			continue;
		dump_vector(get_lower_corner_for_index(h, i));
		printf("..");
		dump_vector(get_upper_corner_for_index(h, i));
		printf("\t%lu\n", value);
	}
}

void run(unsigned int nbins, char * filename) {
	gsl_vector * min;
	gsl_vector * max;
	ndim_histogram * h;
	unsigned int ncolumns;

	debug("looking for number of columns ...");
	ncolumns = get_column_count(filename);
	dump_i("number of columns", ncolumns);
	min = gsl_vector_alloc(ncolumns);
	max = gsl_vector_alloc(ncolumns);

	debug("finding minima/maxima ...");
	dump_s("in file", filename);
	find_min_max(filename, min, max);
	dump_v("minima", min);
	dump_v("maxima", max);
	debug("creating histogram cube ...");
	h = ndim_histogram_alloc(nbins, min, max);

	debug("filling histogram ... ");
	append_to_hist(h, filename);

	output_hist(h);
	ndim_histogram_free(h);
}

void usage(char * progname) {
	fprintf(
			stderr,
			"%s: SYNOPSIS: nbins file...\n"
				"\n"
				"\tnbins\tNumber of bins to use for each dimension\n"
				"\tfile\tfile to include. these should contain float values in one or more rows\n"
				"\n"
				"This program calculates a n-dimensional histogram cube from datafiles.\n",
			progname);
}

int main(int argc, char ** argv) {
	if (argc != 3) {
		usage(argv[0]);
	} else {
		run(atoi(argv[1]), argv[2]);
	}
	return 0;
}
