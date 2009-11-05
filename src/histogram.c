#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "gsl_helper.h"
#include "utils.h"
#include "debug.h"
#include "histogram.h"

gsl_histogram * create_hist(int nbins, double min, double max) {
	gsl_histogram * h;

	h = gsl_histogram_alloc(nbins);
	gsl_histogram_set_ranges_uniform(h, min, max);

	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += (max - min) / 10000;
	return h;
}

gsl_histogram * calc_hist(const gsl_vector * v, int nbins) {
	double max;
	double min;
	unsigned int i;
	double binwidth;
	double sum = 0;
	double val;
	gsl_histogram * h;

	gsl_vector_minmax(v, &min, &max);
	binwidth = (max - min) / nbins;
	dump_d("min", min);
	dump_d("max", max);

	debug("allocating the histogram");
	h = gsl_histogram_alloc(v->size);
	debug("setting range");
	require(gsl_histogram_set_ranges_uniform (h, min, max));

	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += 1;

	debug("summing up");
	for (i = 0; i < v->size; i++) {
		val = gsl_vector_get(v, i);
		sum += val;
		require(gsl_histogram_increment (h, val));
	}
	debug("scaling");
	/* double gsl_histogram_sum (const gsl_histogram * h) */
	require(gsl_histogram_scale (h, 1/sum));
	debug("done");
	return h;
}

void append_to_hists(gsl_histogram ** hists, unsigned int n,
		const char * filename) {
	FILE * input;
	unsigned int i;
	int col;
	double v;
	int line = 0;

	input = openfile(filename);

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

void find_min_max(char * filename, gsl_vector * min, gsl_vector * max) {
	int col;
	int i;
	int n = min->size;
	int line = 0;
	double v;
	FILE * input;
	assert(min->size == max->size);
	input = openfile(filename);
	for (i = 0; i < n; i++) {
		col = fscanf(input, "%lf", &v);
		if (col != 1) {
			fprintf(stderr, "field could not be read: %d, line %d in %s\n", i
					+ 1, line + 1, filename);
			exit(1);
		}
		gsl_vector_set(min, i, v);
		gsl_vector_set(max, i, v);
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
			if (gsl_vector_get(min, i) > v)
				gsl_vector_set(min, i, v);
			if (gsl_vector_get(max, i) < v)
				gsl_vector_set(max, i, v);
		}
		line++;
	}
	dump_i("read lines", line);
	fclose(input);
}

void update_min_max(char * filename, gsl_vector * min, gsl_vector * max) {
	gsl_vector * min2 = dup_vector(min);
	gsl_vector * max2 = dup_vector(max);
	find_min_max(filename, min2, max2);
	max_vector(min, min2);
	min_vector(max, max2);
	gsl_vector_free(min2);
	gsl_vector_free(max2);
}
