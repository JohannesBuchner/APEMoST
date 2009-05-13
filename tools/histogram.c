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

unsigned int get_column_count(char * filename) {
	FILE * input;
	char buf[10000];
	unsigned int count = 0;
	int i = 0;
	int at_whitespace = 1;

	input = openfile(filename);
	fgets(buf, 10000, input);
	while (buf[i] != 0) {
		if (isspace(buf[i])) {
			if (at_whitespace == 0)
				count++;
			at_whitespace = 1;
		} else {
			at_whitespace = 0;
		}
		i++;
	}
	if (at_whitespace == 0)
		count++;
	return count;
}
