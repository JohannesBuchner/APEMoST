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

void usage(char * progname) {
	fprintf(stderr, "%s: SYNOPSIS: [-s] [-1] file1 file2 ...\n"
		"\n"
		"\ts\tuse a square sum\n"
		"\t1\tat the end, sum all columns up to one value\n"
		"\n"
		"This program sums up the values of the given files. They should contain \n"
		"float values in one or more columns\n"
		"\n", progname);
}

void sum_up(gsl_vector * sum, unsigned int n, const char * filename, int square) {
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
			if (square)
				gsl_vector_add_constant(sum, v * v);
			else
				gsl_vector_add_constant(sum, v);
		}
		line++;
	}
	dump_i("read lines", line);
	fclose(input);
}

void run(char ** filenames, unsigned int filecount, int square, int one_value) {
	unsigned int i;
	gsl_vector * sum;
	unsigned int ncolumns;

	debug("looking for number of columns ...");
	ncolumns = get_column_count(filenames[0]);
	dump_i("number of columns", ncolumns);
	sum = gsl_vector_alloc(ncolumns);
	gsl_vector_set_zero(sum);

	debug("filling histograms ... ");
	for (i = 0; i < filecount; i++) {
		dump_s("with file", filenames[i]);
		sum_up(sum, ncolumns, filenames[i], square);
	}

	if (one_value) {
		if (square)
			printf("%e", calc_vector_squaresum(sum));
		else
			printf("%e", calc_vector_sum(sum));
	} else {
		for (i = 0; i < ncolumns - 1; i++) {
			printf("%e\t", gsl_vector_get(sum, i));
		}
		printf("%e", gsl_vector_get(sum, i));
	}
	printf("\n");
	gsl_vector_free(sum);
}

int main(int argc, char ** argv) {
	int square = 0;
	int onevalue = 0;
	if (argc <= 1) {
		usage(argv[0]);
	} else {
		argv++;
		argc--;

		if (argc > 1 && strcmp(argv[0], "-s") == 0) {
			square = 1;
			argv++;
			argc--;
		}
		if (argc > 1 && strcmp(argv[0], "-1") == 0) {
			onevalue = 1;
			argv++;
			argc--;
		}
		run(argv, argc, square, onevalue);
	}
	return 0;
}
