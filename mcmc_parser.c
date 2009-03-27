#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include "debug.h"

#define NCOLUMNS 5
#define MAX_LINE_LENGTH 256

static FILE * openfile(const char * filename) {
	FILE * input = fopen(filename, "r");
	if (input == NULL) {
		fprintf(stderr, "error opening file %s\n", filename);
		perror("file could not be opened");
		exit(1);
	}
	return input;
}

unsigned int countlines(const char * filename) {
	int nlines = 0;
	int c;
	FILE * input = openfile(filename);
	while (1) {
		c = fgetc(input);
		if (c == '\n') {
			nlines++;
		}
		if (c == EOF)
			break;
	}
	assert (fclose(input) == 0);
	return nlines;
}

static int strnlen(char * s, int maxlen) {
	int i;
	for (i = 0; i < maxlen && s[i] != 0; i++)
		;
	return i;
}

/**
 * returns 0 on success.
 */
static int load_parameter(mcmc * m, FILE * input, int i) {
	int col = 0;
	double start;
	double min;
	double max;
	double step;
	char * descr = (char*) calloc(MAX_LINE_LENGTH, sizeof(char));
	dump_i("parsing line", i);

	col = fscanf(input, "%lf\t%lf\t%lf\t%s\t%lf\n", &start, &min, &max, descr,
			&step);
	if (col != 5) {
		fprintf(stderr, "only %d fields matched.\n", col);
		return 1;
	}
	if (step < 0 || step > 1) {
		fprintf(stderr,
				"start (column 1) must be between 0 and 1. currently: %f.\n",
				start);
		return 1;
	}
	if (!(descr != NULL && strnlen(descr, MAX_LINE_LENGTH) > 0 && strnlen(
			descr, MAX_LINE_LENGTH) < MAX_LINE_LENGTH)) {
		fprintf(stderr, "description invalid: %s\n", descr);
		return 1;
	}
	debug("setting values");
	gsl_vector_set(m->params, i, start);
	gsl_vector_set(m->params_best, i, start);
	gsl_vector_set(m->params_min, i, min);
	gsl_vector_set(m->params_max, i, max);
	m->params_descr[i] = descr;
	gsl_vector_set(m->params_step, i, min + start * (max - min));
	debug("setting values done.");
	return 0;
}

static int load_datapoint(mcmc * m, FILE * input, int i) {
	double x;
	double y;
	int col;

	col = fscanf(input, "%lf%lf", &x, &y);
	if (col != 2) {
		fprintf(stderr, "only %d fields matched.\n", col);
		return 1;
	}
	gsl_vector_set(m->x_dat, i, x);
	gsl_vector_set(m->y_dat, i, y);
	return 0;
}

static void load_data(mcmc * m, const char * filename) {
	FILE * input;
	int i = 0;
	int npoints = countlines(filename);
	dump_i("lines", npoints);

	m->x_dat = gsl_vector_alloc(npoints);
	assert(m->x_dat != NULL);
	m->y_dat = gsl_vector_alloc(npoints);
	assert(m->y_dat != NULL);
	m->model = gsl_vector_alloc(npoints);
	assert(m->model != NULL);

	input = openfile(filename);
	for (i = 0; i < npoints; i++) {
		if (load_datapoint(m, input, i) != 0) {
			fprintf(stderr, "Line %d of %s is of incorrect format.", i + 1,
					filename);
			exit(1);
		}
	}
	dump_i("loaded data points", npoints);

}

char * strdup(const char * s) {
	char *buf = calloc(strlen(s) + 1, sizeof(char));
	if (buf != NULL)
		strcpy(buf, s);
	return buf;
}

mcmc * mcmc_load(const char * filename) {
	mcmc * m;
	FILE * input;
	char datafilename[MAX_LINE_LENGTH];
	char datafilepath[MAX_LINE_LENGTH];
	char * datadir;
	int currentline = 0;
	const int pretext = 1;

	int nlines;
	nlines = countlines(filename);
	dump_i("number of lines", nlines);

	m = mcmc_init(nlines - pretext);
	input = openfile(filename);
	fscanf(input, "%s\n", datafilename);
	currentline++;

	while (currentline < nlines) {
		if (load_parameter(m, input, currentline - pretext) != 0) {
			fprintf(stderr, "Line %d of %s is of incorrect format.",
					currentline + 1, filename);
			exit(1);
		}
		currentline++;
	}
	assert (fclose(input) == 0);

	datadir = dirname(strdup(filename));
	sprintf(datafilepath, "%s/%s", datadir, datafilename);
	dump_s("looking for data in file", datafilepath);

	load_data(m, datafilepath);
	return m;
}

