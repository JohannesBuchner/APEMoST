#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include "debug.h"
#include "utils.h"

#define NCOLUMNS 5
#define MAX_LINE_LENGTH 256

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
	gsl_vector_set(m->params_step, i, step * (max - min));
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
	gsl_vector_set((gsl_vector *) m->x_dat, i, x);
	gsl_vector_set((gsl_vector *) m->y_dat, i, y);
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
			fprintf(stderr, "Line %d of %s is of incorrect format.\n", i + 1,
					filename);
			exit(1);
		}
	}
	dump_i("loaded data points", npoints);

}

char * my_strdup(const char * s) {
	char *buf = calloc(strlen(s) + 1, sizeof(char));
	if (buf != NULL)
		strcpy(buf, s);
	return buf;
}

void mcmc_load_data(mcmc * m, const char * datafilename) {
	dump_s("looking for data in file", datafilename);
	load_data(m, datafilename);
	mcmc_check(m);
}

void mcmc_reuse_data(mcmc * m, const mcmc * m_orig) {
	debug("reusing data from other struct");
	m->x_dat = m_orig->x_dat;
	assert(m->x_dat != NULL);
	m->y_dat = m_orig->y_dat;
	assert(m->y_dat != NULL);
	m->model = gsl_vector_alloc(m->y_dat->size);
	assert(m->model != NULL);
	mcmc_check(m);
}

mcmc * mcmc_load_params(const char * filename) {
	mcmc * m;
	FILE * input;
	int currentline = 0;
	const int pretext = 0;
	int r;

	int nlines;
	nlines = countlines(filename);
	dump_i("number of lines", nlines);

	m = mcmc_init(nlines - pretext);
	input = openfile(filename);
	while (currentline < nlines) {
		if (load_parameter(m, input, currentline - pretext) != 0) {
			fprintf(stderr, "Line %d of %s is of incorrect format.\n",
					currentline + 1, filename);
			exit(1);
		}
		currentline++;
	}
	r = fclose(input);
	assert (r == 0);
	return m;
}
mcmc * mcmc_load(const char * filename, const char * datafilename) {
	mcmc * m = mcmc_load_params(filename);
	mcmc_load_data(m, datafilename);
	return m;
}

