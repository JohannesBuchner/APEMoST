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

#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include "debug.h"
#include "utils.h"

#define MAX_LINE_LENGTH 256
#ifdef DEBUGPARSER
#define IFDEBUGPARSER if(1)
#else
#define IFDEBUGPARSER if(0)
#endif

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
	char * descr = (char*) mem_calloc(MAX_LINE_LENGTH, sizeof(char));
	IFDEBUGPARSER
	dump_i("parsing line", i);

	col = fscanf(input, "%lf\t%lf\t%lf\t%s\t%lf\n", &start, &min, &max, descr,
			&step);
	if (col != 5) {
		fprintf(stderr, "only %d fields matched.\n", col);
		return 1;
	}
	if (!(descr != NULL && strnlen(descr, MAX_LINE_LENGTH) > 0 && strnlen(
			descr, MAX_LINE_LENGTH) < MAX_LINE_LENGTH)) {
		fprintf(stderr, "description invalid: %s\n", descr);
		return 1;
	}
	IFDEBUGPARSER
	debug("setting values");
	gsl_vector_set(m->params, i, start);
	gsl_vector_set(m->params_best, i, start);
	if (min > max) {
		fprintf(stderr, "min(%f) < max(%f)\n", min, max);
		return 1;
	}
	if (start > max) {
		fprintf(stderr, "start(%f) > max(%f)\n", start, max);
		return 1;
	}
	if (start < min) {
		fprintf(stderr, "start(%f) < min(%f)\n", start, min);
		return 1;
	}
	if (step < 0) {
		step = (max - min) * 0.1;
		dump_d("using auto step size, 10% of parameter space", step);
	}
	gsl_vector_set(m->params_min, i, min);
	gsl_vector_set(m->params_max, i, max);
	m->params_descr[i] = descr;
	gsl_vector_set(m->params_step, i, step /* * (max - min) */);
	IFDEBUGPARSER
	debug("setting values done.");
	return 0;
}

static void load_data(mcmc * m, const char * filename) {
	FILE * input;
	int npoints = countlines(filename);
	int ndim = get_column_count(filename);
	gsl_matrix * data;
	IFDEBUGPARSER
	dump_i("lines", npoints);
	IFDEBUGPARSER
	dump_i("dimensions", ndim);

	data = gsl_matrix_alloc(npoints, ndim);

	input = openfile(filename);
	if (gsl_matrix_fscanf(input, data) != 0) {
		fprintf(stderr,
				"error reading input data. Perhaps inconsistent format?\n");
		fprintf(stderr, "tried to read %d x %d.\n", ndim, npoints);
		exit(3);
	}
	assert(fclose(input) == 0);

	m->data = data;

	IFDEBUGPARSER
	dump_i("loaded data points", npoints);
}

char * my_strdup(const char * s) {
	char *buf = mem_calloc(strlen(s) + 1, sizeof(char));
	if (buf != NULL)
		strcpy(buf, s);
	return buf;
}

void mcmc_load_data(mcmc * m, const char * datafilename) {
	IFDEBUGPARSER
	dump_s("looking for data in file", datafilename);
	load_data(m, datafilename);
	mcmc_check(m);
	IFDEBUGPARSER
	debug("loading data successful")
}

void mcmc_reuse_data(mcmc * m, const mcmc * m_orig) {
	IFDEBUGPARSER
	debug("reusing data from other struct");
	assert(m_orig->data != NULL);
	m->data = m_orig->data;
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
	IFDEBUGPARSER
	dump_i("number of lines", nlines);

	m = mcmc_init(nlines - pretext);
	IFDEBUGPARSER
	debug("opening params file");
	input = openfile(filename);
	while (currentline < nlines) {
		IFDEBUGPARSER
		dump_i("reading parameter", currentline);
		if (load_parameter(m, input, currentline - pretext) != 0) {
			fprintf(stderr, "Line %d of %s is of incorrect format.\n",
					currentline + 1, filename);
			exit(1);
		}
		currentline++;
	}
	IFDEBUGPARSER
	debug("closing params file");
	r = fclose(input);
	assert (r == 0);
	IFDEBUGPARSER
	debug("finished reading params file");
	return m;
}
mcmc * mcmc_load(const char * filename, const char * datafilename) {
	mcmc * m = mcmc_load_params(filename);
	mcmc_load_data(m, datafilename);
	return m;
}

