#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "debug.h"

/**
 * if you find the precision too low, increase here.
 * if you find the program too slow, decrease here.
 */
#define DUMP_FORMAT "%.3e"

static void mcmc_dump(const gsl_vector * x_dat, const gsl_vector * y_dat,
		const char * filename) {
	unsigned long i;
	double x, y;
	FILE * output;
	assert(x_dat->size == y_dat->size);
	output = fopen(filename, "w");
	assert(output != NULL);
	dump_s("writing dump to file", filename);
	for (i = 0; i < x_dat->size; i++) {
		x = gsl_vector_get(x_dat, i);
		y = gsl_vector_get(y_dat, i);
		fprintf(output, DUMP_FORMAT "\t" DUMP_FORMAT "\n", x, y);
	}
	assert (fclose(output) == 0);
	debug("writing dump to file done");
}

void mcmc_dump_model(mcmc * m) {
	mcmc_dump(m->x_dat, m->model, "model.dat.dump");
}

void mcmc_dump_y_dat(mcmc * m) {
	mcmc_dump(m->x_dat, m->y_dat, "y_dat.dat.dump");
}

void mcmc_dump_probabilities(mcmc * m, int n_values) {
	unsigned int i;
	unsigned int j;
	FILE ** files     = (FILE**) calloc(m->n_par, sizeof(FILE*));
	char ** filenames = (char**) calloc(m->n_par, sizeof(char*));
	for (i = 0; i < m->n_par; i++) {
		filenames[i] = (char*)calloc(strlen(m->params_descr[i]) + strlen(".prob.dump") + 1, sizeof(char));
		sprintf(filenames[i], "%s.prob.dump", m->params_descr[i]);
		dump_s("writing probability/distribution to file", filenames[i]);
		files[i] = fopen(filenames[i], "w");
		assert(files[i] != NULL);
	}

	j = 0;
	if (n_values > 0 && (unsigned int)n_values < m->n_iter)
		j = m->n_iter - n_values;
	dump_i("starting at iteration", j);
	dump_ul("until iteration", m->n_iter);
	for (; j < m->n_iter; j++) {
		dump_i("iteration", j)
		for (i = 0; i < m->n_par; i++) {
			dump_s("variable", m->params_descr[i]);
			fprintf(files[i], DUMP_FORMAT "\n", gsl_vector_get(
					m->params_distr[j], i));
		}
	}
	debug("done writing. closing ...");
	for (i = 0; i < m->n_par; i++) {
		free(filenames[i]);
		assert(fclose(files[i]) == 0);
	}
	free(filenames);
	free(files);
	debug("done.");
}