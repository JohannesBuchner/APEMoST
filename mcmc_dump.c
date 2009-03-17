#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "debug.h"

/**
 * if you find the precision too low, increase here.
 * if you find the program too slow, decrease here.
 */
#define DUMP_PRECISION ".3"

void mcmc_dump(const gsl_vector * x_dat, const gsl_vector * y_dat,
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
		fprintf(output, "%" DUMP_PRECISION "e\t%" DUMP_PRECISION "e\n", x, y);
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

void mcmc_dump_probabilities(mcmc * m, unsigned int n_values) {
	(void)m;
	(void)n_values;
	/*get_hist(m, 0, 0);*/
}

