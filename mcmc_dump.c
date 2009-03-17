#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "debug.h"


void mcmc_dump(const gsl_vector * x_dat, const gsl_vector * y_dat,
		const char * filename) {
	unsigned long i;
	double x, y;
	FILE * output;
	assert(x_dat->size == y_dat->size);
	output = fopen(filename, "w");
	dump_s("writing dump to file", filename);
	for (i = 0; i < x_dat->size; i++) {
		x = gsl_vector_get(x_dat, i);
		y = gsl_vector_get(y_dat, i);
		fprintf(output, "%f\t%f\n", x, y);
	}
	if (fclose(output) != 0) {
		assert(0);
	}
	debug("writing dump to file done");
}


void mcmc_dump_model(mcmc * m) {
	mcmc_dump(m->x_dat, m->model, "model.dat.dump");
}

void mcmc_dump_y_dat(mcmc * m) {
	mcmc_dump(m->x_dat, m->y_dat, "y_dat.dat.dump");
}
