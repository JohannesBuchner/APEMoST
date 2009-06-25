#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "debug.h"

/**
 * if you find the precision too low, increase here.
 * if you find the program too slow, decrease here.
 */
#define DUMP_FORMAT "%.15e"
#ifdef NODUMP
#define ASSURE_DUMP_ENABLED return;
#else
#define ASSURE_DUMP_ENABLED
#endif

static void mcmc_dump(const gsl_vector * x_dat, const gsl_vector * y_dat,
		const char * filename) {
	unsigned long i;
	int r;
	double x, y;
	FILE * output;
	ASSURE_DUMP_ENABLED;
	assert(x_dat->size == y_dat->size);
	output = fopen(filename, "w");
	assert(output != NULL);
	IFVERBOSE
		dump_s("writing dump to file", filename);
	for (i = 0; i < x_dat->size; i++) {
		x = gsl_vector_get(x_dat, i);
		y = gsl_vector_get(y_dat, i);
		fprintf(output, DUMP_FORMAT "\t" DUMP_FORMAT "\n", x, y);
	}
	r = fclose(output);
	assert (r == 0);
	IFVERBOSE
		debug("writing dump to file done");
}

void mcmc_dump_model(const mcmc * m) {
	mcmc_dump(m->x_dat, m->model, "model.dat.dump");
}

void mcmc_dump_y_dat(const mcmc * m) {
	mcmc_dump(m->x_dat, m->y_dat, "y_dat.dat.dump");
}

void mcmc_open_dump_files(mcmc * m, const char * suffix) {
	unsigned int i;
	char ** filenames = (char**) mem_calloc(m->n_par, sizeof(char*));

	m->files = (FILE**) mem_calloc(m->n_par, sizeof(FILE*));
	ASSURE_DUMP_ENABLED;
	for (i = 0; i < get_n_par(m); i++) {
		filenames[i] = (char*) mem_calloc(strlen(m->params_descr[i]) + strlen(
						suffix) + strlen(".prob.dump") + 1, sizeof(char));
		sprintf(filenames[i], "%s%s.prob.dump", m->params_descr[i], suffix);
		IFVERBOSE
			dump_s("writing probability/distribution to file", filenames[i]);
		m->files[i] = fopen(filenames[i], "w");
		assert(m->files[i] != NULL);
		mem_free(filenames[i]);
	}
	mem_free(filenames);
}

void mcmc_dump_current(const mcmc * m) {
	unsigned int i;
	if (m->files == NULL)
		return;
	for (i = 0; i < get_n_par(m); i++) {
		if (m->files[i] == NULL)
			continue;
		fprintf(m->files[i], DUMP_FORMAT "\n", gsl_vector_get(m->params, i));
	}
}

void mcmc_dump_close(mcmc * m) {
	unsigned int i;
	int r;
	if (m->files == NULL)
		return;
	for (i = 0; i < get_n_par(m); i++) {
		if (m->files[i] == NULL)
			continue;
		r = fclose(m->files[i]);
		assert(r == 0);
	}
	mem_free(m->files);
	m->files = NULL;
}
void mcmc_dump_flush(const mcmc * m) {
	unsigned int i;
	if (m->files == NULL)
		return;
	for (i = 0; i < get_n_par(m); i++) {
		if (m->files[i] == NULL)
			continue;
		fflush(m->files[i]);
	}
}
