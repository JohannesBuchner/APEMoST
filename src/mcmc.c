#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include "debug.h"

void init_seed(mcmc * m) {
	gsl_rng_env_setup();
	m->random = gsl_rng_alloc(gsl_rng_default);
}

mcmc * mcmc_init(const unsigned int n_pars) {
	mcmc * m;
	IFSEGV
		debug("allocating mcmc struct");
	m = (mcmc*) mem_malloc(sizeof(mcmc));
	assert(m != NULL);
	m->n_iter = 0;
	m->n_par = n_pars;
	m->accept = 0;
	m->reject = 0;
	m->prob = -1e+10;
	m->prob_best = -1e+10;
	m->files = NULL;

	init_seed(m);

	m->params = gsl_vector_alloc(m->n_par);
	assert(m->params != NULL);
	m->params_best = gsl_vector_alloc(m->n_par);
	assert(m->params_best != NULL);

	m->params_accepts
			= (unsigned long*) mem_calloc(m->n_par, sizeof(unsigned long));
	assert(m->params_accepts != NULL);
	m->params_rejects
			= (unsigned long*) mem_calloc(m->n_par, sizeof(unsigned long));
	assert(m->params_rejects != NULL);
	m->params_step = gsl_vector_calloc(m->n_par);
	assert(m->params_step != NULL);
	m->params_min = gsl_vector_calloc(m->n_par);
	assert(m->params_min != NULL);
	m->params_max = gsl_vector_calloc(m->n_par);
	assert(m->params_max != NULL);

	m->params_descr = (const char**) mem_calloc(m->n_par, sizeof(char*));

	m->data = NULL;
	IFSEGV
		debug("allocating mcmc struct done");
	return m;
}

mcmc * mcmc_free(mcmc * m) {
	unsigned int i;

	mcmc_dump_close(m);
	gsl_rng_free(m->random);
	IFSEGV
		debug("freeing params");
	gsl_vector_free(m->params);
	IFSEGV
		debug("freeing params_best");
	gsl_vector_free(m->params_best);

	IFSEGV
		debug("freeing params_descr");
	for (i = 0; i < get_n_par(m); i++) {
		mem_free(m->params_descr[i]);
	}
	mem_free(m->params_descr);

	IFSEGV
		debug("freeing accepts/rejects");
	mem_free(m->params_accepts);
	mem_free(m->params_rejects);
	IFSEGV
		debug("freeing step/min/max");
	gsl_vector_free(m->params_step);
	gsl_vector_free(m->params_min);
	gsl_vector_free(m->params_max);
	if (m->data != NULL)
		gsl_matrix_free((gsl_matrix*) m->data);
	mem_free(m);
	m = NULL;
	return NULL;
}

void mcmc_check(const mcmc * m) {
	(void) m;
	assert(m != NULL);
	assert(m->n_par > 0);
	assert(m->data != NULL);
	assert(m->data->size2 > 0);
	assert(m->params != NULL);
	assert(m->params->size == m->n_par);
	assert(m->params_best != NULL);
	assert(m->params_best->size == m->n_par);
	assert(m->params_step != NULL);
}

