#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include "debug.h"

#define free(p) { IFSEGV dump_p("about to free", (void*)p); free(p); }

double get_random_number() {
	gsl_rng * random;
	double v;
	
	gsl_rng_env_setup();
	random = gsl_rng_alloc(gsl_rng_default);
	v = gsl_rng_uniform(random);
	gsl_rng_free(random);
	
	return v;
}

#define ALLOCATION_CHUNKS 1
/**
 * for allocation, we don't want to call alloc too often, rather grow in steps
 */
unsigned long get_new_size(unsigned long new_iter) {
	return ((new_iter+1)/ALLOCATION_CHUNKS + 1)*ALLOCATION_CHUNKS;
}

/**
 * make space available as set in m->size
 */
void resize(mcmc * m, unsigned long new_size) {
	unsigned long orig_size = m->size;
	unsigned long i;
	/* TODO: we can not allocate more than the int-space
	 * if we have more iterations than that, we need to move our data to the
	 * disk.
	 */
	dump_ul("resizing params_distr to", new_size);
	if(new_size < orig_size) {
		debug("shrinking -> freeing vectors");
		for(i = orig_size;i > new_size; i--) {
			dump_ul("freeing vector", i - 1);
			gsl_vector_free(m->params_distr[i - 1]);
		}
	}
	m->params_distr = realloc(m->params_distr, new_size);
	
	if(new_size > orig_size) {
		debug("growing -> allocating vectors");
		for(i = orig_size; i < new_size; i++) {
			dump_ul("allocating vector", i);
			m->params_distr[i] = gsl_vector_alloc(m->n_par);
		}
	}
	m->size = new_size;
	debug("done resizing");
}

void prepare_iter(mcmc * m, unsigned long iter) {
	unsigned long new_size = get_new_size(iter);
	if(m->size != new_size) {
		resize(m, new_size);
	}
}

void init_seed(mcmc * m) {
	/* TODO: find right type and initialization for random generator */
	/* TODO: fill arrays with initial values */
	get_random_number();
	m->seed = NULL; 
}

mcmc * mcmc_init(unsigned int n_pars) {
	mcmc * m;
	IFSEGV debug("allocating mcmc struct");
	m = (mcmc*)malloc(sizeof(mcmc));
	assert(m != NULL);
	m->iter = 0;
	m->size = 0;
	m->n_par = n_pars;
	m->accept = 0;
	m->reject = 0;
	m->prob = -1e+10;
	m->prob_best = -1e+10;
	
	init_seed(m);
	
	m->params = gsl_vector_alloc(m->n_par);
	assert(m->params != NULL);
	m->params_best = gsl_vector_alloc(m->n_par);
	assert(m->params_best != NULL);
	m->params_distr = NULL;
	prepare_iter(m, 0);
	assert(m->params_distr != NULL);
	
	m->params_accepts = (long*)calloc(m->n_par, sizeof(long));
	assert(m->params_accepts != NULL);
	m->params_rejects = (long*)calloc(m->n_par, sizeof(long));
	assert(m->params_rejects != NULL);
	m->params_step   = gsl_vector_calloc(m->n_par);
	assert(m->params_step != NULL);
	m->params_min    = gsl_vector_calloc(m->n_par);
	assert(m->params_min != NULL);
	m->params_max    = gsl_vector_calloc(m->n_par);
	assert(m->params_max != NULL);

	m->params_descr = (char**)calloc(m->n_par, sizeof(char*));;
	m->x_dat = NULL;
	m->y_dat = NULL;
	m->model = NULL;
	IFSEGV debug("allocating mcmc struct done");
	return m;
}

void mcmc_free(mcmc * m) {
	unsigned int i;
	
	m->seed = NULL; /* TODO */
	IFSEGV debug("freeing params");
	gsl_vector_free(m->params);
	IFSEGV debug("freeing params_best");
	gsl_vector_free(m->params_best);
	IFSEGV debug("freeing params_distr");
	resize(m, 0);
	
	IFSEGV debug("freeing params_descr");
	for (i = 0; i < m->n_par; i++) {
		free(m->params_descr[i]);
	}
	free(m->params_descr);
	
	IFSEGV debug("freeing accepts/rejects");
	free(m->params_accepts);
	free(m->params_rejects);
	IFSEGV debug("freeing step/min/max");
	gsl_vector_free(m->params_step);
	gsl_vector_free(m->params_min);
	gsl_vector_free(m->params_max);
	free(m->x_dat);
	free(m->y_dat);
	free(m->model);
	free(m);
}


extern void calc_model(mcmc * m);
extern void calc_model_for(mcmc * m, int i);

extern double calc_prob(mcmc * m);

void add_values(mcmc * m, int n_iter) {
	/* TODO */
	(void)m;
	n_iter += 0;
}

void write2files(mcmc * m) {
	(void)m;
	/* TODO */
}

