#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>

double get_random_number() {
	gsl_rng * random;
	double v;
	
	gsl_rng_env_setup();
	random = gsl_rng_alloc(gsl_rng_default);
	v = gsl_rng_uniform(random);
	gsl_rng_free(random);
	
	return v;
}

void mcmc_init(mcmc m) {
    m.n_par = 0;
    m.accept = 0;
    m.reject = 0;
    m.prob = 0;
    m.prob_best = 0;
    m.seed = NULL;
    m.params = NULL;
    m.params_best = NULL;
    m.params_distr = NULL;
    m.params_distr_buf = NULL;
    m.params_descr = NULL;
    m.params_priors = NULL;
    m.params_ar = NULL;
    m.params_step = NULL;
    m.params_minmax = NULL;
    m.x_dat = NULL;
    m.y_dat = NULL;
    m.model = NULL;
    
    get_random_number();
    m.prob = -1e+10;
    m.prob_best = -1e+10;
}

/* TODO: cleanup */

/* TODO: on setting the parameter length, preallocate all 2D arrays. It is a 
 *       slight coding issue when the 2D-arrays remain NULL-pointers.
 */

/* TODO: change from i-1 addressing to i */

/* TODO: all setters seem to need copy functions, check if we really need those
 *       This can probably be done by the caller (using gsl_vector_alloc() and
 *       gsl_vector_memcpy())
 */

/* TODO: check if we need getters. so far they look straight-forward, so 
 *       I'd let people access the struct directly.
 */

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

void setup(mcmc * m, const char * filename) {
	(void)m;
	(void)filename;
	/* TODO */
}

gsl_vector * get_params_ar(mcmc * m) {
	unsigned int i;
	double sum = 0;
	gsl_vector * r = gsl_vector_alloc(m->n_par);
	for (i = 0; i < m->n_par; i++) {
		sum += calc_vector_sum(m->params_ar[i]);
		gsl_vector_set(r, i, gsl_vector_get(m->params_ar[i], 0));
	}
	gsl_vector_scale(r, 1/sum);
	return r;
}

double get_params_ar_for(mcmc * m, int i) {
	return gsl_vector_get(m->params_ar[i], 0) / calc_vector_sum(m->params_ar[i]);
}

gsl_histogram * get_hist(mcmc * m, int index, int nbins) {
	return calc_hist(m->params_distr, index, nbins);
}

void set_params_ar(mcmc * m, gsl_vector ** new_params_ar) {
	m->params_ar = new_params_ar;
}

void set_params_ar_for(mcmc * m, gsl_vector * new_params_ar, int i) {
	m->params_ar[i] = new_params_ar;
}

void set_prob_best(mcmc * m, double new_prob_best) {
	m->prob_best = new_prob_best;
}

void set_minmax_for(mcmc * m, gsl_vector * new_minmax, int i) {
	m->params_minmax[i] = new_minmax;
}

void set_minmax(mcmc * m, gsl_vector ** new_minmax) {
	m->params_minmax = new_minmax;
}

void set_model(mcmc * m, gsl_vector * new_model) {
	m->model = new_model;
}

void set_n_par(mcmc * m, int new_n_par) {
	m->n_par = new_n_par;
}

void set_params_best(mcmc * m, gsl_vector * new_params_best) {
	m->params_best = new_params_best;
}

void set_params_for(mcmc * m, double new_param, int i) {
	gsl_vector_set(m->params, i, new_param);
}

void set_params(mcmc * m, gsl_vector * new_params) {
	m->params = new_params;
}

void set_params_descr_all(mcmc * m, char ** new_par_descr) {
	m->params_descr = new_par_descr;
}

void set_params_descr_for(mcmc * m, char * new_par_descr, int i) {
	m->params_descr[i-1] = new_par_descr;
}

void set_seed(mcmc * m, gsl_vector * new_seed) {
	m->seed = new_seed;
}

void set_probability(mcmc * m, double new_prob) {
	m->prob = new_prob;
}

void set_x(mcmc * m, gsl_vector * new_x) {
	/*if(m->x_dat != NULL)
		gsl_vector_free(m->x_dat);*/
	m->x_dat = new_x;
}

/* TODO: check if we need that
void set_x_copy(mcmc * m, gsl_vector * new_x) {
	gsl_vector_memcpy(m->x_dat, new_x);
}*/

void set_y(mcmc * m, gsl_vector * new_y) {
	/*if(m->x_dat != NULL)
		gsl_vector_free(m->x_dat);*/
	m->y_dat = new_y;
}

/* TODO: check if we need that
void set_y_copy(mcmc * m, gsl_vector * new_y) {
	gsl_vector_memcpy(m->y_dat, new_y);
}*/

void set_steps_for(mcmc * m, double new_steps, int i) {
	m->params_step[i-1] = new_steps;
}

void free_gsl_vector_array(gsl_vector ** arr) {
	int i = 0;
	if(arr != NULL) {
		while(arr[i]!=NULL)
			gsl_vector_free(arr[i]);
	}
}

/*
gsl_vector ** copy_gsl_vector_array(gsl_vector ** arr, const gsl_vector ** src, size_t n) {
	int i = 0;
	assert(src != NULL);
	if(arr == NULL) {
		arr = (gsl_vector *)calloc(n, sizeof(gsl_vector *));
		assert(arr != NULL);
		for(i = 0; i < n; i++) {
			arr[i] = gsl_vector_alloc(n);
			gsl_vector_memcpy(arr[i], src[i]);
		}
		return arr;
	}else{
		for(i = 0; i < n; i++) {
			gsl_vector_memcpy(arr[i], src[i]);
		}
	}
}*/

void set_steps_all(mcmc * m, double * new_steps) {
	unsigned int i;
	for(i = 1; i < m->n_par + 1; i++) {
		set_steps_for(m, new_steps[i+1], i);
	}
}

