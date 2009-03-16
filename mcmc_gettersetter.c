#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>

/* TODO: check if we need getters. so far they look straight-forward, so 
 *       I'd let people access the struct directly.
 */

long get_params_accepts(mcmc * m) {
	unsigned int i;
	long sum = 0;
	for (i = 0; i < m->n_par; i++) {
		sum += m->params_accepts[i];
	}
	return sum;
}
long get_params_rejects(mcmc * m) {
	unsigned int i;
	long sum = 0;
	for (i = 0; i < m->n_par; i++) {
		sum += m->params_rejects[i];
	}
	return sum;
}
long get_params_accepts_for(mcmc * m, int i) {
	return m->params_accepts[i];
}
long get_params_rejects_for(mcmc * m, int i) {
	return m->params_rejects[i];
}

gsl_histogram * get_hist(mcmc * m, int index, int nbins) {
	return calc_hist(m->params_distr, index, nbins);
}


void set_params_accepts_for(mcmc * m, long new_params_accept, int i) {
	m->params_accepts[i] = new_params_accept;
}
void set_params_rejects_for(mcmc * m, long new_params_reject, int i) {
	m->params_rejects[i] = new_params_reject;
}

void set_prob_best(mcmc * m, double new_prob_best) {
	m->prob_best = new_prob_best;
}

void set_minmax_for(mcmc * m, double new_min, double new_max, int i) {
	gsl_vector_set(m->params_min, i, new_min);
	gsl_vector_set(m->params_max, i, new_max);
}

void set_steps_for(mcmc * m, double new_step, int i) {
	gsl_vector_set(m->params_step, i, new_step);
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
	m->params_descr[i] = new_par_descr;
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

void free_gsl_vector_array(gsl_vector ** arr) {
	int i = 0;
	if(arr != NULL) {
		while(arr[i]!=NULL)
			gsl_vector_free(arr[i]);
	}
}

/*
gsl_vector * alloc_gsl_vector_array(unsigned int size) {
	int i = 0;
	arr = (gsl_vector *)calloc(size, sizeof(gsl_vector *));
	assert(arr != NULL);
	for(i = 0; i < size; i++) {
		arr[i] = gsl_vector_alloc(n);
		gsl_vector_memcpy(arr[i], src[i]);
	}
	return arr;
}
*/
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

