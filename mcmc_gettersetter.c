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

const char ** get_params_descr(mcmc * m) {
	return m->params_descr;
}

gsl_histogram * get_hist(mcmc * m, int index, int nbins) {
	return calc_hist(m->params_distr[index], nbins);
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

double get_prob(mcmc * m) {
	return m->prob;
}

double get_prob_best(mcmc * m) {
	return m->prob_best;
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

int get_n_par(mcmc * m) {
	return m->n_par;
}

void set_params_best(mcmc * m, gsl_vector * new_params_best) {
	m->params_best = new_params_best;
}

void set_params_for(mcmc * m, double new_param, int i) {
	gsl_vector_set(m->params, i, new_param);
}

gsl_vector * get_params(mcmc * m) {
	return m->params;
}

gsl_vector * get_params_best(mcmc * m) {
	return m->params_best;
}

void set_params(mcmc * m, gsl_vector * new_params) {
	m->params = new_params;
}

void set_params_descr_all(mcmc * m, const char ** new_par_descr) {
	m->params_descr = new_par_descr;
}

void set_params_descr_for(mcmc * m, const char * new_par_descr, int i) {
	m->params_descr[i] = new_par_descr;
}

gsl_rng * get_random(mcmc * m) {
	return m->random;
}

void set_random(mcmc * m, gsl_rng * newrandom) {
	m->random = newrandom;
}

void set_prob(mcmc * m, double new_prob) {
	m->prob = new_prob;
}
gsl_vector * get_x(mcmc * m) {
	return m->x_dat;
}
gsl_vector * get_y(mcmc * m) {
	return m->y_dat;
}
void set_x(mcmc * m, gsl_vector * new_x) {
	m->x_dat = new_x;
}

void set_y(mcmc * m, gsl_vector * new_y) {
	m->y_dat = new_y;
}

void free_gsl_vector_array(gsl_vector ** arr) {
	int i = 0;
	if (arr != NULL) {
		while (arr[i] != NULL)
			gsl_vector_free(arr[i]);
	}
}

void set_steps_all(mcmc * m, double * new_steps) {
	unsigned int i;
	for (i = 1; i < m->n_par + 1; i++) {
		set_steps_for(m, new_steps[i + 1], i);
	}
}

