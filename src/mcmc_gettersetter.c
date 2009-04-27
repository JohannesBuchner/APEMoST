#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

long get_params_accepts_sum(const mcmc * m) {
	unsigned int i;
	long sum = 0;
	for (i = 0; i < get_n_par(m); i++) {
		sum += m->params_accepts[i];
	}
	return sum;
}

long get_params_rejects_sum(const mcmc * m) {
	unsigned int i;
	long sum = 0;
	for (i = 0; i < get_n_par(m); i++) {
		sum += m->params_rejects[i];
	}
	return sum;
}

long get_params_accepts_for(const mcmc * m, const int i) {
	return m->params_accepts[i];
}

long get_params_rejects_for(const mcmc * m, const int i) {
	return m->params_rejects[i];
}

gsl_vector * get_vector_from_array(const long * array, const unsigned int size) {
	unsigned int i;
	gsl_vector * v = gsl_vector_alloc(size);
	assert(v != NULL);
	for (i = 0; i < size; i++) {
		gsl_vector_set(v, i, array[i]);
	}
	return v;
}

gsl_vector * get_accept_rate(const mcmc * m) {
	gsl_vector * a = get_vector_from_array(m->params_accepts, get_n_par(m));
	gsl_vector * r = get_vector_from_array(m->params_rejects, get_n_par(m));
	gsl_vector * sum = r;
	gsl_vector_add(sum, a);
	gsl_vector_div(a, sum);
	gsl_vector_free(sum);
	return a;
}

const char ** get_params_descr(const mcmc * m) {
	return m->params_descr;
}

gsl_histogram * get_hist(const mcmc * m, int index, int nbins) {
	return calc_hist(m->params_distr[index], nbins);
}

void inc_params_accepts_for(mcmc * m, const unsigned int i) {
	m->params_accepts[i]++;
}
void inc_params_rejects_for(mcmc * m, const unsigned int i) {
	m->params_rejects[i]++;
}
void inc_params_accepts(mcmc * m) {
	unsigned int i;
	m->accept++;
	for (i = 0; i < get_n_par(m); i++)
		inc_params_accepts_for(m, i);
}
void inc_params_rejects(mcmc * m) {
	unsigned int i;
	m->reject++;
	for (i = 0; i < get_n_par(m); i++)
		inc_params_rejects_for(m, i);
}

void set_params_accepts_for(mcmc * m, const long new_params_accept,
		const unsigned int i) {
	m->params_accepts[i] = new_params_accept;
}
void set_params_rejects_for(mcmc * m, const long new_params_reject,
		const unsigned int i) {
	m->params_rejects[i] = new_params_reject;
}
void reset_accept_rejects(mcmc * m) {
	unsigned int i;
	for (i = 0; i < get_n_par(m); i++) {
		set_params_accepts_for(m, 0, i);
		set_params_rejects_for(m, 0, i);
	}
	m->reject = 0;
	m->accept = 0;
}

void set_prob_best(mcmc * m, double new_prob_best) {
	m->prob_best = new_prob_best;
}

double get_prob(const mcmc * m) {
	return m->prob;
}

double get_prob_best(const mcmc * m) {
	return m->prob_best;
}

void set_minmax_for(mcmc * m, const double new_min, const double new_max,
		const unsigned int i) {
	gsl_vector_set(m->params_min, i, new_min);
	gsl_vector_set(m->params_max, i, new_max);
}

void set_steps_for(mcmc * m, const double new_step, const int i) {
	gsl_vector_set(m->params_step, i, new_step);
}

void set_model(mcmc * m, gsl_vector * new_model) {
	gsl_vector_free(m->model);
	m->model = new_model;
}

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Fix the number of parameters to the given value
 *
 * Setting this might lead the compiler to more optimization, a theory that
 * has not been proven for the program.
 * If in doubt, do not set, as your program then can be used for any
 * problems regardless of number of parameters.
 */
#define N_PARAMETERS
#endif

unsigned int get_n_par(const mcmc * m) {
#ifdef N_PARAMETERS
	(void)m;
	return N_PARAMETERS;
#else
#define N_PARAMETERS
	return m->n_par;
#endif
}

void set_params_best(mcmc * m, gsl_vector * new_params_best) {
	m->params_best = new_params_best;
}

void set_params_for(mcmc * m, const double new_param, const unsigned int i) {
	assert(i < m->n_par);
	gsl_vector_set(m->params, i, new_param);
}

gsl_vector * get_params(const mcmc * m) {
	return m->params;
}

gsl_vector * get_params_best(const mcmc * m) {
	return m->params_best;
}

void set_params(mcmc * m, gsl_vector * new_params) {
	assert(m->n_par == new_params->size);
	gsl_vector_free(m->params);
	m->params = new_params;
}

void set_params_descr_all(mcmc * m, const char ** new_par_descr) {
	m->params_descr = new_par_descr;
}

void set_params_descr_for(mcmc * m, const char * new_par_descr,
		const unsigned int i) {
	m->params_descr[i] = new_par_descr;
}

gsl_rng * get_random(const mcmc * m) {
	return m->random;
}

void set_random(mcmc * m, gsl_rng * newrandom) {
	m->random = newrandom;
}

void set_prob(mcmc * m, const double new_prob) {
	m->prob = new_prob;
}
const gsl_vector * get_x(const mcmc * m) {
	return m->x_dat;
}
const gsl_vector * get_y(const mcmc * m) {
	return m->y_dat;
}
void set_x(mcmc * m, const gsl_vector * new_x) {
	m->x_dat = new_x;
}

void set_y(mcmc * m, const gsl_vector * new_y) {
	m->y_dat = new_y;
}

gsl_vector * get_steps(const mcmc * m) {
	return m->params_step;
}

void free_gsl_vector_array(gsl_vector ** arr) {
	int i = 0;
	if (arr != NULL) {
		while (arr[i] != NULL)
			gsl_vector_free(arr[i]);
	}
}

void set_steps_all(mcmc * m, const double * new_steps) {
	unsigned int i;
	for (i = 0; i < get_n_par(m); i++) {
		set_steps_for(m, new_steps[i], i);
	}
}

double get_next_uniform_plusminus_random(const mcmc * m) {
	return 2* get_next_uniform_random (m) - 1;
}
double get_next_uniform_random(const mcmc * m) {
	return gsl_rng_uniform(get_random(m));
}
double get_next_gauss_random(const mcmc * m, double sigma) {
	return gsl_ran_gaussian(get_random(m), sigma);
}
double get_next_alog_urandom(const mcmc * m) {
	return gsl_sf_log(get_next_uniform_random(m));
}
