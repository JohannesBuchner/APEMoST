/*
    APEMoST - Automated Parameter Estimation and Model Selection Toolkit
    Copyright (C) 2009  Johannes Buchner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_sf.h>

unsigned long get_params_accepts_sum(const mcmc * m) {
	unsigned int i;
	unsigned long sum = 0;
	for (i = 0; i < get_n_par(m); i++) {
		sum += m->params_accepts[i];
	}
	return sum;
}

unsigned long get_params_rejects_sum(const mcmc * m) {
	unsigned int i;
	unsigned long sum = 0;
	for (i = 0; i < get_n_par(m); i++) {
		sum += m->params_rejects[i];
	}
	return sum;
}

unsigned long get_params_accepts_for(const mcmc * m, const unsigned int i) {
	return m->params_accepts[i];
}

unsigned long get_params_rejects_for(const mcmc * m, const unsigned int i) {
	return m->params_rejects[i];
}
double get_accept_rate_for(const mcmc * m, const unsigned int i) {
	return get_params_accepts_for(m, i) / (double) (
			get_params_accepts_for(m, i) + get_params_rejects_for(m, i));
}

unsigned long get_params_accepts_global(const mcmc * m) {
	return m->accept;
}

unsigned long get_params_rejects_global(const mcmc * m) {
	return m->reject;
}
double get_accept_rate_global(const mcmc * m) {
	return get_params_accepts_global(m) / (double)
			(get_params_accepts_global(m) + get_params_rejects_global(m));
}

gsl_vector * get_vector_from_array(const unsigned long * array, const unsigned int size)
{
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

void set_prob_best(mcmc * m, const double new_prob_best) {
	m->prob_best = new_prob_best;
}

double get_prob(const mcmc * m) {
	return m->prob;
}

double get_prior(const mcmc * m) {
	return m->prior;
}

double get_prob_best(const mcmc * m) {
	return m->prob_best;
}

void set_minmax_for(mcmc * m, const double new_min, const double new_max,
		const unsigned int i) {
	gsl_vector_set(m->params_min, i, new_min);
	gsl_vector_set(m->params_max, i, new_max);
}

void set_steps_for(mcmc * m, const double new_step, const unsigned int i) {
	gsl_vector_set(m->params_step, i, new_step);
}

void set_steps_for_normalized(mcmc * m, const double new_step, const unsigned int i) {
	gsl_vector_set(m->params_step, i, new_step *
		(get_params_max_for(m, i) - get_params_min_for(m, i)));
}

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Fix the number of parameters to the given value
 *
 * Setting this might leads the compiler to more optimization.
 * If in doubt, do not set, as your program then can be used for any
 * problems regardless of number of parameters.
 *
 * my favorite trick: automatically count the lines in the current file
 * CCFLAGS="-DN_PARAMETERS="$(cat params|wc -l)"
 */
#define N_PARAMETERS
#endif

#ifdef N_PARAMETERS
#define get_n_par(m) N_PARAMETERS
#else
unsigned int get_n_par(const mcmc * m) {
	return m->n_par;
}
#endif

void set_params_best(mcmc * m, const gsl_vector * new_params_best) {
	gsl_vector_memcpy(m->params_best, new_params_best);
}

void set_params_for(mcmc * m, const double new_param, const unsigned int i) {
	assert(i < m->n_par);
	gsl_vector_set(m->params, i, new_param);
}

gsl_vector * get_params(const mcmc * m) {
	return m->params;
}

double get_params_for(const mcmc * m, const unsigned int i) {
	return gsl_vector_get(m->params, i);
}

gsl_vector * get_params_min(const mcmc * m) {
	return m->params_min;
}

double get_params_min_for(const mcmc * m, const unsigned int i) {
	return gsl_vector_get(m->params_min, i);
}
gsl_vector * get_params_max(const mcmc * m) {
	return m->params_max;
}
double get_params_max_for(const mcmc * m, const unsigned int i) {
	return gsl_vector_get(m->params_max, i);
}

gsl_vector * get_params_best(const mcmc * m) {
	return m->params_best;
}
double get_params_best_for(const mcmc * m, const unsigned int i) {
	return gsl_vector_get(m->params_best, i);
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
void set_prior(mcmc * m, const double new_prior) {
	m->prior = new_prior;
}
const gsl_matrix * get_data(const mcmc * m) {
	return m->data;
}
void set_data(mcmc * m, const gsl_matrix * new_data) {
	m->data = new_data;
}

gsl_vector * get_steps(const mcmc * m) {
	return m->params_step;
}

double get_steps_for(const mcmc * m, const unsigned int i) {
	return gsl_vector_get(m->params_step, i);
}
double get_steps_for_normalized(const mcmc * m, const unsigned int i) {
	return get_steps_for(m, i) / (get_params_max_for(m, i) -
			get_params_min_for(m, i));
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

double get_next_random_jump(const mcmc * m, const double sigma) {
#ifdef PROPOSAL_LOGISTIC
	return gsl_ran_logistic(get_random(m), sigma);
#elif defined PROPOSAL_UNIFORM
	return gsl_ran_flat(get_random(m), -sigma, sigma);
#else
/**
 * You can choose a proposal distribution. The default is a gaussian 
 * distribution.
 *
 * Set PROPOSAL_LOGISTIC, PROPOSAL_UNIFORM if you want to use a different
 * proposal distribution.
 */
#define PROPOSAL
	return gsl_ran_gaussian(get_random(m), sigma);
#endif
}
double get_next_alog_urandom(const mcmc * m) {
	return gsl_sf_log(get_next_uniform_random(m));
}

