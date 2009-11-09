#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

#ifndef HMIN
#define HMIN 1e-6
#endif

static double calc_prior(const mcmc * m) {
	unsigned int i;
	const gsl_vector * params = m->params;
	double prior = 0;
	double mode_height;
	for (i = 3; i < get_n_par(m); i += 2) {
		mode_height = gsl_vector_get(params, i + 1);
		prior += gsl_sf_log(mode_height + HMIN);
	}
	i = (get_n_par(m) - 3) / 2;
	return -prior / i;
}

void calc_model(mcmc * m, const gsl_vector * old_values) {
	unsigned int i;
	unsigned int j;
	double y = 0;
	double mode_height;
	double mode_freq;
	double distance;
	double freq;
	const gsl_vector * params = m->params;
	double prob = gsl_vector_get(params, 1);
	double lifetime = gsl_vector_get(params, 0);
	double vrot = gsl_vector_get(params, 2);
	set_prior(m, calc_prior(m));

	(void) old_values;
	assert((get_n_par(m) - 3) % 2 == 0);
	for (i = 0; i < m->data->size1; i++) {
		y = 0;
		freq = gsl_matrix_get(m->data, i, 0);

		j = 3;
		mode_freq = gsl_vector_get(params, j);
		mode_height = gsl_vector_get(params, j + 1);
		distance = mode_freq - freq;
		y += mode_height / (1 + pow(2 * M_PI * distance * lifetime, 2));

		j = 5;
		mode_freq = gsl_vector_get(params, j);
		mode_height = gsl_vector_get(params, j + 1);
		distance = mode_freq - freq + -1 * vrot;
		y += mode_height / (1 + pow(2 * M_PI * distance * lifetime, 2));
		distance = mode_freq - freq;
		y += mode_height / (1 + pow(2 * M_PI * distance * lifetime, 2));
		distance = mode_freq - freq + 1 * vrot;
		y += mode_height / (1 + pow(2 * M_PI * distance * lifetime, 2));

		prob += gsl_sf_log(y) + gsl_matrix_get(m->data, i, 1) / y;
	}

	set_prob(m, get_prior(m) + -get_beta(m) * prob);
}

void calc_model_for(mcmc * m, const unsigned int j, const double old_value) {
	(void) j;
	(void) old_value;
	calc_model(m, NULL);
	return;
}

