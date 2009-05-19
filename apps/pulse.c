#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

#ifndef HMIN
#define HMIN 1e-6
#endif

double apply_formula(const mcmc * m, const unsigned int j) {
	unsigned int i = 0;
	double freq = gsl_vector_get(m->x_dat, j);
	/* i = 0 is lifetime. i is frequency, i+1 is height. */
	double lifetime = gsl_vector_get(get_params(m), 0);
	double y = 0;
	double mode_height;
	double mode_freq;
	double distance;
	for (i = 1; i < get_n_par(m); i += 2) {
		mode_height = gsl_vector_get(get_params(m), i);
		mode_freq = gsl_vector_get(get_params(m), i + 1);
		distance = mode_freq - freq;
		y += mode_height / (1 + 4* distance * distance * M_PI * lifetime);
	}
	gsl_vector_set(m->model, j, y);
	return y;
}

double calc_prior(const mcmc * m) {
	unsigned int i;
	double prior = 0;
	double height;
	/* i = 0 is lifetime. i is frequency, i+1 is height. */
	for (i = 1; i < get_n_par(m); i += 2) {
		height = gsl_vector_get(get_params(m), i + 1);
		prior += gsl_sf_log(height + HMIN);
	}
	return prior;
}

void calc_model(mcmc * m, const gsl_vector * old_values) {
	unsigned int i;
	double prior = calc_prior(m);
	double prob = 0;
	double power_model;
	(void) old_values;

	for (i = 0; i < m->x_dat->size; i++) {
		power_model = apply_formula(m, i);
		prob += gsl_sf_log(power_model) + gsl_vector_get(m->y_dat, i)
				/ power_model;
	}

	set_prob(m, get_beta(m) * -1 * (prob + prior));
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void) i;
	(void) old_value;

	calc_model(m, NULL);
}

