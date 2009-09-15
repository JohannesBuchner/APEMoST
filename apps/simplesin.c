#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

#ifndef SIGMA
#define SIGMA 0.5
#endif

double apply_formula(mcmc * m, unsigned int i, double amplitude, double frequency,
		double phase) {
	double x = gsl_matrix_get(m->data, i, 0);
	double y = amplitude * gsl_sf_sin(2.0 * M_PI * frequency * x + phase);
	return y;
}

void calc_model(mcmc * m, const gsl_vector * old_values) {
	unsigned int i;
	double amplitude = gsl_vector_get(m->params, 0);
	double frequency = gsl_vector_get(m->params, 1);
	double phase     = gsl_vector_get(m->params, 2);
	double y;
	double deltay;
	double square_sum = 0;

	(void) old_values;
	/*dump_v("recalculating model for parameter values", m->params);*/
	for (i = 0; i < m->data->size1; i++) {
		y = gsl_matrix_get(m->data, i, 1);
		deltay = apply_formula(m, i, amplitude, frequency, phase) - y;
		square_sum += deltay * deltay;
	}
	set_prob(m, get_beta(m) * square_sum / (-2 * SIGMA * SIGMA));
	/*debug("model done");*/
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void) i;
	(void) old_value;

	calc_model(m, NULL);
}

