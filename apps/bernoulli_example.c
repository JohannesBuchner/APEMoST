#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

#define SIGMA 2

void calc_model(mcmc * m, const gsl_vector * old_values) {
	unsigned int i;
	unsigned int j;
	unsigned int n_par = get_n_par(m);
	double eta_i;
	double p_i;
	double l_i;
	double prior = 0;
	double prob = 0;

	(void) old_values;
	for (j = 1; j < m->data->size2; j++) {
		prior += -pow(get_params_for(m, j) / SIGMA, 2) / 2;
	}
	set_prior(m, prior);

	assert(n_par == m->data->size2 - 1 + 1);

	/* eta = [1,x] . params (matrix product) */
	for (i = 0; i < m->data->size1; i++) {
		eta_i = get_params_for(m, 0);
		for (j = 1; j < n_par; j++) {
			eta_i += gsl_matrix_get(m->data, i, j) * get_params_for(m, j);
		}

		if (eta_i > 0)
			p_i = 1 / (1 + exp(-eta_i));
		else
			p_i = exp(eta_i) / (1 + exp(eta_i));

		if (gsl_matrix_get(m->data, i, 0) == 0)
			l_i = log(1 - p_i);
		else
			l_i = log(p_i);
		prob += l_i;
	}

	set_prob(m, prior + get_beta(m) * prob);
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void) i;
	(void) old_value;

	calc_model(m, NULL);
}

