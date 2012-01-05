#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"

struct Problem {
	double (*Prior)(mcmc * m, const gsl_vector * old_values);
	double (*LogLike)(mcmc * m, const gsl_vector * old_values);
} p;

void set_function(
	double (*LogLike)(mcmc * m, const gsl_vector * old_values),
	double (*Prior)(mcmc * m, const gsl_vector * old_values)
) {
	p.LogLike = LogLike;
	p.Prior = Prior;
}

void calc_model(mcmc * m, const gsl_vector * old_values) {
	set_prior(m, p.Prior(m, old_values));
	set_prob(m, get_beta(m) * (get_prior(m) + p.LogLike(m,old_values)));
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void) i;
	(void) old_value;

	calc_model(m, NULL);
}

