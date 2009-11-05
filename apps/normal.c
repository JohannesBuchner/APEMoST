#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

#define APOS 0.1
#define BPOS 0.2
#define CPOS -9

void calc_model(mcmc * m, const gsl_vector * old_values) {
	double x = gsl_vector_get(m->params, 0);
	double a;
	double b;
	double c;

	(void) old_values;
	*prior = 0;
	a = -pow( (x-APOS)*400, 2)/2;
	b = -pow( (x-BPOS)*400, 2)/2;
	c = -pow( (x-CPOS)*400, 2)/2;

	if(a > b && a > c)
		set_prob(m, get_beta(m) * a);
	else if(b > c)
		set_prob(m, get_beta(m) * b);
	else
		set_prob(m, get_beta(m) * c);
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void) i;
	(void) old_value;

	calc_model(m, NULL);
}

