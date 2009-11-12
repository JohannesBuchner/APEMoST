#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

void calc_model(mcmc * m, const gsl_vector * old_values) {
	double x = gsl_vector_get(m->params, 0);
	double a;
	double b = 0;
	double sigma;
	double pos;
	double height;
	unsigned int i;
	
	(void) old_values;
	
	for( i = 0; i < 10; i++) {
		pos = exp(i);
		height = 10*pow(1.0, i);
		sigma = i;
		if(i % 2 == 0)
			a = - sigma*pow((x - pos) / sigma, 2) / 2 + height;
		else
			if(x > pos)
				a = - height * (x - pos) / sigma + height;
			else
				a = - height * (pos - x) / sigma + height;
		if (a > b)
			b = a;
	}
	set_prob(m, get_beta(m) * b);
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void) i;
	(void) old_value;

	calc_model(m, NULL);
}

