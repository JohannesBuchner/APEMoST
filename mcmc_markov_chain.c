#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"
#include <math.h>

void restart_from_best(mcmc * m) {
    set_params(m, dup_vector(get_params_best(m)));
    set_prob(m, -1E7);
}

void reset_accept_rejects(mcmc * m) {
	unsigned int i;
	for (i = 0; i < m->n_par; i++) {
		set_params_accepts_for(m, 0, i);
		set_params_rejects_for(m, 0, i);
	}
}

double abs_double(double x) {
	if(x < 0)
		return -x;
	else
		return x;
}

/**
 * aim a acceptance rate between 20 and 30%
 * @param rat_limit average acceptance rates for individual parameters to be achieved
 * @param burn_in_iterations number of burn-in iterations
 * @param iter_limit number of iterations for step width calibration
 * @param mul factor for adjusting the step width during calibration
 * @param adjust_step gives the factor with which to adjust the stepwidths after burn-in
 */
void markov_chain_calibrate(mcmc * m, unsigned int burn_in_iterations, double rat_limit,
		unsigned int iter_limit, double mul, double adjust_step) {
	unsigned int i;

	int flag;
	int cont = 0;
	gsl_vector * old_accepts;
	gsl_vector * old_rejects;
	gsl_vector * accept_rate;
	gsl_vector * reject_rate;
	gsl_vector * old_steps;
	double delta_reject_accept_t;

	unsigned long iter;
	unsigned long subiter;

	debug("Beginning calibration of MCMC ...");

	if(rat_limit < 0)
		rat_limit = pow(0.25, 1.0 / m->n_par);

	debug("Starting burn-in ...");
	for(iter = 0; iter < burn_in_iterations; iter++) {
		markov_chain_step(m);
		if(iter % 200 == 200 - 1) {
			dump_ul("\tBurn-in Iteration", iter);
		}
		mcmc_check_best(m);
		if (iter == burn_in_iterations/2) {
			debug("Re-initializing burn-in ...");
			restart_from_best(m);
		}
	}
	gsl_vector_scale(m->params_step, adjust_step);
	set_params(m, dup_vector(get_params_best(m)));
	debug("Burn-in done.");

	debug("Calibrating step widths ...(set cont=1 to abort)");
	reset_accept_rejects(m);

	while(iter < iter_limit) {
		for (i = 0; i < m->n_par; i++) {
			markov_chain_step(m); /* TODO:, i+1, calc_index=1); */
			mcmc_check_best(m);
		}
		if (iter % 200 == 200 - 1) {
			accept_rate = get_accept_rate(m);
			reject_rate = get_reject_rate(m);

			dump_v("steps", get_params(m));
			if(flag == 0) {
				flag = 1;
				old_accepts = accept_rate;
				old_rejects = get_reject_rate(m);
				old_steps = dup_vector(get_steps(m));
			}
			for(i = 0; i < m->n_par; i++) {
				dump_ul("----------------------------------------- iteration", iter);
				dump_s("Acceptance rates for", m->params_descr[i]);
				dump_v("acceptance rate: accepts", get_accept_rate(m));
				dump_v("acceptance rate: rejects", get_reject_rate(m));
				/* TODO: print diff to old_* variables */
			}
			old_accepts = accept_rate;
			old_rejects = reject_rate;
			old_steps = dup_vector(get_steps(m));

			for(i = 0; i < m->n_par; i++) {
				if(gsl_vector_get(accept_rate, i) > rat_limit + 0.05) {
					gsl_vector_scale(m->params_step, 1.0/mul);
				}
				if(gsl_vector_get(accept_rate, i) < rat_limit - 0.05) {
					gsl_vector_scale(m->params_step, mul);
				}
				dump_v("steps", m->params_step);
			}
			restart_from_best(m);
			reset_accept_rejects(m);
			for(subiter = 0; subiter < 200; subiter++) {
				markov_chain_step(m);
				mcmc_check_best(m);
			}
			dump_v("Global acceptance rate", get_accept_rate(m));
			delta_reject_accept_t = m->accept/(m->accept + m->reject) - 0.23;
			if (abs_double(delta_reject_accept_t) < 0.01) {
				cont = 1;
			} else {
				if (delta_reject_accept_t < 0) {
					rat_limit /= 0.99;
				} else {
					rat_limit *= 0.99;
				}
			}
		}
		iter++;
	}
}
void markov_chain_step(mcmc * m) {
	(void)m;

}
