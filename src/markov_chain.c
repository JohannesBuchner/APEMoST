#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"
#include <math.h>
#include <gsl/gsl_sf.h>

void restart_from_best(mcmc * m) {
	set_params(m, dup_vector(get_params_best(m)));
	set_prob(m, -1E7);
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

double abs_double(double x) {
	if (x < 0)
		return -x;
	else
		return x;
}

void markov_chain_calibrate(mcmc * m, unsigned int burn_in_iterations,
		double rat_limit, unsigned int iter_limit, double mul,
		double adjust_step) {
	/* we aim a acceptance rate between 20 and 30% */
	unsigned int i;

	int reached_perfection = 0;
	gsl_vector * old_accepts = NULL;
	gsl_vector * accept_rate = NULL;
	gsl_vector * old_steps = NULL;
	double delta_reject_accept_t;

	unsigned long iter = 0;
	unsigned long subiter;
	int nchecks_without_rescaling = 0;
	int rescaled;

	if (rat_limit < 0)
		rat_limit = pow(0.25, 1.0 / get_n_par(m));

	debug("Beginning calibration of MCMC ...");
	debug("Starting burn-in ...");
	mcmc_check(m);
	for (iter = 0; iter < burn_in_iterations / 2; iter++) {
		for (subiter = 0; subiter < 200; subiter++) {
			markov_chain_step(m, 0);
		}
		iter += subiter;
		dump_ul("\tBurn-in Iteration", iter);
		mcmc_check_best(m);
	}
	debug("Re-initializing burn-in ...");
	restart_from_best(m);
	for (; iter < burn_in_iterations; iter++) {
		for (subiter = 0; subiter < 200; subiter++) {
			markov_chain_step(m, 0);
		}
		iter += subiter;
		dump_ul("\tBurn-in Iteration", iter);
		mcmc_check_best(m);
	}
	debug("Burn-in done, adjusting steps ...");
	mcmc_check(m);
	gsl_vector_scale(m->params_step, adjust_step);
	set_params(m, dup_vector(get_params_best(m)));
	debug("Burn-in done.");

	debug("Calibrating step widths ...(set cont=1 to abort)");
	reset_accept_rejects(m);

	while (iter < iter_limit) {
		for (i = 0; i < get_n_par(m); i++) {
			markov_chain_step_for(m, i, 1);
			mcmc_check_best(m);
		}
		iter++;
		if (iter % 200 == 0) {
			accept_rate = get_accept_rate(m);

			if (old_accepts == NULL) {
				old_accepts = accept_rate;
				old_steps = dup_vector(get_steps(m));
			}
			dump_ul("------------------------------------------------ iteration", iter);
			dump_v("params", get_params(m));
			dump_v("acceptance rate: ", get_accept_rate(m));
			dump_v("steps", get_steps(m));

			for (i = 0; i < get_n_par(m); i++) {
				/*dump_i_s("Acceptance rates for", i, m->params_descr[i]);*/
				/* TODO: print diff to old_* variables */
			}
			/*gsl_vector_free(accept_rate);*/
			old_accepts = accept_rate;
			/*gsl_vector_free(old_steps);*/
			old_steps = dup_vector(get_steps(m));

			rescaled = 0;
			for (i = 0; i < get_n_par(m); i++) {
				IFDEBUG
					printf(
							"\t\tneeded acceptance rate: <%f, >%f; got %f for %i",
							rat_limit + 0.05, rat_limit - 0.05, gsl_vector_get(
									accept_rate, i), i);
				if (gsl_vector_get(accept_rate, i) > rat_limit + 0.05) {
					gsl_vector_set(m->params_step, i, gsl_vector_get(
							m->params_step, i) / mul);
					IFDEBUG
						printf("\t scaling up   ^");
					rescaled = 1;
				}
				if (gsl_vector_get(accept_rate, i) < rat_limit - 0.05) {
					gsl_vector_set(m->params_step, i, gsl_vector_get(
							m->params_step, i) * mul);
					IFDEBUG
						printf("\t scaling down v");
					rescaled = 1;
				}
				IFDEBUG
					printf("\n");
				dump_v("steps", m->params_step);
			}
			if (rescaled == 0)
				nchecks_without_rescaling++;
			restart_from_best(m);
			reset_accept_rejects(m);
			for (subiter = 0; subiter < 200; subiter++) {
				markov_chain_step(m, 0);
				mcmc_check_best(m);
			}
			dump_v("New overall accept rate after reset", get_accept_rate(m));
			delta_reject_accept_t = m->accept * 1.0 / (m->accept + m->reject)
					- 0.23;
			dump_d("Compared to desired rate", delta_reject_accept_t);
			if (abs_double(delta_reject_accept_t) < 0.01) {
				reached_perfection = 1;
				debug("calibration reached the desired acceptance rate");
			} else {
				reached_perfection = 0;
				if (delta_reject_accept_t < 0) {
					rat_limit /= 0.99;
				} else {
					rat_limit *= 0.99;
				}
			}
			if (nchecks_without_rescaling > NO_RESCALING_LIMIT
					&& reached_perfection == 1 && rescaled == 0) {
				debug("quitting calibration because we did not need to rescale for several times");
				break;
			}
		}
	}
	reset_accept_rejects(m);
	debug("calibration of markov-chain done.");
}

double mod_double(double x, double div) {
	while (1) {
		if (x >= div)
			x -= div;
		else if (x < 0)
			x += div;
		else
			return x;
	}
}

void do_step_for(mcmc * m, unsigned int i) {
	double step = gsl_vector_get(m->params_step, i);
	double old_value = gsl_vector_get(m->params, i);
	double new_value = old_value + get_next_gauss_random(m, step);
	double max = gsl_vector_get(m->params_max, i);
	double min = gsl_vector_get(m->params_min, i);
	/* dump_d("Jumping from", old_value); */

	if (new_value > max)
		new_value = max - mod_double(new_value - max, max - min);
	else if (new_value < min)
		new_value = min + mod_double(min - new_value, max - min);
	/* dump_d("To", new_value); */
	set_params_for(m, new_value, i);
}

void do_step(mcmc * m) {
	unsigned int i;
	for (i = 0; i < get_n_par(m); i++) {
		do_step_for(m, i);
	}
}

/**
 * @returns 1 if accept, 0 if rejecting
 */
int check_accept(mcmc * m, double prob_old) {
	double prob_new = get_prob(m);
	double prob_still_accept;

	/* shortcut */
	if (prob_new == prob_old) {
		return 1;
	}

	if (prob_new > prob_old) {
		IFVERBOSE
			dump_d("accepting improvement of", prob_new - prob_old);
		return 1;
	} else {
		prob_still_accept = get_next_alog_urandom(m);
		if (prob_still_accept < (prob_new - prob_old)) {
			IFVERBOSE {
				dump_d("accepting probability", prob_still_accept);
				dump_d("accepting worsening of", prob_new - prob_old);
			}
			return 1;
		} else {
			IFVERBOSE
				dump_d("rejecting worsening of", prob_new - prob_old);
			return 0;
		}
	}
}

void revert(mcmc * m, gsl_vector * old_model, double prob_old) {
	set_model(m, old_model);
	set_prob(m, prob_old);
}

void markov_chain_step_for(mcmc * m, unsigned int index, int calc_index) {
	double prob_old = get_prob(m);
	gsl_vector * old_model = dup_vector(m->model);
	double old_value = gsl_vector_get(m->params, index);

	mcmc_check(m);
	do_step_for(m, index);

	if (calc_index <= 1) {
		calc_model_for(m, index, old_value);
	}

	if (check_accept(m, prob_old) == 1) {
		inc_params_accepts_for(m, index);
		gsl_vector_free(old_model);
	} else {
		revert(m, old_model, prob_old);
		set_params_for(m, old_value, index);
		inc_params_rejects_for(m, index);
	}
}

#ifndef MINIMAL_STEPWIDTH
#define MINIMAL_STEPWIDTH 0.000001
#endif
#ifndef MAXIMAL_STEPWIDTH
#define MAXIMAL_STEPWIDTH 0.75
#endif

void rmw_adapt_stepwidth(mcmc * m, double prob_old) {
	unsigned int i;
	double step;
	double min;
	double scale;
	double max;
	double alpha = exp(get_prob(m) - prob_old);
	if (alpha > 1)
		alpha = 1;
	for (i = 0; i < get_n_par(m); i++) {
		scale = (gsl_vector_get(m->params_max, i)
				- gsl_vector_get(m->params_min, i));
		min = MINIMAL_STEPWIDTH * scale;
		max = MAXIMAL_STEPWIDTH * scale;

		step = gsl_vector_get(get_steps(m), i);
		step += get_next_uniform_random(m) / sqrt(m->n_iter) * (alpha - 0.234)
				* scale;
		;
		if (step < min)
			step = min;
		if (step > max)
			step = max;
		gsl_vector_set(get_steps(m), i, step);
	}
}

void markov_chain_step(mcmc * m, int calc_index) {
	double prob_old = get_prob(m);
	gsl_vector * old_model = dup_vector(m->model);
	gsl_vector * old_values = dup_vector(m->params);

	mcmc_check(m);
	do_step(m);

	if (calc_index <= 1) {
		calc_model(m, old_values);
	}

	if (check_accept(m, prob_old) == 1) {
		inc_params_accepts(m);
		gsl_vector_free(old_model);
		gsl_vector_free(old_values);
	} else {
		revert(m, old_model, prob_old);
		set_params(m, old_values);
		inc_params_rejects(m);
	}
}
