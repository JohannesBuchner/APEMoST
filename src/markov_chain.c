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

#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"
#include <gsl/gsl_sf.h>

void restart_from_best(mcmc * m) {
	set_params(m, dup_vector(get_params_best(m)));
	set_prob(m, get_prob_best(m));
}

void burn_in(mcmc * m, const unsigned int burn_in_iterations) {
	unsigned long iter;
	unsigned long subiter;

	gsl_vector * original_steps = get_steps(m);
	m->params_step = dup_vector(m->params_max);
	gsl_vector_sub(m->params_step, m->params_min);
	gsl_vector_scale(m->params_step, 0.1);

	debug("Beginning calibration of MCMC ...");
	debug("Starting burn-in ...");
	mcmc_check(m);
	for (iter = 0; iter < burn_in_iterations / 2; ) {
		for (subiter = 0; subiter < 200; subiter++) {
			markov_chain_step(m);
		}
		iter += subiter;
		dump_ul("\tBurn-in Iteration", iter);
		IFVERBOSE {
			dump_v("stepwidth", get_steps(m));
			dump_v("params", get_params(m));
		}
		mcmc_check_best(m);
	}
	debug("Re-initializing burn-in ...");
	restart_from_best(m);
	gsl_vector_scale(m->params_step, 0.5);
	for (; iter < burn_in_iterations; ) {
		for (subiter = 0; subiter < 200; subiter++) {
			markov_chain_step(m);
		}
		iter += subiter;
		dump_ul("\tBurn-in Iteration", iter);
		IFVERBOSE {
			dump_v("stepwidth", get_steps(m));
			dump_v("params", get_params(m));
		}
		mcmc_check_best(m);
	}
	debug("Burn-in done, adjusting steps ...");
	gsl_vector_memcpy(get_steps(m), original_steps);
	gsl_vector_free(original_steps);
	mcmc_check(m);
	debug("Burn-in done.");

}

void clear_bit(char * bitfield, unsigned int i) {
	bitfield[i / 8] &= ~(1 << (i % 8));
}
void set_bit(char * bitfield, unsigned int i) {
	bitfield[i / 8] |= (1 << (i % 8));
}
int get_bit(char * bitfield, unsigned int i) {
	return bitfield[i / 8] & (1 << (i % 8));
}

#ifndef ACCURACY_DEVIATION_FACTOR
/**
 * How good should the acceptance rate be calculated in dependence of
 * deviation from the desired value?
 * accuracy = factor * deviation
 */
#define ACCURACY_DEVIATION_FACTOR 0.25
#endif
/**
 * Get acceptance rate.
 * The closer the acceptance rate is to the desired acceptance rate, the more
 * accurately will it be assessed.
 *
 * @param m
 * @param param
 * @param desired_acceptance_rate
 * @param min_accuracy you can request a upper limit on the accuracy, e.g. 1%.
 * 		any calculation will have at most the accuracy of 1% then. Otherwise put 0 here.
 * @param max_accuracy you can request a lower limit on the accuracy, e.g. 0.1%
 * 		any calculation will have at least the accuracy of 0.1% then. Otherwise put 1 here.
 *
 * @param acceptance_rate here the a/r gets stored
 * @param accuracy here the accuracy gets stored
 *
 * @return iterations used
 */
unsigned int assess_acceptance_rate(mcmc * m, unsigned int param,
		double desired_acceptance_rate, double min_accuracy,
		double max_accuracy, double * acceptance_rate, double * accuracy) {
	unsigned int i = 0;
	unsigned int j;
	unsigned int n = 40;
	unsigned int n_par = get_n_par(m);
	unsigned long accepts = 0;
	double stdev = 0;
	unsigned int maxdev = 0;
	double accept_rate;
	double required_accuracy = min_accuracy;
	char * acceptslog = NULL;
	/*
	 gsl_vector * start_params = dup_vector(get_params(m));
	 double start_prob = get_prob(m);
	 */
	reset_accept_rejects(m);

	while (1) {
		IFVERBOSE
			printf("calculating %d steps.\n", n);
		acceptslog = (char*) realloc(acceptslog, n * sizeof(char));
		assert(acceptslog != NULL);

		for (; i < n; i++) {
			if (param < n_par) {
				accepts = get_params_accepts_for(m, param);
				/*restart_from_best(m);*/
				markov_chain_step_for(m, param);
				mcmc_check_best(m);
				/*set_prob(m, start_prob);
				 set_params(m, dup_vector(start_params));*/
				if (accepts == get_params_accepts_for(m, param)) {
					/* had a reject -> set bit to 0 */
					clear_bit(acceptslog, i);
				} else {
					/* had a accept -> set bit to 1 */
					set_bit(acceptslog, i);
				}
			} else {
				accepts = get_params_accepts_global(m);
				markov_chain_step(m);
				mcmc_check_best(m);
				if (accepts == get_params_accepts_global(m)) {
					clear_bit(acceptslog, i);
				} else {
					set_bit(acceptslog, i);
				}

			}
		}
		accept_rate = accepts / (double) n;
		IFVERBOSE
			printf("accept rate: %f (%lu/%d)\n", accept_rate, accepts, n);

		/* get max deviation */
		accepts = 0;
		stdev = 0;
		maxdev = 0 + 1;
		for (j = 0; j < n; j++) {
			if (get_bit(acceptslog, j) != 0) {
				accepts++;
			}
			stdev += pow(accepts - accept_rate * j, 2);
			if (abs(accepts - accept_rate * j) > maxdev) {
				maxdev = abs(accepts - accept_rate * j);
			}
		}
		stdev = sqrt(stdev / n) * 2;

		/*
		 * if we are way off, we don't need to be that accurate.
		 * if we are close, we want to be more accurate
		 * 30% could also be 10% to be more cautious
		 */
		required_accuracy = abs_double(accept_rate - desired_acceptance_rate)
				* ACCURACY_DEVIATION_FACTOR;
		if (required_accuracy < 0.005)
			required_accuracy = 0.005;
		if (required_accuracy < min_accuracy) {
			required_accuracy = min_accuracy;
		}
		if (required_accuracy > max_accuracy) {
			required_accuracy = max_accuracy;
		}

		/*
		 * we assume we have a deviation of maxdev at the end.
		 * how many values do we need to get below required_accuracy
		 */
		*acceptance_rate = accept_rate;
		*accuracy = maxdev / 1. / n;
		IFVERBOSE
			printf("accuracy wanted: %f, got: %f\n", required_accuracy,
					*accuracy);

		if (*accuracy <= required_accuracy) {
			break;
		}
		/*
		 * we need (maxdev / min_accuracy) values to achieve min_accuracy
		 */
		assert(maxdev / required_accuracy >= n);
		n = ((unsigned int) ((maxdev / 1. / required_accuracy) / 8) + 1) * 8;
	}
	return n;
}

void do_step_for(mcmc * m, const unsigned int i) {
	const double step = gsl_vector_get(m->params_step, i);
	const double old_value = gsl_vector_get(m->params, i);
	double new_value;
	const double max = gsl_vector_get(m->params_max, i);
	const double min = gsl_vector_get(m->params_min, i);
	/* dump_d("Jumping from", old_value); */

#if CIRCULAR_PARAMS == 0
	do {
		new_value = old_value + get_next_random_jump(m, step);
		IFVERBOSE
			printf("Value borders reached; looking for new starting point"
				" for %d \n", i);
	} while (new_value > max || new_value < min);
#else
	unsigned int j = 0;
	/** circular parameters **/
	unsigned int parameters[] = {CIRCULAR_PARAMS, 0};
	
	new_value = old_value + get_next_random_jump(m, step);
	
	if (new_value > max || new_value < min) {
		while (1) {
			if (parameters[j] == 0) {
				/* non-circular parameter */
				do {
					new_value = old_value + get_next_random_jump(m, step);
				}while (new_value > max || new_value < min);
				break;
			}
			if (parameters[j] == i + 1) {
				/* circular parameter */
				new_value = min + mod_double(new_value - min, max - min);
				break;
			}
			j++;
		}
	}
#endif
	assert(new_value <= max);
	assert(new_value >= min);
	/* dump_d("To", new_value); */
	set_params_for(m, new_value, i);
}

static void do_step(mcmc * m) {
	unsigned int i;
	for (i = 0; i < get_n_par(m); i++) {
		do_step_for(m, i);
	}
}

/**
 * @returns 1 if accept, 0 if rejecting
 */
static int check_accept(mcmc * m, const double prob_old) {
	double prob_new = get_prob(m);
	double prob_still_accept;

	/* shortcut */
	if (prob_new == prob_old) {
		return 1;
	}
	IFVERBOSE
		dump_v("suggesting parameter", get_params(m));

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

static void revert(mcmc * m, const double prob_old) {
	set_prob(m, prob_old);
}

void markov_chain_step_for(mcmc * m, const unsigned int index) {
	double prob_old = get_prob(m);
	double old_value = gsl_vector_get(m->params, index);

	mcmc_check(m);
	do_step_for(m, index);

	calc_model_for(m, index, old_value);

	if (check_accept(m, prob_old) == 1) {
		inc_params_accepts_for(m, index);
	} else {
		revert(m, prob_old);
		set_params_for(m, old_value, index);
		inc_params_rejects_for(m, index);
	}
}

#ifndef MINIMAL_STEPWIDTH
#define MINIMAL_STEPWIDTH 0.0000001
#endif
#ifndef MAXIMAL_STEPWIDTH
#define MAXIMAL_STEPWIDTH 1000000
#endif

void rmw_adapt_stepwidth(mcmc * m, const double prob_old) {
	unsigned int i;
	double step;
	double min;
	double scale;
	double max;
	double alpha = exp(get_prob(m) - prob_old);
	if (alpha > 1)
		alpha = 1;
	for (i = 0; i < get_n_par(m); i++) {
		scale = (gsl_vector_get(m->params_max, i) - gsl_vector_get(
				m->params_min, i));
		min = MINIMAL_STEPWIDTH * scale;
		max = MAXIMAL_STEPWIDTH * scale;

		step = gsl_vector_get(get_steps(m), i);
		step += get_next_uniform_random(m) / sqrt(m->n_iter) * (alpha
				- TARGET_ACCEPTANCE_RATE) * scale;
		;
		if (step < min)
			step = min;
		if (step > max)
			step = max;
		gsl_vector_set(get_steps(m), i, step);
	}
}

void markov_chain_step(mcmc * m) {
	double prob_old = get_prob(m);
	gsl_vector * old_values = dup_vector(m->params);

	mcmc_check(m);
	do_step(m);

	calc_model(m, old_values);

	if (check_accept(m, prob_old) == 1) {
		inc_params_accepts(m);
		gsl_vector_free(old_values);
	} else {
		revert(m, prob_old);
		set_params(m, old_values);
		inc_params_rejects(m);
	}
}
