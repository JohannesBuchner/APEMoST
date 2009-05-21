#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

#ifndef HMIN
#define HMIN 1e-6
#endif

static gsl_vector ** power_profile = NULL;

static double calc_prior(const mcmc * m) {
	unsigned int i;
	const gsl_vector * params = m->params;
	double prior = 0;
	double mode_height;
	for (i = 1; i < get_n_par(m); i += 2) {
		mode_height = gsl_vector_get(params, i + 1);
		prior += gsl_sf_log(mode_height + HMIN);
	}
	return prior;
}

void calc_model(mcmc * m, const gsl_vector * old_values) {
	unsigned int i;
	unsigned int j;
	double prob = 0;
	double y = 0;
	double mode_height;
	double mode_freq;
	double distance;
	double freq;
	const gsl_vector * params = m->params;
	double prior = calc_prior(m);
	double lifetime = gsl_vector_get(params, 0);

	(void) old_values;
	if (power_profile != NULL)
		power_profile = NULL;

	for (i = 0; i < m->x_dat->size; i++) {
		y = 0;
		freq = gsl_vector_get(m->x_dat, i);
		for (j = 1; j < get_n_par(m); j += 2) {
			mode_freq = gsl_vector_get(params, j);
			mode_height = gsl_vector_get(params, j + 1);
			distance = mode_freq - freq;
			y += mode_height / (1 + 4 * distance * distance * M_PI * lifetime);
		}
		gsl_vector_set(m->model, i, y);

		prob += gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
	}

	set_prob(m, get_beta(m) * -1 * (prob + prior));
}

void calc_model_for(mcmc * m, const unsigned int j, const double old_value) {
#ifdef NO_PARTIAL
	(void)j;
	(void)old_value;
	calc_model(m, NULL);
	return;
#else
	unsigned int i;
	unsigned int k;
	double diffprob = 0;
	double y;
	double prob = -1 * get_prob(m) / get_beta(m);
	double mode_height;
	double freq;
	double mode_freq;
	double distance;
	gsl_vector * params = m->params;
	double lifetime = gsl_vector_get(params, 0);
	int mode = (j - 1) / 2;

	if (power_profile == NULL) {
		power_profile = mem_calloc(get_n_par(m) - 1, sizeof(gsl_vector*));
		assert(power_profile != NULL);
	}

	if (j == 0) {
		/* lifetime changed -> recalculate everything */
		for (i = 1; i < get_n_par(m); i += 2) {
			if (power_profile[(i - 1) / 2] == NULL) {
				power_profile[(i - 1) / 2] = gsl_vector_alloc(m->x_dat->size);
			}
		}
		prob = calc_prior(m);
		for (i = 0; i < m->x_dat->size; i++) {
			y = 0;
			freq = gsl_vector_get(m->x_dat, i);
			for (k = 1; k < get_n_par(m); k += 2) {
				mode_freq = gsl_vector_get(params, k);
				mode_height = gsl_vector_get(params, k + 1);
				distance = mode_freq - freq;
				gsl_vector_set(power_profile[(k - 1) / 2], i, mode_height / (1
						+ 4 * distance * distance * M_PI * lifetime));
				y += gsl_vector_get(power_profile[(k - 1) / 2], i);
			}
			gsl_vector_set(m->model, i, y);

			prob += gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
		}
		set_prob(m, get_beta(m) * -1 * prob);
		return;
	} else if ((j - 1) % 2 == 1) {
		/* a height parameter changed -> recalculate with new factor */
		if (0 && power_profile[mode] != NULL) {
			mode_height = gsl_vector_get(params, j);
			for (i = 0; i < m->x_dat->size; i++) {
				y = gsl_vector_get(m->model, i);
				prob -= gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
				y -= gsl_vector_get(power_profile[mode], i);
				y += mode_height / old_value * gsl_vector_get(
						power_profile[mode], i);
				gsl_vector_set(m->model, i, y);
				prob += gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
			}
			gsl_vector_scale(power_profile[mode], mode_height / old_value);
			set_prob(m, get_beta(m) * -1 * (prob + diffprob));
			return;
		}
	}
	if (power_profile[(j - 1) / 2] == NULL) {
		power_profile[(j - 1) / 2] = gsl_vector_alloc(m->x_dat->size);
	}

	/* a frequency parameter changed. -> recalculate this */
	for (i = 0; i < m->x_dat->size; i++) {
		y = gsl_vector_get(m->model, i);
		diffprob -= gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
		freq = gsl_vector_get(m->x_dat, i);
		mode_freq = gsl_vector_get(params, 1 + mode * 2);
		mode_height = gsl_vector_get(params, 1 + mode * 2 + 1);
		distance = mode_freq - freq;
		gsl_vector_set(power_profile[mode], i, 1 / (1 + 4 * distance
				* distance * M_PI * lifetime));
		y += mode_height * gsl_vector_get(power_profile[mode], i);
		if ((j - 1) % 2 == 0) {
			y -= old_value / (1 + 4 * distance * distance * M_PI * lifetime);
		} else {
			distance = old_value - freq;
			y -= mode_height / (1 + 4 * distance * distance * M_PI * lifetime);
		}
		gsl_vector_set(m->model, i, y);
		diffprob += gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
	}
	set_prob(m, get_beta(m) * -1 * (prob + diffprob));
#endif
}

