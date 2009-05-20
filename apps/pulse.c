#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "debug.h"

static gsl_vector ** power_profile = NULL;

double calc_prior(const mcmc * m, double hmin) {
	unsigned int i;
	const gsl_vector * params = m->params;
	double prior = 0;
	double mode_height;
	for (i = 2; i < get_n_par(m); i += 2) {
		mode_height = gsl_vector_get(params, i + 1);
		prior += gsl_sf_log(mode_height + hmin);
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
	double hmin = pow(10, gsl_vector_get(params, 0));
	double prior = calc_prior(m, hmin);
	double lifetime = gsl_vector_get(params, 1);

	(void) old_values;
	if (power_profile != NULL)
		power_profile = NULL;
	debug("normal calculation");

	for (i = 0; i < m->x_dat->size; i++) {
		y = 0;
		freq = gsl_vector_get(m->x_dat, i);
		for (j = 2; j < get_n_par(m); j += 2) {
			mode_height = gsl_vector_get(params, j);
			mode_freq = gsl_vector_get(params, j + 1);
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
	double hmin;
	double distance;
	gsl_vector * params = m->params;
	double lifetime = gsl_vector_get(params, 1);

	if (power_profile == NULL) {
		power_profile = mem_calloc(get_n_par(m)-2, sizeof(gsl_vector*));
		assert(power_profile != NULL);
	}

	dump_i("parameter changed", j);
	if (j == 0) {
		/* hmin changed -> recalculate prior */
		hmin = pow(10, gsl_vector_get(params, 0));
		for (i = 2; i < get_n_par(m); i += 2) {
			mode_height = gsl_vector_get(params, i + 1);
			diffprob += gsl_sf_log(mode_height + hmin);
			diffprob -= gsl_sf_log(mode_height + pow(10, old_value));
		}
		set_prob(m, get_beta(m) * -1 * (prob + diffprob));
		return;
	} else if (j == 1) {
		/* lifetime changed -> recalculate everything */
		for (i = 2; i < get_n_par(m); i += 2) {
			if (power_profile[(i - 2) / 2] == NULL) {
				power_profile[(i - 2) / 2] = gsl_vector_alloc(m->x_dat->size);
			}
		}
		hmin = pow(10, gsl_vector_get(params, 0));
		prob = calc_prior(m, hmin);
		for (i = 0; i < m->x_dat->size; i++) {
			y = 0;
			freq = gsl_vector_get(m->x_dat, i);
			for (k = 2; k < get_n_par(m); k += 2) {
				mode_height = gsl_vector_get(params, k);
				mode_freq = gsl_vector_get(params, k + 1);
				distance = mode_freq - freq;
				gsl_vector_set(power_profile[(k - 2) / 2], i, mode_height / (1
						+ 4 * distance * distance * M_PI * lifetime));
				y += gsl_vector_get(power_profile[(k - 2) / 2], i);
			}
			gsl_vector_set(m->model, i, y);

			prob += gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
		}
		set_prob(m, get_beta(m) * -1 * prob);
		return;
	} else if ((j - 2) % 2 == 0) {
		/* a height parameter changed -> recalculate with new factor */
		if (power_profile[(j - 2) / 2] != NULL) {
			mode_height = gsl_vector_get(params, j);
			for (i = 0; i < m->x_dat->size; i++) {
				y = gsl_vector_get(m->model, i);
				prob -= gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
				y -= gsl_vector_get(power_profile[(j - 2) / 2], i);
				y += mode_height / old_value * gsl_vector_get(power_profile[(j
						- 2) / 2], i);
				gsl_vector_set(m->model, i, y);
				prob += gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
			}
			gsl_vector_scale(power_profile[(j - 2) / 2], mode_height
					/ old_value);
			set_prob(m, get_beta(m) * -1 * (prob + diffprob));
			return;
		}
	}
	if (power_profile[(j - 2) / 2] == NULL) {
		power_profile[(j - 2) / 2] = gsl_vector_alloc(m->x_dat->size);
	}

	/* a frequency parameter changed. -> recalculate this */
	for (i = 0; i < m->x_dat->size; i++) {
		y = gsl_vector_get(m->model, i);
		diffprob -= gsl_sf_log(y) + gsl_vector_get(m->y_dat, i) / y;
		freq = gsl_vector_get(m->x_dat, i);
		mode_height = gsl_vector_get(params, (j / 2) * 2);
		mode_freq = gsl_vector_get(params, (j / 2) * 2 + 1);
		distance = mode_freq - freq;
		gsl_vector_set(power_profile[(j - 2) / 2], i, 1 / (1 + 4 * distance
				* distance * M_PI * lifetime));
		y += mode_height * gsl_vector_get(power_profile[(j - 2) / 2], i);
		if ((j - 2) % 2 == 0) {
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

