#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include "gsl_helper.h"
#include "debug.h"

gsl_rng * rng;

void setup_rng() {
	gsl_rng_env_setup();
	rng = gsl_rng_alloc(gsl_rng_default);
}
/* this is a singleton */
gsl_rng * get_rng_instance() {
	return rng;
}

gsl_vector * get_random_uniform_vector(unsigned int size) {
	unsigned int j;
	gsl_vector * v = gsl_vector_alloc(size);
	for (j = 0; j < size; j++) {
		gsl_vector_set(v, j, gsl_rng_uniform(get_rng_instance()));
	}
	return v;
}

/*
 * a expensive function, that can be assumed to be
 * always non-constant in at least one dimension, except on the maximum to
 * be found
 */
double f(gsl_vector * x);

gsl_vector ** cache = NULL;
double * values = NULL;
int cachesize = 0;

double f_cached(gsl_vector * x, double(* intern_f)(gsl_vector*)) {
	int i;

	for (i = 0; i < cachesize; i++) {
		if (calc_same(x, cache[i]) == 1) {
			return values[i];
		}
	}
	cachesize++;
	cache = (gsl_vector**) realloc(cache, cachesize * sizeof(gsl_vector*));
	assert(cache != NULL);
	cache[cachesize - 1] = dup_vector(x);
	values = (double*) realloc(values, cachesize * sizeof(double));
	assert(values != NULL);
	values[cachesize - 1] = intern_f(x);
	return values[cachesize - 1];
}

/* limit the space of v to [0..1] in each dimension */
void limit(gsl_vector * v) {
	gsl_vector * low = gsl_vector_alloc(v->size);
	gsl_vector * high = gsl_vector_alloc(v->size);
	gsl_vector_set_all(low, 0.0);
	gsl_vector_set_all(high, 1.0);
	max_vector(v, low);
	min_vector(v, high);
	gsl_vector_free(low);
	gsl_vector_free(high);
}
#ifndef RANDOM_SCALE_CIRCLE_JUMP
#ifdef ADAPTIVE
#define RANDOM_SCALE_CIRCLE_JUMP 1.83
#else
#define RANDOM_SCALE_CIRCLE_JUMP 1.0
#endif
#endif
#ifndef JUMP_SCALE
#define JUMP_SCALE 1.0
#endif
#ifndef RANDOM_SCALE
#ifdef ADAPTIVE
#define RANDOM_SCALE 0.001
#else
#define RANDOM_SCALE 0.04
#endif
#endif

int detect_circle_jump(gsl_vector * flaps, gsl_vector * probe_values) {
	int i;
	int ndim = flaps->size;
	int possibly_circle_jump = 0;
	for (i = 0; i < ndim; i++) {
		if (gsl_vector_get(probe_values, i) > 0) {
			if (gsl_vector_get(flaps, i) == 1) {
				possibly_circle_jump = 1;
				continue; /* we jump back */
			} else {
				/* a jump forward */
				possibly_circle_jump = 0;
				break;
			}
		} else {
			if (gsl_vector_get(flaps, i) == 0) {
				/* turning around after a jump */
				possibly_circle_jump = 1;
				continue;
			} else {
				/* turning around after turning around */
				if (gsl_vector_get(flaps, i) == 2) {
					possibly_circle_jump = 1;
					continue; /* and we are inconclusive in the others */
				} else {
					possibly_circle_jump = 0;
					break;
				}
			}
		}
	}
	if (i != ndim || ndim == 1) {
		possibly_circle_jump = 0;
	} else if (possibly_circle_jump == 1) {
		debug("circle-jump possible. increased randomness");
	}
	return possibly_circle_jump;
}

/* find a local maximum (climb the hill)
 * diagonals */
int find_local_maximum_multi(unsigned int ndim, double exactness,
		gsl_vector * start) {
	unsigned int i;
	unsigned int count = 0;
	int possibly_circle_jump;
	double current_val;
	gsl_vector * current_probe = gsl_vector_alloc(ndim);
	gsl_vector * next_probe = gsl_vector_alloc(ndim);
	gsl_vector * current_x = dup_vector(start);
	gsl_vector * scales = gsl_vector_alloc(ndim);
	/* did we switch direction in the last move? */
	gsl_vector * flaps = gsl_vector_alloc(ndim);
	gsl_vector * probe_values = gsl_vector_alloc(ndim);
	gsl_vector_set_all(scales, 1.0 / 3);
	gsl_vector_set_all(flaps, 0);
	assert(exactness < 1);

	while (1) {
		dump_v("currently at", current_x)
		current_val = f(current_x);
		count++;
		dump_d("current value", current_val);

		gsl_vector_memcpy(next_probe, current_x);
		gsl_vector_add(next_probe, scales);
		limit(next_probe);
		dump_v("will probe at", next_probe);

		for (i = 0; i < ndim; i++) {
			gsl_vector_memcpy(current_probe, current_x);
			gsl_vector_set(current_probe, i, gsl_vector_get(next_probe, i));

			gsl_vector_set(probe_values, i, f(current_probe) - current_val);
			if (gsl_vector_get(probe_values, i) < 0)
				gsl_vector_set(probe_values, i, 0);
			count++;
		}
		if(gsl_vector_max(probe_values) != 0)
			gsl_vector_scale(probe_values, 1 / gsl_vector_max(probe_values));
		dump_v("probe results", probe_values);
		gsl_vector_memcpy(start, current_x);

		possibly_circle_jump = detect_circle_jump(flaps, probe_values);

		for (i = 0; i < ndim; i++) {
			if (gsl_vector_get(probe_values, i) > 0) {
				dump_i("we jump forward in", i);
				gsl_vector_set(
						current_x,
						i,
						gsl_vector_get(current_x, i)
								+ gsl_vector_get(scales, i) * JUMP_SCALE
										*
#ifdef ADAPTIVE
										gsl_vector_get(probe_values, i) *
#endif
										(1
												+ (gsl_rng_uniform(
														get_rng_instance())
														- 0.5) * 2
														* (RANDOM_SCALE
																+ possibly_circle_jump
																		* RANDOM_SCALE_CIRCLE_JUMP)));
				limit(current_x);
				if (gsl_vector_get(current_x, i) == gsl_vector_get(start, i)) {
					/* we clashed against a wall. That means we are ready to
					 * refine */
					gsl_vector_set(flaps, i, 2);
				} else {
					gsl_vector_set(flaps, i, 0);
				}
			} else {
				if (gsl_vector_get(flaps, i) == 0) {
					dump_i("we turn back in", i);
					gsl_vector_set(flaps, i, 1);
					/* TODO: should we step back a little?
					 * no we can't, otherwise our double-turnback is tainted */
					gsl_vector_set(scales, i, gsl_vector_get(scales, i) * -1);
				} else {
					dump_i("we turned back twice in", i);
					gsl_vector_set(flaps, i, 2);
				}
			}
		}
		if (gsl_vector_min(flaps) == 2) {
			debug("all dimensions are ready, lets refine");
			dump_d("exactness (min)", gsl_vector_min(scales));
			dump_d("exactness (max)", gsl_vector_max(scales));
			dump_d("exactness (desired)", exactness);
			if (gsl_vector_max(scales) < exactness && gsl_vector_min(scales)
					> -exactness) {
				for (i = 0; i < ndim; i++) {
					gsl_vector_memcpy(current_probe, start);
					gsl_vector_set(current_probe, i, gsl_vector_get(
							current_probe, i) + abs(gsl_vector_get(scales, i)));
					assert(f(current_probe) >= current_val);
					gsl_vector_set(current_probe, i, gsl_vector_get(
							current_probe, i) - 2* abs (gsl_vector_get(scales,
							i)));
					assert(f(current_probe) >= current_val);
				}
				gsl_vector_free(scales);
				gsl_vector_free(flaps);
				gsl_vector_free(probe_values);
				gsl_vector_free(current_probe);
				gsl_vector_free(next_probe);
				gsl_vector_free(current_x);
				return count;
			}
			gsl_vector_scale(scales, 0.5);
			gsl_vector_set_all(flaps, 0);
			dump_d("new exactness (min)", gsl_vector_min(scales));
			dump_d("new exactness (max)", gsl_vector_max(scales));
		}
	}
}

/* find a local maximum (climb the hill)
 * one probe at a time */
int find_local_maximum_naive(unsigned int ndim, double exactness,
		gsl_vector * current_x) {
	unsigned int i = 0;
	unsigned int last_i = 0;
	unsigned int j = 0;
	unsigned int count = 0;
	double current_val;
	gsl_vector * current_probe = gsl_vector_alloc(ndim);
	gsl_vector * next_probe = gsl_vector_alloc(ndim);
	gsl_vector * scales = gsl_vector_alloc(ndim);
	/* did we switch direction in the last move? */
	int flaps = 0;
	double probe_value;
	gsl_vector_set_all(scales, 1.0 / 3);

	current_val = f(current_x);
	count++;

	dump_v("currently at", current_x)
	dump_d("current value", current_val);
	while (1) {
		for (j = 0; j < ndim; j++) {
			i = (last_i + j) % ndim;
			gsl_vector_memcpy(current_probe, current_x);
			gsl_vector_set(current_probe, i, gsl_vector_get(current_probe, i)
					+ gsl_vector_get(scales, i));
			limit(current_probe);
			if (calc_same(current_probe, current_x) == 1) {
				dump_i("we clashed a wall in", i);
				gsl_vector_set(scales, i, gsl_vector_get(scales, i) * -1);
				continue;
			}
			dump_v("will probe at", current_probe);
			probe_value = f(current_probe);
			count++;
			if (probe_value > current_val) {
				dump_i("we jump forward in", i);
				current_val = probe_value;
				gsl_vector_memcpy(current_x, current_probe);
				dump_v("currently at", current_x)
				dump_d("current value", current_val);
				break;
			} else {
				dump_i("we turn back in", i);
				gsl_vector_set(scales, i, gsl_vector_get(scales, i) * -1);
			}
		}
#ifdef ROUNDROBIN
		last_i = i;
#else
		last_i = 0;
#endif
		if (j == ndim) {
			if (flaps == 1) {
				debug("all dimensions are ready, lets refine");
				dump_v("currently at", current_x)
				dump_d("exactness (min)", gsl_vector_min(scales));
				dump_d("exactness (max)", gsl_vector_max(scales));
				dump_d("exactness (desired)", exactness);
				if (gsl_vector_max(scales) < exactness
						&& gsl_vector_min(scales) > -exactness) {
					for (i = 0; i < ndim; i++) {
						gsl_vector_memcpy(current_probe, current_x);
						gsl_vector_set(current_probe, i, gsl_vector_get(
								current_probe, i) + abs(gsl_vector_get(scales,
								i)));
						assert(f(current_probe) >= current_val);
						gsl_vector_set(current_probe, i, gsl_vector_get(
								current_probe, i) - 2* abs (gsl_vector_get(
								scales, i)));
						assert(f(current_probe) >= current_val);
					}
					gsl_vector_free(scales);
					gsl_vector_free(current_probe);
					gsl_vector_free(next_probe);
					return count;
				}
				gsl_vector_scale(scales, 0.5);
				flaps = 0;
				dump_d("new exactness (min)", gsl_vector_min(scales));
				dump_d("new exactness (max)", gsl_vector_max(scales));
			} else {
				flaps = 1;
			}
		} else {
			flaps = 0;
		}
	}
	return count;
}

int find_local_maximum(unsigned int ndim, double exactness, gsl_vector * start) {
	assert(ndim > 1);
	assert(exactness > 0 && exactness < 1);
#ifdef NAIVE
	return find_local_maximum_naive(ndim, exactness, start);
#else
	return find_local_maximum_multi(ndim, exactness, start);
#endif
}
