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

gsl_vector ** generate_coefficients(unsigned int dimensions) {
	unsigned int i;
	unsigned int order;
	gsl_vector ** c = (gsl_vector**) calloc(dimensions, sizeof(gsl_vector*));

	for (i = 0; i < dimensions; i++) {
		order = gsl_rng_uniform(get_rng_instance()) * 10 + 1;
		assert(order <= 10);
		c[i] = get_random_uniform_vector(order);
		gsl_vector_add_constant(c[i], -0.5);
		gsl_vector_scale(c[i], 20);
		printf("  coeff %d: ", i);
		dump_vectorln(c[i]);
	}
	return c;
}
gsl_vector ** coefficients;
double my_polynome(gsl_vector * x) {
	double v = 1;
	double w;
	unsigned int i, j;
	for (i = 0; i < x->size; i++) {
		w = 0;
		for (j = 0; j < coefficients[i]->size; j++) {
			w += pow(gsl_vector_get(x, i), j) * gsl_vector_get(coefficients[i],
					j);
		}
		v *= w;
	}
	return v;
}

double get_vector_max(gsl_vector * v) {
	unsigned int i;
	double max;
	max = gsl_vector_get(v, 0);
	for (i = 1; i < v->size; i++) {
		if (max < gsl_vector_get(v, i))
			max = gsl_vector_get(v, i);
	}
	return max;
}
double get_vector_min(gsl_vector * v) {
	unsigned int i;
	double min;
	min = gsl_vector_get(v, 0);
	for (i = 1; i < v->size; i++) {
		if (min > gsl_vector_get(v, i))
			min = gsl_vector_get(v, i);
	}
	return min;
}

/*
 * a expensive function, that can be assumed to be
 * always non-constant in at least one dimension, except on the maximum to
 * be found
 */
double f(gsl_vector * x) {
	return my_polynome(x);
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

/* find a local maximum (climb the hill) */
int find_local_maximum(int ndim, double exactness, gsl_vector * start) {
	int i;
	int count = 0;
	int ignore_candidate;
	double current_val;
	gsl_vector * current_probe = gsl_vector_alloc(ndim);
	gsl_vector * next_probe = gsl_vector_alloc(ndim);
	gsl_vector * current_x = start;
	gsl_vector * scales = gsl_vector_alloc(ndim);
	/* did we switch direction in the last move? */
	gsl_vector * flaps = gsl_vector_alloc(ndim);
	gsl_vector * probe_values = gsl_vector_alloc(ndim);
	gsl_vector_set_all(scales, 1.0 / 3);
	gsl_vector_set_all(flaps, 0);

	while (1) {
		dump_vectorln(current_x);
		/*dump_v("currently at", current_x)*/
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
			count++;
		}
		dump_v("probe results", probe_values);

		/* avoid circle-jumps */

		ignore_candidate = -1;
		for (i = 0; i < ndim; i++) {
			if (gsl_vector_get(probe_values, i) > 0) {
				if (gsl_vector_get(flaps, i) == 1) {
					ignore_candidate = i;
					continue; /* we jump back */
				} else
					break;
			} else {
				if (gsl_vector_get(flaps, i) == 0) {
					continue;
				} else {
					if (gsl_vector_get(flaps, i) == 2) {
						continue; /* and we are inconclusive in the others */
					} else {
						break;
					}
				}
			}
		}
		if (ignore_candidate > 0) {
			dump_i("circle-jump protection in action for", ignore_candidate);
			gsl_vector_set(probe_values, ignore_candidate, -1);
		}

		for (i = 0; i < ndim; i++) {
			if (gsl_vector_get(probe_values, i) > 0) {
				dump_i("we jump forward in", i);
				gsl_vector_set(current_x, i, gsl_vector_get(current_x, i)
						+ gsl_vector_get(scales, i) * 1.0);
				limit(current_x);
				gsl_vector_set(flaps, i, 0);
			} else {
				if (gsl_vector_get(flaps, i) == 0) {
					dump_i("we turn back in", i);
					gsl_vector_set(flaps, i, 1);
					gsl_vector_set(scales, i, gsl_vector_get(scales, i) * -1);
				} else {
					dump_i("we turned back twice in", i);
					gsl_vector_set(flaps, i, 2);
					if (calc_vector_sum(flaps) == 2 * ndim) {
						debug("all dimensions are ready, lets refine");
						if (get_vector_min(scales) < exactness) {
							gsl_vector_free(scales);
							gsl_vector_free(flaps);
							gsl_vector_free(probe_values);
							gsl_vector_free(current_probe);
							gsl_vector_free(next_probe);
							return count;
						}
						gsl_vector_scale(scales, 0.5);
						dump_d("new exactness", get_vector_min(scales));
					}
				}
			}
		}
	}
}

void free_coefficients(unsigned int ndim) {
	unsigned int i;
	for (i = 0; i < ndim; i++) {
		gsl_vector_free(coefficients[i]);
	}
	free(coefficients);
}
int main(int argc, char ** argv) {
	unsigned int i;
	long count = 0;
	unsigned int limit = 100;
	gsl_vector * start;
	if(argc == 2) {
		limit = atoi(argv[1]);
	}
	setup_rng();
	for (i = 0; i < limit; i++) {
		coefficients = generate_coefficients(5);
		start = get_random_uniform_vector(5);
		count += find_local_maximum(5, 0.001, start);
		free_coefficients(5);
		gsl_vector_free(start);
	}
	gsl_rng_free(get_rng_instance());
	printf("%lu steps\n", count);
	return 0;
}
