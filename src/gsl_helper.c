
#include "gsl_helper.h"
#include "debug.h"

gsl_histogram * calc_hist(const gsl_vector * v, int nbins) {
	double max;
	double min;
	unsigned int i;
	double binwidth;
	double sum = 0;
	double val;
	gsl_histogram * h;

	gsl_vector_minmax(v, &min, &max);
	binwidth = (max - min) / nbins;
	dump_d("min", min);
	dump_d("max", max);

	debug("allocating the histogram");
	h = gsl_histogram_alloc(v->size);
	debug("setting range");
	require(gsl_histogram_set_ranges_uniform (h, min, max));

	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += 1;

	debug("summing up");
	for (i = 0; i < v->size; i++) {
		val = gsl_vector_get(v, i);
		sum += val;
		require(gsl_histogram_increment (h, val));
	}
	debug("scaling");
	/* double gsl_histogram_sum (const gsl_histogram * h) */
	require(gsl_histogram_scale (h, 1/sum));
	debug("done");
	return h;
}

double calc_vector_sum(const gsl_vector * v) {
	double sum = 0;
	unsigned int i;
	for (i = 0; i < v->size; i++) {
		sum += gsl_vector_get(v, i);
	}
	return sum;
}

double calc_vector_squaresum(const gsl_vector * v) {
	static double sum = 0;
	static double x;
	static unsigned int i;
	sum = 0;
	for (i = 0; i < v->size; i++) {
		x = gsl_vector_get(v, i);
		sum += x * x;
	}
	return sum;
}

gsl_vector * dup_vector(const gsl_vector * v) {
	gsl_vector * r;
	assert(v != NULL);
	assert(v->size > 0);
	r = gsl_vector_alloc(v->size);
	assert(r != NULL);
	require(gsl_vector_memcpy(r, v));
	return r;
}

gsl_vector * calc_normalized(const gsl_vector * v) {
	double sum = calc_vector_sum(v);
	gsl_vector * r = dup_vector(v);
	require(gsl_vector_scale(r, 1/sum));
	return r;
}

int calc_same(const gsl_vector * a, const gsl_vector * b) {
	unsigned int i;
	assert(a->size == b->size);

	if (a == b)
		return 1;

	for (i = 0; i < a->size; i++) {
		if (gsl_vector_get(a, i) != gsl_vector_get(b, i))
			return 0;
	}
	return 1;
}

void max_vector(gsl_vector * a, const gsl_vector * b) {
	unsigned int i;
	assert(a->size == b->size);

	if (a == b)
		return;

	for (i = 0; i < a->size; i++) {
		if (gsl_vector_get(a, i) < gsl_vector_get(b, i))
			gsl_vector_set(a, i, gsl_vector_get(b, i));
	}
}

void min_vector(gsl_vector * a, const gsl_vector * b) {
	unsigned int i;
	assert(a->size == b->size);

	if (a == b)
		return;

	for (i = 0; i < a->size; i++) {
		if (gsl_vector_get(a, i) > gsl_vector_get(b, i))
			gsl_vector_set(a, i, gsl_vector_get(b, i));
	}
}
