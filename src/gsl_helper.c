#include "gsl_helper.h"
#include "utils.h"
#include "debug.h"
#include <gsl/gsl_linalg.h>

gsl_histogram * create_hist(int nbins, double min, double max) {
	gsl_histogram * h;

	h = gsl_histogram_alloc(nbins);
	gsl_histogram_set_ranges_uniform(h, min, max);

	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += 1;
	return h;
}

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

void append_to_hists(gsl_histogram ** hists, unsigned int n,
		const char * filename) {
	FILE * input;
	unsigned int i;
	int col;
	double v;
	int line = 0;

	input = openfile(filename);

	while (!feof(input)) {
		for (i = 0; i < n; i++) {
			col = fscanf(input, "%lf", &v);
			if (col != 1) {
				if (feof(input))
					break;
				fprintf(stderr, "field could not be read: %d, line %d in %s\n",
						i + 1, line + 1, filename);
				exit(1);
			}
			gsl_histogram_increment(hists[i], v);
		}
		line++;
	}
	dump_i("read lines", line);
	fclose(input);
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

void sort(gsl_vector ** vs, unsigned int nvectors, unsigned int vector_size) {
	unsigned int i = 0;
	unsigned int j = 0;
	unsigned int best = 0;
	double temp;

	for (j = 0; j < vector_size; j++) {
		best = j;
		for (i = j + 1; i < vector_size; i++) {
			if (gsl_vector_get(vs[0], i) > gsl_vector_get(vs[0], best)) {
				best = i;
			}
		}
		if (j != best) {
			/* switch j and best */
			for (i = 0; i < nvectors; i++) {
				temp = gsl_vector_get(vs[i], j);
				gsl_vector_set(vs[i], j, gsl_vector_get(vs[i], best));
				gsl_vector_set(vs[i], best, temp);
			}
		}
		if (j % 1000 == 0)
			IFDEBUG
				printf("sorted: %u/%u\r", j, vector_size);
	}
}

double min_column(const gsl_matrix * m, const unsigned int i) {
	unsigned int j = 1;
	double min = gsl_matrix_get(m, i, 0);
	while (j < m->size2) {
		if (min > gsl_matrix_get(m, i, j))
			min = gsl_matrix_get(m, i, j);
		j++;
	}
	return min;
}
double min_row(const gsl_matrix * m, const unsigned int i) {
	unsigned int j = 1;
	double min = gsl_matrix_get(m, 0, i);
	while (j < m->size2) {
		if (min > gsl_matrix_get(m, j, i))
			min = gsl_matrix_get(m, j, i);
		j++;
	}
	return min;
}
double max_column(const gsl_matrix * m, const unsigned int i) {
	unsigned int j = 1;
	double max = gsl_matrix_get(m, i, 0);
	while (j < m->size2) {
		if (max < gsl_matrix_get(m, i, j))
			max = gsl_matrix_get(m, i, j);
		j++;
	}
	return max;
}
double max_row(const gsl_matrix * m, const unsigned int i) {
	unsigned int j = 1;
	double max = gsl_matrix_get(m, 0, i);
	while (j < m->size2) {
		if (max < gsl_matrix_get(m, j, i))
			max = gsl_matrix_get(m, j, i);
		j++;
	}
	return max;
}

double xbar(const gsl_vector * x) {
	unsigned int i;
	double sum = 0;

	for (i = 0; i < x->size; i++) {
		sum += gsl_vector_get(x, i);
	}
	return sum / x->size;
}
double xbar_j(const gsl_matrix * x, const unsigned int j) {
	unsigned int i;
	double sum = 0;

	for (i = 0; i < x->size2; i++) {
		sum += gsl_matrix_get(x, j, i);
	}
	return sum / x->size1;
}

gsl_vector * linreg_n(const gsl_matrix * x, const gsl_vector * y, double * d,
		const gsl_vector * weights) {
	unsigned int n = y->size;
	unsigned int m = x->size2;
	gsl_matrix * xx = gsl_matrix_alloc(m, m);
	gsl_vector * xy = gsl_vector_alloc(m);
	gsl_vector * k = gsl_vector_alloc(m);
	gsl_permutation * p = gsl_permutation_alloc(m);

	gsl_vector * xbars = gsl_vector_alloc(m);
	double ybar = 0;

	unsigned int i;
	unsigned int j;
	unsigned int l;
	double sum;
	double weightsum = 0;

	int s;

	IFDEBUG
		printf("linear regression of %d datapoints in %d dimensions\n", n, m);

	assert(x->size1 == n);
	assert(n >= m);
	weightsum = 0;
	for (i = 0; i < n; i++) {
		ybar += gsl_vector_get(y, i) * gsl_vector_get(weights, i);
		weightsum += gsl_vector_get(weights, i);
	}
	ybar = ybar / weightsum;
	for (l = 0; l < m; l++) {
		gsl_vector_set(xbars, l, xbar_j(x, l));
	}

	for (l = 0; l < m; l++) {
		for (j = 0; j < m; j++) {
			if (j < l) {
				/* symmetric matrix */
				gsl_matrix_set(xx, l, j, gsl_matrix_get(xx, j, l));
			} else {
				sum = 0;
				for (i = 0; i < n; i++) {
					sum += (gsl_matrix_get(x, i, l) - gsl_vector_get(xbars, l))
							* (gsl_matrix_get(x, i, j) - gsl_vector_get(xbars,
									j));
				}
				gsl_matrix_set(xx, l, j, sum / n);
			}
		}
		sum = 0;
		weightsum = 0;
		for (i = 0; i < n; i++) {
			sum += (gsl_matrix_get(x, i, l) - gsl_vector_get(xbars, l))
					* (gsl_vector_get(y, i) - ybar)
					* gsl_vector_get(weights, i);
			weightsum += gsl_vector_get(weights, i);
		}
		gsl_vector_set(xy, l, sum * n / weightsum);
	}

	gsl_linalg_LU_decomp(xx, p, &s);

	gsl_linalg_LU_solve(xx, p, xy, k);

	/*
	 for (i = 0; i < m; i++) {
	 for (j = 0; j < m; j++) {
	 printf("%g\t", gsl_matrix_get(xx, i, j));
	 }
	 printf("| %g \t|| %g \t|\n", gsl_vector_get(k, i),
	 gsl_vector_get(xy, i));
	 }*/

	*d = xbar(y);
	for (j = 0; j < m; j++) {
		*d -= gsl_vector_get(k, j) * gsl_vector_get(xbars, j);
	}

	gsl_permutation_free(p);
	gsl_vector_free(xy);
	gsl_matrix_free(xx);
	gsl_vector_free(xbars);

	return k;
}

double calc_deviation(const gsl_matrix * x, const gsl_vector * y,
		const gsl_vector * k, const double d, const gsl_vector * weights) {
	double sumsq = 0;
	unsigned int i;
	unsigned int j;
	unsigned int n = x->size1;
	unsigned int m = x->size2;
	double zi;
	double weightsum = 0;

	assert(n == y->size);

	for (i = 0; i < n; i++) {
		zi = d;
		for (j = 0; j < m; j++) {
			zi += gsl_matrix_get(x, i, j) * gsl_vector_get(k, j);
		}

		sumsq += gsl_vector_get(weights, i) * pow(gsl_vector_get(y, i) - zi, 2);
		weightsum += gsl_vector_get(weights, i);
	}
	return sqrt(sumsq / weightsum);
}
