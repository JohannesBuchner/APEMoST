

#include "gsl_helper.h"
#include "debug.h"

/*
 * PTR_NEW:
 *   creates a pointer to a heap variable
 *   argument is copied and a pointer is returned.
 **
 * with foo and bar being pointers
 *   bar = foo lets bar point to the same values as foo
 *   *bar = foo changes the values of bar to the same ones as foo (copy)
 *   * does dereferencing
 **
 * DBLARR:
 *   creates a double/float array/vectro
 *   argument: dimension, up to 8
 **
 * (*var)[a,b]
 *   dereferences pointer, accesses index a of var, and of that accesses index b
 * (*var)[a,*]
 *   full vector of 2D array (*var), var being a pointer to the 2D-array
 **
 * ?INDGEN:
 *   generates a array where each entry has the value of its subscript.
 *   examples: findgen (float), dindgen (double)
 *
 **
 * Histograms:
 *   count stuff that falls into bins.
 *   return the number of elements in the bin for each bin.
 *   usual arguments: nbins, min, max
 */

gsl_histogram * calc_hist(gsl_vector ** param_distr, int index, int nbins) {
	double max;
	double min;
	unsigned int i;
	double binwidth;
	double sum = 0;
	double v;
	gsl_histogram * h;
	gsl_vector * dat = param_distr[index];

	gsl_vector_minmax (dat, &min, &max);
	binwidth = (max - min)/nbins;
	dump_d("min", min);
	dump_d("max", max);

	debug("allocating the histogram");
	h = gsl_histogram_alloc (dat->size);
	debug("setting range");
	require(gsl_histogram_set_ranges_uniform (h, min, max));

	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += 1;

	debug("summing up");
	for(i=0; i<dat->size; i++){
		v = gsl_vector_get (dat, i);
		sum += v;
		require(gsl_histogram_increment (h, v));
	}
	debug("scaling");
	/* double gsl_histogram_sum (const gsl_histogram * h) */
	require(gsl_histogram_scale (h, 1/sum));
	debug("done");
	return h;
}

double calc_vector_sum(gsl_vector * v) {
	double sum = 0;
	unsigned int i;
	for(i = 0; i < v->size; i++) {
		sum += gsl_vector_get(v, i);
	}
	return sum;
}

gsl_vector * dup_vector(gsl_vector * v) {
	gsl_vector * r = gsl_vector_alloc(v->size);
	require(gsl_vector_memcpy(r, v));
	return r;
}

gsl_vector * calc_normalized(gsl_vector * v) {
	double sum = calc_vector_sum(v);
	gsl_vector * r = dup_vector(v);
	require(gsl_vector_scale(r, 1/sum));
	return r;
}


