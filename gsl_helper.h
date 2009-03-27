#ifndef MCMC_GSL_HELPER
#define MCMC_GSL_HELPER

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>


#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif


/**
 * calculate a histogram
 * @param v vector to look at
 * @param nbins number of bins to use for the histogram
 */
gsl_histogram * calc_hist(gsl_vector * v, int nbins);

/**
 * sums the values
 */
double calc_vector_sum(gsl_vector * v);

/**
 * sums the squared values
 */
double calc_vector_squaresum(gsl_vector * v);
/**
 * returns a duplicate.
 *
 * The caller has to free the returned vector.
 */
gsl_vector * dup_vector(gsl_vector * v);

/**
 * normalizes the vector, i.e. the values are scaled so that the sum of
 * all values is 1

 * The caller has to free the returned vector.
 */
gsl_vector * calc_normalized(gsl_vector * v);


#endif
