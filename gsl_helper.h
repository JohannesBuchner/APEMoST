#ifndef MCMC_GSL_HELPER
#define MCMC_GSL_HELPER

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>



/**
 * calculate a histogram
 */
gsl_histogram * calc_hist(gsl_vector ** vs, int index, int nbins);

/**
 * sums the values
 */
double calc_vector_sum(gsl_vector * v);

/**
 * returns a duplicate
 */
gsl_vector * dup_vector(gsl_vector * v);

/**
 * normalizes the vector, i.e. the values are scaled so that the sum of 
 * all values is 1
 */
gsl_vector * calc_normalized(gsl_vector * v);


#endif