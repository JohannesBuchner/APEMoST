/**
 * Here are the "private" methods of the class and helper functions
 */
#ifndef MCMC_INTERNAL_H_
#define MCMC_INTERNAL_H_

#include "mcmc.h"
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf.h>

/**
 * create class
 * \private
 * @param n_pars parameters
 */
mcmc * mcmc_init(const unsigned int n_pars);

/**
 * prepare the calculation of (next) iteration, i.e., allocate space
 * @param m
 * @param iter index of new iteration
 */
void mcmc_prepare_iteration(mcmc * m, const unsigned long iter);

/**
 * @see calc_hist
 */
gsl_histogram * get_hist(mcmc * m, int i, int nbins);

/**
 * count the lines (\n) in the file
 * @param filename
 */
unsigned int countlines(const char * filename);

/**
 * a modulo operator for double values
 */
double mod_double(const double x, const double div);

#endif /* MCMC_INTERNAL_H_ */
