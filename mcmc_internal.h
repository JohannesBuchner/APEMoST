/**
 * Here are the "private" methods of the class
 */
#ifndef MCMC_INTERNAL_H_
#define MCMC_INTERNAL_H_

#include "mcmc.h"
#include <gsl/gsl_histogram.h>

/**
 * create class
 * \private
 * @param n_pars parameters
 */
mcmc * mcmc_init(unsigned int n_pars);

/**
 *
 */
double get_random_number();

/**
 * prepare the calculation of (next) iteration, i.e., allocate space
 * @param m
 * @param iter number of iteration
 */
void mcmc_prepare_iteration(mcmc * m, unsigned long iter);

/**
 * @see calc_hist
 */
gsl_histogram * get_hist(mcmc * m, int i, int nbins);


#endif /* MCMC_INTERNAL_H_ */
