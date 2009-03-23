/**
 * Public methods of the class
 */

#ifndef MCMC
#define MCMC

#include <stdio.h>
#include <stdlib.h>

#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif

#include "mcmc_struct.h"

/**
 * create and initialize a mcmc class using the configuration given in
 * @param filename
 */
mcmc * mcmc_load(const char * filename);
/**
 * frees the memory used by the class
 */
void mcmc_free(mcmc * m);

/**
 * adds the current parameter values to params_distr as nth iteration
 * was: add_values
 * @param m
 * @param n_iter
 */
void mcmc_append_current_parameters(mcmc * m, int n_iter);

void mcmc_dump_model(mcmc * m);
void mcmc_dump_y_dat(mcmc * m);

/**
 * write probability/distribution (params_distr) out to files.
 *
 * The files are named after the parameter names, with .prob.dump appended.
 * You can use those files to build a histogram of where the algorithm has been
 * (and how often) to find out
 *
 * @param m
 * @param n_values use the n last iterations. if negative, all are used.
 */
void mcmc_dump_probabilities(mcmc * m, int n_values);

/**
 * check if a new best value has been found
 */
void mcmc_check_best(mcmc * m);


#include "mcmc_markov_chain.h"

#include "mcmc_gettersetter.h"

#endif

