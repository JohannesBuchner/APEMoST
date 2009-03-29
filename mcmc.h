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
 * @return the created mcmc class
 */
mcmc * mcmc_load(const char * filename);
/**
 * frees the memory used by the class
 *
 * @return NULL for simple assignment <code>x = mcmc_free(x)</code>;
 */
mcmc * mcmc_free(mcmc * m);

/**
 * checks the pointers and dimensions
 */
void mcmc_check(const mcmc * m);

/**
 * adds the current parameter values to params_distr as nth iteration
 * was: add_values
 * @param m
 */
void mcmc_append_current_parameters(mcmc * m);

void mcmc_dump_model(const mcmc * m);
void mcmc_dump_y_dat(const mcmc * m);

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
void mcmc_dump_probabilities(const mcmc * m, int n_values);

/**
 * check if a new best value has been found
 * @param m
 */
void mcmc_check_best(mcmc * m);

#include "mcmc_markov_chain.h"

#include "mcmc_gettersetter.h"

/* calculations done by the application */

/**
 * update the model according to the new parameter values and
 * recalculate the probability for the model
 *
 * @param m
 * @param old_values previous values, or NULL
 */
void calc_model(mcmc * m, const gsl_vector * old_values);
/**
 * update the model as the new parameter value i changed and
 * recalculate the probability for the model
 *
 * @param m
 * @param i index of the parameter value that changed
 * @param old_value previous value of the parameter
 */
void calc_model_for(mcmc * m, const unsigned int i, const double old_value);

#endif

