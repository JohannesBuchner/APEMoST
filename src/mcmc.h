/**
 * Public methods of the class
 */

#ifndef MCMC
#define MCMC

#include <stdio.h>
#include <stdlib.h>
#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Turns additional sanity checks off.
 *
 * One could think that turning this off would improve performance,
 * but tests have shown the modification is not significant.
 */
#define NOASSERT
#endif

#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif

#include "mcmc_struct.h"

/**
 * create and initialize a mcmc class using the configuration given in
 * @param filename  file containing the model parameters
 * @param datafilename x/y observed data
 * @return the created mcmc class
 */
mcmc * mcmc_load(const char * filename, const char * datafilename);

/**
 * create and initialize a mcmc class using the configuration given in
 * @param filename  file containing the model parameters
 * @return the created mcmc class
 */
mcmc * mcmc_load_params(const char * filename);

/**
 * loads the data from the given file as x/y values
 * @param m
 * @param datafilename
 */
void mcmc_load_data(mcmc * m, const char * datafilename);

/**
 * reference to another object for x and y-data.
 *
 * @param m the object to fill
 * @param m_orig the object with loaded data
 */
void mcmc_reuse_data(mcmc * m, const mcmc * m_orig);
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
 * @param suffix added to output filename as a suffix
 */
void mcmc_dump_probabilities(const mcmc * m, int n_values, const char * suffix);

/**
 * check if a new best value has been found
 * @param m
 */
void mcmc_check_best(mcmc * m);

#include "markov_chain.h"

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

