#ifndef _PARALLEL_TEMPERING_BETA_H
#define _PARALLEL_TEMPERING_BETA_H
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "parallel_tempering.h"
#include "parallel_tempering_interaction.h"

void set_beta(mcmc * m, double newbeta);
double get_beta(const mcmc * m);
void inc_swapcount(mcmc * m);
unsigned long get_swapcount(const mcmc * m);
void print_current_positions(const mcmc ** sinmod, const int n_beta);

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Defines how the beta value should be assigned/distributed between the chains.
 *
 * beta = 1 / temperature.
 *
 * You can choose equidistant (linear) alignment of the temperature or beta.
 * Or, what often proves to be a good choice, you can use Chebyshev nodes
 * for the values of the temperature or beta.
 *
 * e.g.: BETA_ALIGNMENT=equidistant_beta <br>
 * also available: equidistant_temperature, chebyshev_beta,
 * chebyshev_temperature, equidistant_stepwidth, chebyshev_stepwidth
 *
 * default: chebyshev_beta
 */
#define BETA_ALIGNMENT
#endif

#ifndef BETA_ALIGNMENT
#define BETA_ALIGNMENT chebyshev_beta
#endif

double get_chain_beta(unsigned int i, unsigned int n_beta, double beta_0);

/**
 * hottest chain stepwidth in units of parameter space
 */
#define BETA_0_STEPWIDTH 1.0

double calc_beta_0(mcmc * m, gsl_vector * stepwidth_factors);

#endif

