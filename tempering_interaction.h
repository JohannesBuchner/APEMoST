#ifndef TEMPERING_INTERACTION_H_
#define TEMPERING_INTERACTION_H_

#include "mcmc.h"

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Enables resetting to the best found value of this chain.
 *
 * A chain is chosen by n_beta * RESET_TO_BEST * urandom() < n_beta
 * Example: Set this to 1000.
 *
 * NB: Do not use, this is cheating the probability process.
 */
#define RESET_TO_BEST
#endif

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * this shouldn't make any difference at the moment ...
 */
#define RANDOMSWAP
#endif

/**
 * does swapping chains and other mixes.
 *
 * @param sinmod
 * @param iter number of iterations that passed. Is a multiple of n_swap.
 * @param n_beta size of sinmod
 */
void tempering_interaction(mcmc ** sinmod, unsigned int n_beta,
		unsigned long iter);

#endif /* TEMPERING_INTERACTION_H_ */
