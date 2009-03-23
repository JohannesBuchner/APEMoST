#ifndef MCMC_MARKOV_CHAIN_H_
#define MCMC_MARKOV_CHAIN_H_

#include "mcmc.h"

#define DEFAULT_MUL 0.85
#define DEFAULT_ADJUST_STEP 0.5
#define DEFAULT_RAT_LIMIT -1 /* is calculated then */

/**
 * create/calibrate the markov-chain
 */
void markov_chain_calibrate(mcmc * m, unsigned int burn_in_iterations, double rat_limit,
		unsigned int iter_limit, double mul, double adjust_step);
/**
 * take a step using the markov-chain
 */
void markov_chain_step(mcmc * m);

#endif /* MCMC_MARKOV_CHAIN_H_ */
