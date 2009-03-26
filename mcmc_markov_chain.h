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
 * @param calc_index 1 if the model should be recalculated, 0 otherwise
 */
void markov_chain_step(mcmc * m, int calc_index);
/**
 * take a step using the markov-chain for the indexth parameter
 * @param index the param to look at
 * @param calc_index 1 if the model should be recalculated, 0 otherwise
 */
void markov_chain_step_for(mcmc * m, unsigned int index, int calc_index);

#endif /* MCMC_MARKOV_CHAIN_H_ */
