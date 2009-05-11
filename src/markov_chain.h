#ifndef MCMC_MARKOV_CHAIN_H_
#define MCMC_MARKOV_CHAIN_H_

#include "mcmc.h"

#define DEFAULT_MUL 0.85
#define DEFAULT_ADJUST_STEP 0.5
#define DEFAULT_RAT_LIMIT -1 /* is calculated then */
#define NO_RESCALING_LIMIT 15

#ifndef TARGET_ACCEPTANCE_RATE
#define TARGET_ACCEPTANCE_RATE 0.23
#endif

#ifndef CIRCULAR_PARAMS
/**
 * Which parameters are circular?
 *
 * e.g. if the first and second parameters are angles, you write
 * CIRCULAR_PARAMS=1,2;
 * Do not use zero for the first parameter.
 *
 * N.B.: if you have a angle, you will want to make it go from 0 to 1
 * or similar and use pi in your formula. This way it will be more exact.
 */
#define CIRCULAR_PARAMS 0
#endif

/**
 * create/calibrate the markov-chain
 *
 * @param m
 * @param rat_limit average acceptance rates for individual parameters to be achieved
 * @param burn_in_iterations number of burn-in iterations
 * @param iter_limit number of iterations for step width calibration
 * @param mul factor for adjusting the step width during calibration
 * @param adjust_step gives the factor with which to adjust the stepwidths after burn-in
 */
void markov_chain_calibrate(mcmc * m, unsigned int burn_in_iterations,
		double rat_limit, unsigned int iter_limit, double mul,
		double adjust_step);
/**
 * take a step using the markov-chain
 * @param m
 * @param calc_index 1 if the model should be recalculated, 0 otherwise
 */
void markov_chain_step(mcmc * m, int calc_index);
/**
 * take a step using the markov-chain for the indexth parameter
 * @param m
 * @param index the param to look at
 * @param calc_index 1 if the model should be recalculated, 0 otherwise
 */
void markov_chain_step_for(mcmc * m, unsigned int index, int calc_index);

/**
 * adapts the step width
 */
void rmw_adapt_stepwidth(mcmc * m, double prob_old);

#endif /* MCMC_MARKOV_CHAIN_H_ */
