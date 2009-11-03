#ifndef MCMC_MARKOV_CHAIN_H_
#define MCMC_MARKOV_CHAIN_H_

#include "mcmc.h"
#include "define_defaults.h"

#define DEFAULT_ADJUST_STEP 0.5
#ifndef NO_RESCALING_LIMIT
#define NO_RESCALING_LIMIT 15
#endif

#ifndef ITER_READJUST
#define ITER_READJUST 200
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
void markov_chain_calibrate(mcmc * m, const unsigned int burn_in_iterations,
		double desired_acceptance_rate, const double max_ar_deviation,
		const unsigned int iter_limit, double mul, const double adjust_step);
/**
 * take a step using the markov-chain
 * @param m
 */
void markov_chain_step(mcmc * m);
/**
 * take a step using the markov-chain for the indexth parameter
 * @param m
 * @param index the param to look at
 */
void markov_chain_step_for(mcmc * m, const unsigned int index);

/**
 * adapts the step width
 */
void rmw_adapt_stepwidth(mcmc * m, double prob_old);

void burn_in(mcmc * m, const unsigned int burn_in_iterations);
#endif /* MCMC_MARKOV_CHAIN_H_ */
