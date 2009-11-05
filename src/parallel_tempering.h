#ifndef PARALLEL_TEMPERING_H_
#define PARALLEL_TEMPERING_H_

#include "mcmc.h"
#include "parallel_tempering_beta.h"

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * set the number of iterations after you want the program to terminate.
 *
 * This is especially useful in benchmarking.
 * Example: Set this to 100000.
 */
#define MAX_ITERATIONS
/**
 * should all chains be dumped?
 *
 * Otherwise, only chain0 is dumped (beta = 1)
 */
#define DUMP_ALL_CHAINS
/**
 * Enabling this will only calibrate the first two chains, not all of them.
 * The stepwidths for the rest of the chains will be predicted.
 */
#define SKIP_CALIBRATE_ALLCHAINS
#endif

#ifndef PRINT_PROB_INTERVAL
/**
 * After how many iterations do you want the program to write out
 * information?
 *
 * Note that this should be a multiple of N_SWAP, otherwise you get
 * weird effects, since this is checked after N_SWAP iterations using
 * <code>(iter % PRINT_PROB_INTERVAL == 0)</code>.
 */
#define PRINT_PROB_INTERVAL 1000
#endif

#define CALIBRATION_FILE "calibration_results"

/** applications can run the follwing functions */

void calibrate_first();

void prepare_and_run_sampler(int append);

void calibrate_rest();

void analyse_marginal_distributions();

void analyse_data_probability();

#endif /* PARALLEL_TEMPERING_H_ */
