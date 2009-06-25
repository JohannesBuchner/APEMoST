#ifndef PARALLEL_TEMPERING_H_
#define PARALLEL_TEMPERING_H_

#include "mcmc.h"

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

typedef struct {
	/**
	 * likelihood of acceptance
	 */
	double beta;

	/**
	 * times this was swapped
	 */
	unsigned long swapcount;

} parallel_tempering_mcmc;

void set_beta(mcmc * m, const double newbeta);

double get_beta(const mcmc * m);

void inc_swapcount(mcmc * m);

unsigned long get_swapcount(const mcmc * m);

/**
 * begin a parallel tempering execution
 */
void parallel_tempering(const char * params_filename,
		const char * data_filename, const int n_beta, const double beta_0,
		const unsigned long burn_in_iterations, const double rat_limit,
		const unsigned long iter_limit, const double mul, const int n_swap);

#endif /* PARALLEL_TEMPERING_H_ */
