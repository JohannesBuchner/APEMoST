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
#endif

#ifndef DUMP_PROB_LENGTH
/**
 * How many probability values should be dumped?
 * Set to -1 if you want to dump all of them.
 */
#define DUMP_PROB_LENGTH  1000*3
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

void set_beta(mcmc * m, double newbeta);

double get_beta(mcmc * m);

void inc_swapcount(mcmc * m);

unsigned long get_swapcount(mcmc * m);

/**
 * begin a parallel tempering execution
 */
void parallel_tempering(const char * params_filename,
		const char * data_filename, int n_beta, double beta_0,
		unsigned long burn_in_iterations, double rat_limit,
		unsigned long iter_limit, double mul, int n_swap);

#endif /* PARALLEL_TEMPERING_H_ */
