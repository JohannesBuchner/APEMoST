/**
 * Number of chains to use for parallel tempering
 *
 * <= 12 is not recommended
 */
#ifndef N_BETA
#define N_BETA 20
#endif
/**
 * Position where the smallest/hottest chain should start.
 * Remember, beta = 1/T, so beta = 0 is infinitely hot.
 *
 * Has to be between 0 and 1. If < 0, it is determined automatically
 * so that the hottest chain will aim at a stepwidth of 1/3 of
 * parameter space.
 */
#ifndef BETA_0
#define BETA_0 -0.001
#endif
/**
 * How many iterations should be used for the burn-in phase?
 */
#ifndef BURN_IN_ITERATIONS
#define BURN_IN_ITERATIONS 10000
#endif
/**
 * How much deviation from the desired acceptance rate is acceptable
 */
#ifndef RAT_LIMIT
#define RAT_LIMIT 0.5
#endif
/**
 * How many iterations should be used after burn-in for adjusting
 * the step widths?
 */
#ifndef ITER_LIMIT
#define ITER_LIMIT 100000
#endif
/**
 * Factor used for scaling the step widths
 */
#ifndef MUL
#define MUL 0.85
#endif
/**
 * After how many iterations should a swap occur?
 *
 * Performance evaluations suggest to set this to 2000/N_BETA.
 * If < 0, this will be done for you.
 */
#ifndef N_SWAP
#define N_SWAP -30
#endif


#ifndef PARAMS_FILENAME
#define PARAMS_FILENAME "params"
#endif
#ifndef DATA_FILENAME
#define DATA_FILENAME "data"
#endif

