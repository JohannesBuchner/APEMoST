#include <unistd.h>
#include <string.h>
#include "mcmc.h"
#include "debug.h"
#include "parallel_tempering.h"
#include "tempering_interaction.h"

/**
 * \mainpage
 *
 * This is the documentation of mcmc, a parameter fitting program.
 *
 * It uses markov chain and Michael Gruberbauers most sophisticated
 * algorithms to check a model with a number of parameters against
 * observed data.
 * You can read more about the theory behind it in
 * Gregory, "Bayesian Logical Data Analysis for the Physical Sciences".
 *
 * \section Compiling
 * See the file INSTALL
 *
 * \section cc-params Compile-time parameters
 * Instead of adjusting parameters in files (which would cause a lot of
 * unnecessary versions), you should pass the parameters to the compiler.
 * This includes decisions on which algorithm to use, as well as fine-
 * tuning values like number of burn-in iterations.
 *
 * Do it like this:<br>
 * <code>$ CCFLAGS="-DN_PARAMETERS=3 -DDEBUG" make simplesin.exe</code>
 *
 * The possible parameters and their values:
 * \subsection tuning Tuning the algorithm-related
 * <ul>
 * <li style="color:red">#BETA_DISTRIBUTION (required)</li>
 * <li>#N_BETA</li>
 * <li>#BETA_0</li>
 * <li>#BURN_IN_ITERATIONS</li>
 * <li>#RAT_LIMIT</li>
 * <li>#ITER_LIMIT</li>
 * <li>#MUL</li>
 * <li>#N_SWAP</li>
 * <li>#SIGMA</li>
 * <li>#N_PARAMETERS</li>
 * </ul>
 * \subsection select Selecting a algorithm
 * <ul>
 * <li>#RANDOMSWAP</li>
 * <li>#RESET_TO_BEST</li>
 * </ul>
 * \subsection Running
 * <ul>
 * <li>#MAX_ITERATIONS</li>
 * <li>#PRINT_PROB_INTERVAL</li>
 * <li>#DUMP_PROB_LENGTH</li>
 * <li>#PRINT_PROB_INTERVAL</li>
 * </ul>
 * \subsection Others
 * <ul>
 * <li>#DEBUG</li>
 * <li>#VERBOSE</li>
 * <li>#SEGV</li>
 * <li>#NOASSERT</li>
 * </ul>
 *
 * \section run-params Runtime parameters
 * At runtime, the program looks for the files #PARAMS_FILENAME and
 * #DATA_FILENAME.
 * #PARAMS_FILENAME should contain the following columns (tab seperated):
 * start-value, minimum, maximum, name, step-size.
 * Example: <br>
 * <code>0.7	0.4	3.0   Amplitude 0.01</code>
 *
 * The file #DATA_FILENAME should contain two columns (x/y) values.
 *
 * \section Results
 * So the program runs and shows you where your chains are. But how do you
 * get the visited datapoints in parameter space?
 * For performance reasons, these are not dumped all the time.
 * The are written out on exit (or termination with Ctrl-C).
 *
 * You can also send the process a signal to dump the datapoints: <br>
 * <code>$ kill -SIGUSR1 processid</code>
 *
 * If your program is called "simplesin.exe", you can do
 *
 * <code>$ kill -SIGUSR1 $(pidof simplesin.exe)</code>
 *
 * @see DUMP_PROB_LENGTH
 *
 */

/**
 * Number of chains to use for parallel tempering
 */
#ifndef N_BETA
#define N_BETA 12/3
#endif
/**
 * Position where the smallest/hottest chain should start.
 * Remember, beta = 1/T, so beta = 0 is infinite heat.
 *
 * Has to be between 0 and 1
 */
#ifndef BETA_0
#define BETA_0 0.001
#endif
/**
 * How many iterations should be used for the burn-in phase?
 */
#ifndef BURN_IN_ITERATIONS
#define BURN_IN_ITERATIONS 10000
#endif
/**
 * TODO: document
 */
#ifndef RAT_LIMIT
#define RAT_LIMIT 0.5
#endif
/**
 * How many iterations should be used after burn-in for adjusting
 * the step widths?
 */
#ifndef ITER_LIMIT
#define ITER_LIMIT 20000
#endif
/**
 * Factor used for scaling the step widths
 */
#ifndef MUL
#define MUL 0.85
#endif
/**
 * After how many iterations should a swap occur?
 */
#ifndef N_SWAP
#define N_SWAP 30
#endif
/**
 * TODO: document
 */
#ifndef SIGMA
#define SIGMA 0.5
#endif

#define PARAMS_FILENAME "params"
#define DATA_FILENAME "data"

void usage(const char * progname);
void check(const char * progname);

int main(int argc, char ** argv) {
	if (argc > 1) {
		if (0 == strcmp(argv[1], "--help") || 0 == strcmp(argv[1], "-h")) {
			usage(argv[0]);
		} else if (0 == strcmp(argv[1], "--check")) {
			check(argv[0]);
		} else {
			fprintf(stderr, "You are doing it wrong.\n");
			usage(argv[0]);
		}
	} else {
		parallel_tempering(PARAMS_FILENAME, DATA_FILENAME, N_BETA, BETA_0,
				BURN_IN_ITERATIONS, RAT_LIMIT, ITER_LIMIT, MUL, N_SWAP);
	}
	return 0;
}

void usage(const char * progname) {
	fprintf(stderr, "%s: SYNAPSIS\n\n"
		"\t-h, --help\tthis clutter\n"
		"\t--check\toutput which parameters and files will be used\n"
		"\t       \tand check if they are there\n"
		"\n"
		"Read the manual on how to setup a working directory and \n"
		"what parameters can be set.\n", progname);
}

void checkfile(char * filename) {
	printf("\t%s\t%s\n", filename, access(filename, R_OK) == 0 ? "found"
			: "not a readable file");
}

#define OUTPUT_PARAMD(P) printf("\t%s: %f\n", #P, P);
#define OUTPUT_PARAMI(P) printf("\t%s: %d\n", #P, P);

void check(const char * progname) {
	printf("%s: Checking environment:\n\n", progname);
	printf("Compile Parameters:\n");
	printf("\tBETA_DISTRIBUTION: %s\n", TOSTRING(BETA_DISTRIBUTION));
	OUTPUT_PARAMI(N_BETA);
	OUTPUT_PARAMD(BETA_0);
	OUTPUT_PARAMI(BURN_IN_ITERATIONS);
	OUTPUT_PARAMD(RAT_LIMIT);
	OUTPUT_PARAMI(ITER_LIMIT);
	OUTPUT_PARAMD(MUL);
	OUTPUT_PARAMI(N_SWAP);
	OUTPUT_PARAMD(SIGMA);

	printf("\nFiles:\n");
	checkfile(PARAMS_FILENAME);
	checkfile(DATA_FILENAME);

	printf("\nDebugging Parameters:\n");
#ifdef DEBUG
	printf("\tDebug output: on\n");
#else
	printf("\tDebug output: off\n");
#endif
#ifdef VERBOSE
	printf("\tVerbose debug output: on\n");
#else
	printf("\tVerbose debug output: off\n");
#endif
#ifdef SEGV
	printf("\tSEGFAULT output: on\n");
#else
	printf("\tSEGFAULT output: off\n");
#endif
#ifdef NOASSERT
	printf("\tRuntime sanity checks: off\n");
#else
	printf("\tRuntime sanity checks: on\n");
#endif
#ifdef N_PARAMETERS
	printf("\tCompiled for fixed number of parameters: %d\n", N_PARAMETERS);
#else
	printf("\tCompiled for variable number of parameters\n");
#endif

}
