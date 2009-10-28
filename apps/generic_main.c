#include <unistd.h>
#include <string.h>
#include "mcmc.h"
#include "debug.h"
#include "parallel_tempering.h"
#include "parallel_tempering_interaction.h"

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
 * \subsection tuning Fine-tuning the algorithm
 * <ul>
 * <li>#BETA_ALIGNMENT</li>
 * <li>#N_BETA</li>
 * <li>#BETA_0</li>
 * <li>#BURN_IN_ITERATIONS</li>
 * <li>#RAT_LIMIT</li>
 * <li>#ITER_LIMIT</li>
 * <li>#MUL</li>
 * <li>#N_SWAP</li>
 * <li>#N_PARAMETERS</li>
 * <li>#SKIP_CALIBRATE_ALLCHAINS</li>
 * </ul>
 * \subsection alg Defining algorithm behaviour
 * <ul>
 * <li>#RANDOMSWAP</li>
 * <li>#ADAPT</li>
 * <li>#RWM</li>
 * </ul>
 * \subsection Running
 * <ul>
 * <li>#MAX_ITERATIONS</li>
 * <li>#DUMP_ALL_CHAINS</li>
 * <li>#PRINT_PROB_INTERVAL</li>
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
 */

#include "define_defaults.h"
char * progname;
void usage();
void check();

int main(int argc, char ** argv) {
	if (argc > 1) {
		progname = argv[0];
		if (0 == strcmp(argv[1], "help") || 0 == strcmp(argv[1], "-h")) {
			usage();
		} else if (0 == strcmp(argv[1], "check")) {
			check();
		} else if (0 == strcmp(argv[1], "calibrate_first")) {
			calibrate_first();
		} else if (0 == strcmp(argv[1], "calibrate_rest")) {
			calibrate_rest();
		} else if (0 == strcmp(argv[1], "run")) {
			prepare_and_run_sampler();
		/*} else if (0 == strcmp(argv[1], "analyze")) {
			analyze();*/
		} else {
			fprintf(stderr, "You are doing it wrong.\n");
			usage(argv[0]);
		}
	} else {
		fprintf(stderr, "You are doing it wrong.\n");
		usage(argv[0]);
	}
	return 0;
}

void usage() {
	fprintf(stderr, "%s: SYNAPSIS\n\n"
		"\t-h, help\tthis clutter\n"
		"\tcheck\toutput which parameters and files will be used\n"
		"\t       \tand check if they are there\n"
		"\tcalibrate_first\tcalibrate first chain (beta = 1)\n"
		"\tcalibrate_rest\tcalibrate remaining chains (beta < 1)\n"
		"\trun\tcreate and dump sampling data\n"
		"\tanalyze\tanalyze the available data\n"
		"\t       \tmarginal distribution, data probability"
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

void check() {
	printf("%s: Checking environment:\n", progname);

	printf("\nFiles:\n");
	checkfile(PARAMS_FILENAME);
	checkfile(DATA_FILENAME);

	printf("\nFine-tuning the algorithm:\n");
	printf("\tBETA_ALIGNMENT: %s\n", TOSTRING(BETA_ALIGNMENT));
	OUTPUT_PARAMI(N_BETA);
	OUTPUT_PARAMD(BETA_0);
	OUTPUT_PARAMI(BURN_IN_ITERATIONS);
	OUTPUT_PARAMD(RAT_LIMIT);
	OUTPUT_PARAMI(ITER_LIMIT);
	OUTPUT_PARAMD(MUL);
	OUTPUT_PARAMI(N_SWAP);

	printf("\nDefining algorithm behaviour:\n");
	printf("\tRANDOMSWAP: Random swapping: ");
#ifdef RANDOMSWAP
	printf("on\n");
#else
	printf("off\n");
#endif
	printf("\tRESET_TO_BEST: Resetting to best: ");
#ifdef RESET_TO_BEST
	printf("on, %d\n", RESET_TO_BEST);
#else
	printf("off\n");
#endif
#ifdef MAX_ITERATIONS
	printf("\tMAX_ITERATIONS: Stops after %d iterations\n", MAX_ITERATIONS);
#else
	printf("\tMAX_ITERATIONS: Run indefinitely long\n");
#endif
	OUTPUT_PARAMI(PRINT_PROB_INTERVAL);

	printf("\nDebugging Parameters:\n");
	printf("\tDEBUG: Debug output: ");
#ifdef DEBUG
	printf("on\n");
#else
	printf("off\n");
#endif
	printf("\tVERBOSE: Verbose debug output: ");
#ifdef VERBOSE
	printf("on\n");
#else
	printf("off\n");
#endif
	printf("\tSEGV: Segfault detecting output: ");
#ifdef SEGV
	printf("on\n");
#else
	printf("off\n");
#endif
	printf("\tNOASSERT: Runtime sanity checks: ");
#ifdef NOASSERT
	printf("off\n");
#else
	printf("on\n");
#endif
#ifdef N_PARAMETERS
	printf("\tN_PARAMETERS: Compiled for fixed number of parameters: %d\n", N_PARAMETERS);
#else
	printf("\tN_PARAMETERS: Compiled for variable number of parameters\n");
#endif

}
