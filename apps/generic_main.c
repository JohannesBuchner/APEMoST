/*
    APEMoST - Automated Parameter Estimation and Model Selection Toolkit
    Copyright (C) 2009  Johannes Buchner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <unistd.h>
#include <string.h>
#include "mcmc.h"
#include "debug.h"
#include "parallel_tempering.h"
#include "parallel_tempering_interaction.h"

/**
 * \mainpage
 *
 * This is the code documentation of APEMoST, the Automated Parameter 
 * Estimation and Model Selection Toolkit.
 *
 * You can find more information about it <a href="../../index.html">in
 * its documentation</a>.
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
 * <li>#TARGET_ACCEPTANCE_RATE</li>
 * <li>#MAX_AR_DEVIATION</li>
 * <li>#ITER_LIMIT</li>
 * <li>#MUL</li>
 * <li>#N_SWAP</li>
 * <li>#N_PARAMETERS</li>
 * <li>#SKIP_CALIBRATE_ALLCHAINS</li>
 * <li>#PROPOSAL</li>
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
 * </ul>
 * \subsection Analyzing
 * <ul>
 * <li>#GNUPLOT_STYLE</li>
 * <li>#HISTOGRAMS_MINMAX</li>
 * </ul>
 * \subsection Others
 * <ul>
 * <li>#DEBUG</li>
 * <li>#VERBOSE</li>
 * <li>#SEGV</li>
 * <li>#NOASSERT</li>
 * <li>#WITHOUT_GARBAGE_COLLECTOR</li>
 * </ul>
 *
 * \section run-params Runtime parameters
 * At runtime, the program looks for the files #PARAMS_FILENAME and
 * #DATA_FILENAME.
 * 
 * Really, read the manual.
 */

#include "define_defaults.h"

#ifndef MAX_ITERATIONS
/**
 * set the number of iterations after you want the program to terminate.
 *
 * This is especially useful in benchmarking.
 * Example: Set this to 100000.
 * 
 * 0 means run indefinitely
 */
#define MAX_ITERATIONS 0
#endif

char * progname;
void usage();
void help_phase(char * phase);
void check();

int main(int argc, char ** argv) {
	progname = argv[0];
	if (argc > 1) {
		if (0 == strcmp(argv[1], "help") || 0 == strcmp(argv[1], "-h")) {
			if (argc == 3)
				help_phase(argv[2]);
			else
				usage();
		} else if (0 == strcmp(argv[1], "check")) {
			check();
		} else if (0 == strcmp(argv[1], "calibrate_first")) {
			calibrate_first();
		} else if (0 == strcmp(argv[1], "calibrate_rest")) {
			calibrate_rest();
		} else if (0 == strcmp(argv[1], "run")) {
			if (argc == 3 && strcmp(argv[2], "--append") == 0)
				prepare_and_run_sampler(MAX_ITERATIONS, 1);
			else if (argc == 2)
				prepare_and_run_sampler(MAX_ITERATIONS, 0);
			else {
				fprintf(stderr, "You are doing it wrong.\n");
				fprintf(stderr, "Did you want to write --append?\n");
				usage();
			}
		} else if (0 == strcmp(argv[1], "analyse")) {
			if (argc == 3 && strcmp(argv[2], "marginal") == 0)
				analyse_marginal_distributions();
			else if (argc == 3 && strcmp(argv[2], "model") == 0)
				analyse_marginal_distributions();
			else if (argc == 2) {
				analyse_marginal_distributions();
				analyse_data_probability();
			} else {
				fprintf(stderr, "You are doing it wrong.\n");
				fprintf(stderr, "Did you want to write 'analyse marginal' or "
					"'analyse model'?\n");
				usage();
			}
		} else {
			fprintf(stderr, "You are doing it wrong.\n");
			usage();
		}
	} else {
		fprintf(stderr, "No phase specified.\n");
		usage();
	}
	return 0;
}

void usage() {
	fprintf(stderr, "SYNAPSIS: %s <phase> <...>\n\n", progname);
	fprintf(stderr, "\t-h, help\tthis clutter\n"
		"\tphase is one of: \n"
		"\t\tcheck          \toutput which parameters and files will be used\n"
		"\t\t               \tand check if they are there\n"
		"\t\tcalibrate_first\tcalibrate first chain (beta = 1)\n"
		"\t\tcalibrate_rest \tcalibrate remaining chains (beta < 1)\n"
		"\t\trun [--append] \tcreate and dump sampling data\n"
		"\t\t               \twithout adding --append, existing data is overwritten\n"
		"\t\tanalyse        \tanalyse the available data probability\n"
		"\t\thelp <phase>   \tprint more information about a phase\n"
		"\n");
	fprintf(stderr,
			"Read the manual on how to setup a working directory and \n"
				"what parameters can be set.\n\n"
				"Remember: Scientific simulations should be repeatable, unbiased, \n"
				"          statistically sound and rigorous.\n"
			"\n"
			"APEMoST  Copyright (C) 2009  Johannes Buchner\n"
			"This program comes with ABSOLUTELY NO WARRANTY; for details see the file LICENSE.\n"
			"This is free software, and you are welcome to redistribute it\n"
			"under certain conditions; see the file LICENSE.\n");
}

void help_phase(char * phase) {
	if (0 == strcmp(phase, "check")) {
		printf("Phase 'check'\n\n"
			"Prerequisites: \n"
			"\tnone\n"
			"Provides: \n"
			"Does:\n"
			"\tChecks if data and parameter files are available\n"
			"\tPrints compile information such as defined algorithm "
			"behaviour)\n");
	} else if (0 == strcmp(phase, "calibrate_first")) {
		printf("Phase 'calibrate_first'\n\n"
			"Prerequisites: \n"
			"\tparameters file " PARAMS_FILENAME "\n"
		"\tdata file " DATA_FILENAME "\n"
		"\tBURN_IN_ITERATIONS\n"
		"Provides: \n"
		"\tstepwidths for chain 0 (beta = 1)\n"
		"\tstart parameters for chain 0 (beta = 1)\n"
		"Does:\n"
		"\tBurn-in and step-width calibration for first chain\n"
		"\twrites the resulting stepwidths and current parameter values\n"
		"\t\tout to the file " CALIBRATION_FILE " and a new params file\n"
		"\t\tat " PARAMS_FILENAME "_suggested\n"
		"\n");
	} else if (0 == strcmp(phase, "calibrate_rest")) {
		printf("Phase 'calibrate_rest'\n\n"
			"Prerequisites: \n"
			"\tparameters file " PARAMS_FILENAME "\n"
		"\tdata file " DATA_FILENAME "\n"
		"\tBURN_IN_ITERATIONS\n"
		"\tBETA_0\n"
		"\tN_BETA\n"
		"\tBETA_ALIGNMENT\n"
		"\tstepwidths and start parameters of chain 0 (beta = 1)\n"
		"\t\t(read from first line of file "CALIBRATION_FILE ")\n");
		printf("Provides: \n"
			"\tstepwidths for all chains\n"
			"\tstart parameters for chains\n"
			"\tbeta values for all chains\n"
			"Does:\n"
			"\tBurn-in and step-width calibration for the remaining chains (beta < 1)\n"
			"\twrites the resulting stepwidths, betas and current parameter values\n"
			"\t\tout to the file " CALIBRATION_FILE "\n"
		"\n");
	} else if (0 == strcmp(phase, "run")) {
		printf(
				"Phase 'run'\n\n"
					"Prerequisites: \n"
					"\tparameters file " PARAMS_FILENAME "\n"
				"\tdata file " DATA_FILENAME "\n"
				"\tN_BETA\n"
				"\tstepwidths, betas and start parameters of all chains\n"
				"\t\t(read from file "CALIBRATION_FILE ")\n"
				"Provides: \n"
				"\tdata dumps (probabilities and visited parameters)\n"
				"Does:\n"
				"\tRun the MCMC engine and continously write out the data of\n"
				"\tthe visited parameters and probabilities to files. \n"
				"Options:\n"
				"\t--append\tcauses to append to the existing dump files rather than overwrite\n"
				"\n");
	} else if (0 == strcmp(phase, "analyse")) {
		printf("Phase 'analyse'\n\n"
			"Prerequisites: \n"
			"\tparameters file " PARAMS_FILENAME "\n"
		"\tdata file " DATA_FILENAME "\n"
		"\tN_BETA\n"
		"\tbetas of all chains\n"
		"\t\t(read from file "CALIBRATION_FILE ")\n"
		"\tprobabilities data dumps\n"
		"\tvisited values data dumps\n"
		"Provides: \n"
		"\tmodel probability\n"
		"Does:\n"
		"\tAnalyses the model probability of the visited values.\n"
		"\tThe resulting value can be used to compare to other models.\n"
		"\tAnalyses the marginal distribution by making a histograms of all\n"
		"\tvisited values.\n");
		printf("Options:\n"
			"\tmarginal\tcalculate marginal distribution only\n"
			"\tmodel\tcalculate model probability only\n"
			"\n");
	} else {
		printf("No help available for unknown phase.\n");
		usage();
	}
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
	OUTPUT_PARAMD(TARGET_ACCEPTANCE_RATE);
	OUTPUT_PARAMD(MAX_AR_DEVIATION);
	OUTPUT_PARAMI(ITER_LIMIT);
	OUTPUT_PARAMD(MUL);
	OUTPUT_PARAMI(N_SWAP);
#ifdef SKIP_CALIBRATE_ALLCHAINS
	printf("\tSKIP_CALIBRATE_ALLCHAINS: enabled (calibrating only 2 chains)\n");
#else
	printf("\tSKIP_CALIBRATE_ALLCHAINS: disabled (calibrating all chains)\n");
#endif
#ifdef PROPOSAL_UNIFORM
	printf("\tPROPOSAL: uniform proposal distribution\n");
#elif defined PROPOSAL_LOGISTIC
	printf("\tPROPOSAL: logistic proposal distribution\n");
#else
	printf("\tPROPOSAL: gaussian/normal proposal distribution\n");
#endif

#ifdef N_PARAMETERS
	printf("\tN_PARAMETERS: Compiled for fixed number of parameters: %d\n", N_PARAMETERS);
#else
	printf("\tN_PARAMETERS: Compiled for variable number of parameters\n");
#endif

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
#ifdef WITHOUT_GARBAGE_COLLECTOR
	printf("WITHOUT_GARBAGE_COLLECTOR: no garbage collector\n");
#else
	printf("WITHOUT_GARBAGE_COLLECTOR: garbage collector enabled\n");
#endif

}
