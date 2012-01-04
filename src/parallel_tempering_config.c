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

#include <omp.h>

#include "mcmc.h"
#include "parallel_tempering_config.h"
#include "debug.h"
#include "define_defaults.h"
#include "gsl_helper.h"
#include "utils.h"

void write_params_file(mcmc * m) {
	unsigned int i;
	FILE * f = fopen(PARAMS_FILENAME "_suggested", "w");
	if (f != NULL) {
		for (i = 0; i < get_n_par(m); i++) {
			fprintf(
					f,
					DUMP_FORMAT "\t" DUMP_FORMAT "\t" DUMP_FORMAT "\t%s\t" DUMP_FORMAT "\n",
					gsl_vector_get(get_params_best(m), i), gsl_vector_get(
							get_params_min(m), i), gsl_vector_get(
							get_params_max(m), i), get_params_descr(m)[i],
					gsl_vector_get(get_steps(m), i));
		}
		fclose(f);
		printf("new suggested parameters file has been written\n");
	} else {
		fprintf(stderr,
				"Could not write to file " PARAMS_FILENAME "_suggested\n");
	}
}

void write_calibration_summary(mcmc ** chains, unsigned int n_chains) {
	unsigned int i;
	unsigned int j;
	double beta_0 = get_beta(chains[n_chains - 1]);
	unsigned int n_pars = get_n_par(chains[0]);
	gsl_vector * steps;

	FILE * f = fopen("calibration_summary", "w");
	if (f != NULL) {
		fprintf(f, "Summary of calibrations\n");
		fprintf(f, "\nBETA TABLE\n");
		fprintf(f, "Chain # | Calculated | Calibrated\n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "Chain %d | " DUMP_FORMAT " | %f\n", i, get_chain_beta(
					i, n_chains, beta_0), get_beta(chains[i]));
		}
		fprintf(f, "\nSTEPWIDTH TABLE\n");
		fprintf(f, "Chain # | Calibrated stepwidths... \n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "%d", i);
			for (j = 0; j < n_pars; j++) {
				fprintf(f, "\t" DUMP_FORMAT, get_steps_for(chains[i], j));
			}
			fprintf(f, "\n");
		}
		fprintf(f, "\nSTEPWIDTH ESTIMATE TABLE\n");
		fprintf(f, "If you find that the estimate deviates much or "
			"systematically from the ");
		fprintf(f, "calibrated stepwidths, please notify the authors.\n");
		fprintf(f, "Chain # | Calculated stepwidths... \n");
		for (i = 0; i < n_chains; i++) {
			fprintf(f, "%d", i);
			steps = dup_vector(get_steps(chains[0]));
			gsl_vector_scale(steps, pow(get_beta(chains[i]), -0.5));
			for (j = 0; j < n_pars; j++) {
				fprintf(f, "\t" DUMP_FORMAT, gsl_vector_get(steps, j));
			}
			gsl_vector_free(steps);
			fprintf(f, "\n");
		}
		fclose(f);
		printf("calibration summary has been written\n");
	} else {
		fprintf(stderr, "Could not write to file calibration_summary\n");
	}
}
mcmc ** setup_chains() {
	unsigned int i;
	mcmc ** chains;
	const char * params_filename = PARAMS_FILENAME;
	const char * data_filename = DATA_FILENAME;
	const unsigned int n_beta = N_BETA;
	chains = (mcmc**) mem_calloc(n_beta, sizeof(mcmc*));
	assert(chains != NULL);

	printf("Initializing %d chains\n", n_beta);
	for (i = 0; i < n_beta; i++) {
		chains[i] = mcmc_load_params(params_filename);
		if (i == 0) {
			mcmc_load_data(chains[i], data_filename);
			chains[i]->additional_data
					= mem_malloc(sizeof(parallel_tempering_mcmc));
			set_beta(chains[i], 1.0);
		} else {
			mcmc_reuse_data(chains[i], chains[0]);
		}
		mcmc_check(chains[i]);
		chains[i]->additional_data = mem_malloc(
				sizeof(parallel_tempering_mcmc));
		set_beta(chains[i], 1);
	}
	return chains;
}

/**
 * read betas, stepwidths and start values of all chains
 *
 * @return lines read
 **/
void read_calibration_file(mcmc ** chains, unsigned int n_chains) {
	unsigned int i;
	unsigned int j;
	unsigned int err = 0;
	unsigned int n_par = get_n_par(chains[0]);
	double v;
	FILE * f;

	f = fopen(CALIBRATION_FILE, "r");
	if (f == NULL) {
		fprintf(f, "could not read calibration file %s\n", CALIBRATION_FILE);
		exit(1);
	}

	for (i = 0; i < n_chains && !feof(f); i++) {
		if (fscanf(f, "%lf", &v) != 1)
			err++;
		set_beta(chains[i], v);
		for (j = 0; j < n_par && !feof(f) && err == 0; j++) {
			if (fscanf(f, "%lf", &v) != 1)
				err++;
			set_steps_for(chains[i], v, j);
		}
		for (j = 0; j < n_par && !feof(f) && err == 0; j++) {
			if (fscanf(f, "%lf", &v) != 1)
				err++;
			set_params_for(chains[i], v, j);
		}
		if (feof(f) && j < n_par) {
			fprintf(f, "could not read %d chain calibrations. \nError with "
				"line %d.\n", n_chains, i + 1);
			exit(1);
		}
		set_params_best(chains[i], get_params(chains[i]));
	}

	fclose(f);
}

void write_calibrations_file(mcmc ** chains, const unsigned int n_chains) {
	FILE * f;
	unsigned int i;
	unsigned int j;
	unsigned int n_par = get_n_par(chains[0]);

	f = fopen(CALIBRATION_FILE, "w");
	if (f == NULL) {
		perror("error writing to calibration results file");
		exit(1);
	}
	for (j = 0; j < n_chains; j++) {
		fprintf(f, DUMP_FORMAT, get_beta(chains[j]));
		for (i = 0; i < n_par; i++) {
			fprintf(f, "\t" DUMP_FORMAT, get_steps_for(chains[j], i));
		}
		for (i = 0; i < n_par; i++) {
			fprintf(f, "\t" DUMP_FORMAT, get_params_for(chains[j], i));
		}
		fprintf(f, "\n");
	}
	fclose(f);
	printf("wrote calibration results for %d chains to %s\n", n_chains,
			CALIBRATION_FILE);
}


