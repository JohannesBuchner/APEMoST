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
#include "define_defaults.h"

int main(int argc, char ** argv) {
	mcmc * m;
	int r;
	unsigned int i;
	double v;

	if (argc != 1) {
		fprintf(stderr, "SYNOPSIS: %s\n"
			"\n"
			"This program evaluates the model for the parameters read from \n"
			"stdin. Enter the parameters in the same order as defined in %s.\n"
			"The log-posterior probability (and the prior) will be printed.\n"
			"\n"
			"APEMoST  Copyright (C) 2009  Johannes Buchner\n"
			"This program comes with ABSOLUTELY NO WARRANTY; for details see the file LICENSE.\n"
			"This is free software, and you are welcome to redistribute it\n"
			"under certain conditions; see the file LICENSE.\n"
			"", argv[0], PARAMS_FILENAME);
		exit(1);
	}
	m = mcmc_load_params(PARAMS_FILENAME);
	mcmc_load_data(m, DATA_FILENAME);
	mcmc_check(m);
	m->additional_data = mem_malloc(sizeof(parallel_tempering_mcmc));
	set_beta(m, 1.0);
	while (!feof(stdin)) {
		r = 1;
		for (i = 0; i < get_n_par(m) && r == 1; i++) {
			r = scanf("%lf", &v);
			if (r == 1) {
				set_params_for(m, v, i);
				assert(v <= get_params_max_for(m, i));
				assert(v >= get_params_min_for(m, i));
			}
		}
		if(i == get_n_par(m)) {
			calc_model(m, NULL);
			printf(DUMP_FORMAT "\t" DUMP_FORMAT "\n", get_prob(m), get_prior(m));
		}
	}
	return 0;
}
