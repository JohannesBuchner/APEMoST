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
	unsigned int n;
	unsigned int p;
	unsigned int i;
	unsigned int j;
	mcmc * m;
	double prob;

	if(argc != 2 + 1) {
		fprintf(stderr, "SYNOPSIS: %s <npartialcalc> <ncalc>\n"
				"\n"
				"This program calculates the model the given number of times\n"
				"to allow benchmarking of the model calculations.\n"
				"\n"
				"\tnpartialcalc\tnumber of calls to calc_model_for\n"
				"\tncalc\tnumber of calls to calc_model\n"
				"\n"
				"APEMoST  Copyright (C) 2009  Johannes Buchner\n"
				"This program comes with ABSOLUTELY NO WARRANTY; for details see the file LICENSE.\n"
				"This is free software, and you are welcome to redistribute it\n"
				"under certain conditions; see the file LICENSE.\n"
				"", argv[0]);
		exit(1);
	}
	assert(atoi(argv[1]) >= 0);
	assert(atoi(argv[2]) >= 0);
	n = atoi(argv[1]);
	p = atoi(argv[2]);

	m = mcmc_load_params(PARAMS_FILENAME);
	mcmc_load_data(m, DATA_FILENAME);
	mcmc_check(m);
	m->additional_data = mem_malloc(sizeof(parallel_tempering_mcmc));
	set_beta(m, 1.0);
	calc_model(m, NULL);
	prob = get_prob(m);
	for (i = 0; i < n;) {
		for (j = 0; j < get_n_par(m); j++) {
			calc_model_for(m, j, gsl_vector_get(get_params(m), j));
			i++;
			if ((prob - get_prob(m)) / (prob + get_prob(m)) > 1e-7) {
				dump_d("original prob", prob);
				dump_d("new prob", get_prob(m));
				assert(prob == get_prob(m));
			}
		}
	}
	for (; i < n + p; i++) {
		calc_model(m, NULL);
		assert(prob == get_prob(m));
	}
	return 0;
}
