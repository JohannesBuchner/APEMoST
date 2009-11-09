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
