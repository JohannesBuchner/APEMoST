#include <unistd.h>
#include <string.h>
#include "mcmc.h"
#include "debug.h"
#include "parallel_tempering.h"
#include "parallel_tempering_interaction.h"

#define PARAMS_FILENAME "params"
#define DATA_FILENAME "data"

int main(int argc, char ** argv) {
	unsigned int n;
	unsigned int p;
	unsigned int i;
	unsigned int j;
	mcmc * m;
	double prob;

	assert(argc == 2 + 1);
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
