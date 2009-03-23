#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"

/*#define m->n_par 5*/

/*
 * TODO: write out to file, because memory will get too big.
 */
void mcmc_append_current_parameters(mcmc * m, int n_iter) {
	mcmc_prepare_iteration(m, n_iter);
	require(gsl_vector_memcpy(m->params_distr[n_iter], m->params));
	m->n_iter++;
}

void mcmc_check_best(mcmc * m) {

	if(m->prob_best < m->prob) {
		m->prob_best = m->prob;
		gsl_vector_free(m->params_best);
		m->params_best = dup_vector(m->params);
		mcmc_dump_model(m);
	}

}
