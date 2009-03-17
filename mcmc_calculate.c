#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "debug.h"

/*#define m->n_par 5*/

/*
 * TODO: write out to file, because memory will get too big.
 */
void mcmc_append_current_parameters(mcmc * m, int n_iter) {
	mcmc_prepare_iteration(m, n_iter);
	require(gsl_vector_memcpy(m->params_distr[n_iter], m->params));
}
