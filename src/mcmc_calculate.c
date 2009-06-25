#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"

/*#define m->n_par 5*/

void mcmc_append_current_parameters(mcmc * m) {
	mcmc_dump_current(m);
	m->n_iter++;
}

void mcmc_check_best(mcmc * m) {
	if (m->prob > m->prob_best) {
		dump_v("found a better solution", m->params);
		m->prob_best = m->prob;
		set_params_best(m, m->params);
		mcmc_dump_model(m);
	}
}

