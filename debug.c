#include "debug.h"

/**
 * Dump the mcmc structure
 */
void dump(mcmc * m) {
	unsigned int i;
	dump_p("dumping m", (void*)m);
	dump_i("n_par", m->n_par);
	dump_ul("accept", m->accept);
	dump_ul("reject", m->reject);
	dump_d("prob", m->prob);
	dump_d("prob_best", m->prob_best);
	for(i = 0; i < m->n_par; i++) {
		dump_i("parameter", i);
		if(m->params_descr != NULL)
			dump_i_s("\tparameter name", i, m->params_descr[i]);
		dump_d("\tvalue", gsl_vector_get(m->params, i));
		dump_d("\tbest", gsl_vector_get(m->params_best, i));
		dump_ul("\taccepts", m->params_accepts[i]);
		dump_ul("\trejects", m->params_rejects[i]);
		dump_d("\tstep-size", gsl_vector_get(m->params_step, i));
		dump_d("\tmin", gsl_vector_get(m->params_min, i));
		dump_d("\tmax", gsl_vector_get(m->params_max, i));
	}
	dump_p("x_dat",  (void*) m->x_dat);
	if(m->x_dat != NULL)
		dump_ul("x-size", m->x_dat->size);
	dump_p("y_dat",  (void*)m->y_dat);
	dump_p("model",  (void*)m->model);
	dump_ul("size", m->size);
	dump_ul("iter", m->n_iter);
	debug("dumping m done ---- ");
}

void require(const int returncode) {
	if(returncode != 0) {
		fprintf(stderr, "a gsl call returned with an error: %s\n", gsl_strerror(returncode));
	}
}
