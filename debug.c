#include "debug.h"

void dump_vector(gsl_vector * v) {
	unsigned int i;
	printf("Vector%ud[", (unsigned int)v->size);
	for(i = 0; i < v->size - 1; i++) {
		printf("%f;", gsl_vector_get(v, i));
	}
	printf("%f]\n", gsl_vector_get(v, v->size - 1));
}

/**
 * Dump the mcmc structure
 */
void dump(mcmc * m) {
	unsigned int i;
	IFDEBUG {
		dump_p("dumping m", (void*)m);
		dump_i("n_par", get_n_par(m));
		dump_ul("accept", m->accept);
		dump_ul("reject", m->reject);
		dump_d("prob", m->prob);
		dump_d("prob_best", m->prob_best);
		for(i = 0; i < get_n_par(m); i++) {
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
			dump_size("x-size", m->x_dat->size);
		dump_p("y_dat",  (void*)m->y_dat);
		dump_p("model",  (void*)m->model);
		dump_size("size", m->size);
		dump_ul("iter", m->n_iter);
		debug("dumping m done ---- ");
	}
}

#ifndef require
void require(const int returncode) {
	if(returncode != 0) {
		fprintf(stderr, "a gsl call returned with an error: %s\n", gsl_strerror(returncode));
	}
}
#endif
