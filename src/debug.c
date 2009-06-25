#include "debug.h"

void dump_vector(const gsl_vector * v) {
	unsigned int i;
	printf("Vector%ud[", (unsigned int) v->size);
	for (i = 0; i < v->size - 1; i++) {
		printf("%f;", gsl_vector_get(v, i));
	}
	printf("%f]", gsl_vector_get(v, v->size - 1));
}

void dump_vectorln(const gsl_vector * v) {
	dump_vector(v);
	printf("\n");
}

/**
 * Dump the mcmc structure
 */
void dump_mcmc(const mcmc * m) {
	unsigned int i;
	IFDEBUG {
		debug("dumping m      ---- ");
		IFSEGV
			dump_p("dumping m@%p\n", (void*)m);
		printf("\t\tn_par=%d; a/r=%lu/%lu prob/best=%f/%f\n", get_n_par(m),
				m->accept, m->reject, m->prob, m->prob_best);
		for (i = 0; i < get_n_par(m); i++) {
			if (m->params_descr != NULL)
				dump_i_s("parameter", i, m->params_descr[i]);
			dump_ul("\taccepts", m->params_accepts[i]);
			dump_ul("\trejects", m->params_rejects[i]);
		}
		dump_v("values", m->params);
		dump_v("best", m->params_best);
		dump_v("min", m->params_min);
		dump_v("max", m->params_max);
		dump_v("step-size", m->params_step);

		IFSEGV
			dump_p("x_dat", (void*) m->x_dat);
		if (m->x_dat != NULL)
			dump_size("x-size", m->x_dat->size);
		IFSEGV
			dump_p("y_dat", (void*)m->y_dat);
		IFSEGV
			dump_p("model", (void*)m->model);
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
