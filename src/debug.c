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
