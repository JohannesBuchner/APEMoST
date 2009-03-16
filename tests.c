
/* we do testing of the inner works */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcmc.h"
#include "debug.h"
#include "gsl_helper.h"

#define DUMPONFAIL 1

#define ASSERTDUMP(condition, text, a) { \
		if(!(condition)) { \
			printf("  ASSERT FAILED: %s\n",text); \
			if(DUMPONFAIL && a != NULL){ \
				printf("    dump: \n"); \
				dump(a); \
			} \
			return 1; \
		} else { \
			printf("  subtest ok: %s\n",text); \
		} \
	}

#define ASSERT(condition, text) ASSERTDUMP(condition, text, NULL)
#define ASSERTEQUALI(result, expected, text) { \
		if(expected != result) { \
			printf("  ASSERT EQUAL FAILED: expected: %d; got: %d; %s\n", expected, result, text); \
			return 1; \
		} else { \
			printf("  subtest ok: %s\n", text); \
		} \
	}
#define ASSERTEQUALD(result, expected, text) { \
		if(expected == result || expected - result < (expected+result)*0.001 ) { \
			printf("  subtest ok: %s\n", text); \
		} else { \
			printf("  ASSERT EQUAL FAILED: expected: %e; got: %e; %s\n", expected, result, text); \
			return 1; \
		} \
	}
int count_tests();

/* example test function. return value is 0 iff succeeded.*/
int test_tests(void){
	return (count_tests()>0 ? 0 : 1);
}

int test_hist(void){
	gsl_histogram * h;
	
	/*gsl_histogram * get_hist(gsl_vector * vs, int index, int nbins);*/
	gsl_vector * vs = gsl_vector_calloc (3);
	gsl_vector_set (vs, 0, 3);
	
	ASSERTEQUALI((int)vs->size, 3, "setup");
	h = calc_hist(&vs, 0, 3);
	ASSERTEQUALD(gsl_histogram_min(h), 0.0, "lower bound");
	ASSERTEQUALD(gsl_histogram_max(h), 3.0, "upper bound");
	ASSERTEQUALI((int)gsl_histogram_bins(h), 3, "nbins");

	ASSERTEQUALD(gsl_histogram_get (h, 0), 0.666667, "bin:0 -> 1 elem");
	ASSERTEQUALD(gsl_histogram_get (h, 1), 0.0, "bin:1 -> 0 elem");
	ASSERTEQUALD(gsl_histogram_get (h, 2), 0.333333, "bin:2 -> 2 elem");
	gsl_histogram_free (h);

	return 0;
}


int test_create(void){
	mcmc * m;
	debug("test-create");
	m = mcmc_init(3);
	ASSERTEQUALI(m->n_par, 3, "number of parameters");
	dump(m);
	debug("freeing");
	mcmc_free(m);
	return 0;
}


int test_load(void){
	mcmc * m = mcmc_load("tests/testinput1");
	ASSERT(m!=NULL, "loaded");
	ASSERTEQUALI(m->n_par, 3, "number of parameters");
	ASSERTEQUALD(gsl_vector_get(m->params, 0),      0.7,  "start");
	ASSERTEQUALD(gsl_vector_get(m->params_min, 0),  0.4,  "min");
	ASSERTEQUALD(gsl_vector_get(m->params_max, 0),  3.0,  "max");
	ASSERTEQUALD(gsl_vector_get(m->params_step, 0), 0.01, "step");
	ASSERT(strcmp(m->params_descr[0], "Amplitude")==0, "description");
	ASSERTEQUALD(gsl_vector_get(m->params, 1),      5.2,  "start");
	ASSERTEQUALD(gsl_vector_get(m->params_min, 1),  4.0,  "min");
	ASSERTEQUALD(gsl_vector_get(m->params_max, 1),  24.0, "max");
	ASSERTEQUALD(gsl_vector_get(m->params_step, 1), 0.01, "step");
	ASSERT(strcmp(m->params_descr[1], "Frequenz")==0, "description");
	ASSERTEQUALD(gsl_vector_get(m->params, 2),      5.4,  "start");
	ASSERTEQUALD(gsl_vector_get(m->params_min, 2),  0.0,  "min");
	ASSERTEQUALD(gsl_vector_get(m->params_max, 2),  6.1,  "max");
	ASSERTEQUALD(gsl_vector_get(m->params_step, 2), 0.01, "step");
	ASSERT(strcmp(m->params_descr[2], "Phase")==0, "description");

	ASSERT(m->x_dat->size == 1522, "x_dat size");
	ASSERT(m->y_dat->size == m->x_dat->size, "y_dat size");
	ASSERT(m->model->size == m->x_dat->size, "model size");
	ASSERTEQUALD(gsl_vector_get(m->x_dat, 0),       1.7355217099999998,  "x 0");
	ASSERTEQUALD(gsl_vector_get(m->x_dat, 1),       1.7356002600000000,  "x 1");
	ASSERTEQUALD(gsl_vector_get(m->x_dat, 1521),       46.8043750000000003,  "x last");
	ASSERTEQUALD(gsl_vector_get(m->y_dat, 0),       0.4731745866773314,  "y 0");
	ASSERTEQUALD(gsl_vector_get(m->y_dat, 1),      -0.9900871130450772,  "y 1");
	ASSERTEQUALD(gsl_vector_get(m->y_dat, 1521),   -0.3527955490681067,  "y last");

	mcmc_free(m);
	return 0;
}


int test_resize(void){
	mcmc * m = mcmc_init(3);
	int i = 0;
	int j = 0;
	for(j = 0; j < 5; j++) {
		for(; i % 1024 != 1023; i++){
			prepare_iter(m, i);
		}
		i++;
		dump_i("tested iterations", i);
	}
	mcmc_free(m);
	return 0;
}


int test_write(void){
	mcmc * m = mcmc_load("tests/testinput1");
	
	mcmc_free(m);
	return 0;
}




/* register of all tests */
int (*tests_registration[])(void)  = {
	/* this is test 1 */ /*test_tests, */
	test_hist,
	test_create,
	test_load,
	test_resize,
	test_write,
	
	/* register more tests before here */
	NULL,
};


