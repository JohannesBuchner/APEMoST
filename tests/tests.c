
/* we do testing of the inner works */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mcmc.h"
#include "mcmc_internal.h"
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
int test_tests(void) {
	return (count_tests() > 0 ? 0 : 1);
}

int test_hist(void) {
	gsl_histogram * h;

	/*gsl_histogram * get_hist(gsl_vector * vs, int index, int nbins);*/
	gsl_vector * vs = gsl_vector_calloc(3);
	gsl_vector_set(vs, 0, 3);

	ASSERTEQUALI((int)vs->size, 3, "setup");
	h = calc_hist(vs, 3);
	ASSERTEQUALD(gsl_histogram_min(h), 0.0, "lower bound");
	ASSERTEQUALD(gsl_histogram_max(h), 3.0, "upper bound");
	ASSERTEQUALI((int)gsl_histogram_bins(h), 3, "nbins");

	ASSERTEQUALD(gsl_histogram_get (h, 0), 0.666667, "bin:0 -> 1 elem");
	ASSERTEQUALD(gsl_histogram_get (h, 1), 0.0, "bin:1 -> 0 elem");
	ASSERTEQUALD(gsl_histogram_get (h, 2), 0.333333, "bin:2 -> 2 elem");
	gsl_histogram_free(h);
	gsl_vector_free(vs);
	return 0;
}


int test_create(void) {
	mcmc * m;
	debug("test-create");
	m = mcmc_init(3);
	ASSERTEQUALI(m->n_par, 3, "number of parameters");
	dump(m);
	debug("freeing");
	m = mcmc_free(m);
	return 0;
}


int test_alloc(void) {
	unsigned int amount = 1000;
	unsigned int i;
	mcmc * m = mcmc_init(3);
	debug("allocating a lot\n");

	for(i = 0; i < amount; i++) {
		mcmc_prepare_iteration(m, i);
	}
	m = mcmc_free(m);
	return 0;
}

int test_load(void) {
	mcmc * m = mcmc_load("tests/testinput1", "tests/testlc.dat");
	ASSERT(m != NULL, "loaded");
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

	m = mcmc_free(m);
	return 0;
}


int test_resize(void) {
	mcmc * m = mcmc_init(3);
	int i = 0;
	int j = 0;
	for (j = 0; j < 5; j++) {
		for (; i % 64 != 63; i++) {
			mcmc_prepare_iteration(m, i);
		}
		i++;
		dump_i("tested iterations", i);
	}
	m = mcmc_free(m);
	return 0;
}

int test_append(void) {
	mcmc * m = mcmc_init(3);
	int i;
	for (i = 0; i < 10; i++) {
		mcmc_append_current_parameters(m);
		ASSERTEQUALI(i + 1, (int)m->n_iter, "number of iterations");
	}
	m = mcmc_free(m);
	return 0;
}


int test_random(void) {
	mcmc * m = mcmc_load("tests/testinput1", "tests/testlc.dat");
	double v = 0, last_v;
	int i;
	for (i = 0; i < 10; i++) {
		last_v = v;
		v = get_next_uniform_random(m);
		ASSERT( v != last_v, "v != last_v");
		ASSERT( v >= 0, "v >= 0");
		ASSERT( v <= 1, "v <= 1");
	}
	ASSERTEQUALD(gsl_sf_log(1E-100), -1E100, "ln 0");
	ASSERTEQUALD(gsl_sf_log(1.0), 0.0, "ln 1");
	ASSERT(get_next_alog_urandom(m) <= 0.0, "ln 1");
	mcmc_free(m);
	return 0;
}

int test_mod(void) {
	ASSERTEQUALD(mod_double(3.14, 3.00), 0.14, "mod ints");
	ASSERTEQUALD(mod_double(3.14, 1.30), 0.54, "mod doubles");
	ASSERTEQUALD(mod_double(-3.14, 1.30), 0.76, "mod doubles");
	ASSERTEQUALD(mod_double(0, 1.30), 0.00, "mod doubles");
	ASSERTEQUALD(mod_double(6000.3214, 1.1324), 0.8662, "mod doubles");
	ASSERTEQUALD(mod_double(-6000.3214, 1.1324), 0.2662, "mod doubles");
	return 0;
}

int test_write(void) {
	mcmc * m = mcmc_load("tests/testinput1", "tests/testlc.dat");
	debug("lets cheat and say we got the sum x-data + y-data as model");
	gsl_vector_memcpy(m->model, m->y_dat);
	gsl_vector_add(m->model, m->x_dat);
	mcmc_dump_model(m);
	m = mcmc_free(m);
	return 0;
}

int test_write_prob(void) {
	mcmc * m = mcmc_load("tests/testinput1", "tests/testlc.dat");
	debug("add starting points, ...");
	mcmc_append_current_parameters(m);
	mcmc_check(m);
	debug("maxima ...");
	require(gsl_vector_memcpy(m->params, m->params_max));
	mcmc_append_current_parameters(m);
	mcmc_check(m);
	debug("and minima as visited parameter values");
	require(gsl_vector_memcpy(m->params, m->params_max));
	mcmc_append_current_parameters(m);
	mcmc_check(m);
	mcmc_dump_probabilities(m, 1, "");
	ASSERTEQUALI(countlines("Amplitude.prob.dump"), 1, "" );
	ASSERTEQUALI(countlines("Frequenz.prob.dump"),  1, "" );
	mcmc_dump_probabilities(m, -1, "");
	ASSERTEQUALI(countlines("Amplitude.prob.dump"), 3, "" );
	ASSERTEQUALI(countlines("Frequenz.prob.dump"),  3, "" );
	m = mcmc_free(m);
	return 0;
}


void calc_prob(mcmc * m) {
	(void)m;
}
void calc_model(mcmc * m, const gsl_vector * old_values) {
	(void)m;
	(void)old_values;
}
void calc_model_for(mcmc * m, const unsigned int index, const double old_value) {
	(void)m;
	(void)index;
	(void)old_value;
}


/* register of all tests */
int (*tests_registration[])(void)  = {
	/* this is test 1 */ /*test_tests, */
	test_hist,
	test_create,
	test_alloc,
	test_load,
	test_resize,
	test_append,
	test_random,
	test_mod,
	test_write,
	test_write_prob,

	/* register more tests before here */
	NULL,
};


