
/* we do testing of the inner works */
#include <stdio.h>
#include <stdlib.h>
#include "mcmc.h"
#include "gsl_helper.h"

void dump(void* a) {
	(void)a;
}

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
	h = calc_hist(&vs, 1, 3);
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
	mcmc m;
	mcmc_init(m);
	
	
	return 0;
}



/* register of all tests */
int (*tests_registration[])(void)  = {
	test_tests, /* this is test 1 */
	test_hist,
	
	/* register more tests before here */
	NULL,
};


