#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>

#include <assert.h>
#include <string.h>
#include <ctype.h>
#include <signal.h>
#include "gsl_helper.h"
#include "debug.h"

int run = 1;

void ctrl_c_handler(int signalnr) {
	(void)signalnr;
	run = 0;
}

int main(void) {
	gsl_rng * rand;
	/* we want to prevent partial line output. */
	signal(SIGINT, ctrl_c_handler);

	gsl_rng_env_setup();
	rand = gsl_rng_alloc(gsl_rng_default);
	while (run) {
		printf("%e\n", gsl_rng_uniform(rand) );
	}
	return 0;
}
