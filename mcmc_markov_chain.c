#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"

void markov_chain_calibrate(mcmc * m, double rat_limit, unsigned int burn_in_iterations,
		unsigned int iter_limit, double mul) {
	(void)m;
	rat_limit+=burn_in_iterations+iter_limit + mul;
}
void markov_chain_step(mcmc * m) {
	(void)m;

}
