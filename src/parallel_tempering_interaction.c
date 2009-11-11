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

#include "parallel_tempering_interaction.h"
#include "parallel_tempering.h"
#include "debug.h"
#include "mcmc_internal.h"
#include "gsl_helper.h"

static int check_swap_probability(mcmc * a, mcmc * b) {
	double a_beta, b_beta;
	double a_prob, b_prob;
	double r, c;
	/* this code is probably not threadsafe -- start */
	a_prob = get_prob(a);
	b_prob = get_prob(b);
	/* this code is probably not threadsafe -- end */
	a_beta = get_beta(a);
	b_beta = get_beta(b);
	r = a_beta * b_prob / b_beta + b_beta * a_prob / a_beta - (a_prob + b_prob);
	c = get_next_alog_urandom(a);
	if (r > c) {
		return 1;
	} else {
		return 0;
	}
}

/**
 * chooses a chain to swap by random. Swaps occur with probability 1/n_swap
 */
int parallel_tempering_decide_swap_random(mcmc ** chains, int n_beta,
		int n_swap) {
	double swap_probability;
	int a, b;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	IFVERBOSE
		debug("checking if we do a swap");
	swap_probability = get_next_uniform_random(chains[0]);
	if (swap_probability < 1.0 / n_swap) {
		a = (int) (n_beta * 1000 * get_next_uniform_random(chains[0]))
				% (n_beta - 1);
		b = (a + 1) % n_beta;
		if (check_swap_probability(chains[a], chains[b]) == 1)
			return a;
	}
	return -1;
}

/**
 * chooses one chain after another to swap. Swaps occur every n_swap.
 */
int parallel_tempering_decide_swap_nonrandom(mcmc ** chains, int n_beta,
		int n_swap, int iter) {
	int a;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	if (iter % n_swap == 0) {
		a = iter / n_swap % (n_beta - 1);
		if (check_swap_probability(chains[a], chains[a + 1]) == 1)
			return a;
	}
	return -1;
}

/**
 * chooses one chain to swap by random. Swap occurs every time.
 */
int parallel_tempering_decide_swap_now(mcmc ** chains, int n_beta) {
	int a;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	a = (int) (n_beta * 1000 * get_next_uniform_random(chains[0])) % (n_beta
			- 1);
	if (check_swap_probability(chains[a], chains[a + 1]) == 1)
		return a;
	return -1;
}

static void parallel_tempering_do_swap(mcmc ** chains, int n_beta, int a) {
	double r;
	int b;
	gsl_vector * temp;
	assert(a < n_beta - 1);
	b = a + 1;
	IFDEBUG
		printf("swapping %d with %d\n", a, b);
	temp = dup_vector(get_params(chains[a]));
	set_params(chains[a], dup_vector(get_params(chains[b])));
	set_params(chains[b], temp);

	r = get_prob_best(chains[a]);
	if (r > get_prob_best(chains[b])) {
		set_prob_best(chains[b], r);
		set_params_best(chains[b], get_params_best(chains[a]));
	} else {
		r = get_prob_best(chains[b]);
		set_prob_best(chains[a], r);
		set_params_best(chains[a], get_params_best(chains[b]));
	}

	mcmc_check(chains[a]);
	mcmc_check(chains[b]);
}

void tempering_interaction(mcmc ** chains, unsigned int n_beta,
		unsigned long iter) {
	int candidate;
	(void) iter;

#ifdef RANDOMSWAP
	candidate = parallel_tempering_decide_swap_random(chains, n_beta, 1);
#else
	candidate = parallel_tempering_decide_swap_now(chains, n_beta);
	/*candidate = parallel_tempering_decide_swap_nonrandom(chains, n_beta, n_swap, iter);*/
#endif
	if (candidate != -1) {
		/* wait for threads to reach iteration */
		parallel_tempering_do_swap(chains, n_beta, candidate);
		inc_swapcount(chains[candidate]);
	}
}

