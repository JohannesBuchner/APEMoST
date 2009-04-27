#include "parallel_tempering_interaction.h"
#include "parallel_tempering.h"
#include "debug.h"
#include "mcmc_internal.h"
#include "gsl_helper.h"

int check_swap_probability(mcmc * a, mcmc * b) {
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
int parallel_tempering_decide_swap_random(mcmc ** sinmod, int n_beta,
		int n_swap) {
	double swap_probability;
	int a, b;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	IFVERBOSE
		debug("checking if we do a swap");
	swap_probability = get_next_uniform_random(sinmod[0]);
	if (swap_probability < 1.0 / n_swap) {
		a = (int) (n_beta * 1000 * get_next_uniform_random(sinmod[0])) % (n_beta - 1);
		b = (a + 1) % n_beta;
		if (check_swap_probability(sinmod[a], sinmod[b]) == 1)
			return a;
	}
	return -1;
}

/**
 * chooses one chain after another to swap. Swaps occur every n_swap.
 */
int parallel_tempering_decide_swap_nonrandom(mcmc ** sinmod, int n_beta,
		int n_swap, int iter) {
	int a;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	if (iter % n_swap == 0) {
		a = iter / n_swap % (n_beta - 1);
		if (check_swap_probability(sinmod[a], sinmod[a + 1]) == 1)
			return a;
	}
	return -1;
}

/**
 * chooses one chain to swap by random. Swap occurs every time.
 */
int parallel_tempering_decide_swap_now(mcmc ** sinmod, int n_beta) {
	int a;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	a = (int) (n_beta * 1000 * get_next_uniform_random(sinmod[0])) % (n_beta - 1);
	if (check_swap_probability(sinmod[a], sinmod[a + 1]) == 1)
		return a;
	return -1;
}

void parallel_tempering_do_swap(mcmc ** sinmod, int n_beta, int a) {
	double r;
	int b;
	gsl_vector * temp;
	assert(a < n_beta - 1);
	b = a + 1;
	IFDEBUG
		printf("swapping %d with %d\n", a, b);
	temp = dup_vector(get_params(sinmod[a]));
	set_params(sinmod[a], dup_vector(get_params(sinmod[b])));
	set_params(sinmod[b], temp);

	r = get_prob_best(sinmod[a]);
	if (r > get_prob_best(sinmod[b])) {
		set_prob_best(sinmod[b], r);
		set_params_best(sinmod[b], dup_vector(get_params_best(sinmod[a])));
	} else {
		r = get_prob_best(sinmod[b]);
		set_prob_best(sinmod[a], r);
		set_params_best(sinmod[a], dup_vector(get_params_best(sinmod[b])));
	}

	mcmc_check(sinmod[a]);
	mcmc_check(sinmod[b]);
}

void reset_to_best(mcmc ** sinmod, int n_beta) {

#ifdef RESET_TO_BEST
	int a;
	a = (int) (n_beta * RESET_TO_BEST * get_next_uniform_random(sinmod[0]));
	if (a < n_beta) {
		printf("Resetting chain %d to", a);
		dump_vector(get_params_best(sinmod[a]));
		set_prob(sinmod[a], get_prob_best(sinmod[a]));
		set_params(sinmod[a], dup_vector(get_params_best(sinmod[a])));
		calc_model(sinmod[a], NULL);
	}
#else
	(void)sinmod;
	(void)n_beta;
#endif
}

void tempering_interaction(mcmc ** sinmod, unsigned int n_beta,
		unsigned long iter) {
	int candidate;
	(void) iter;

#ifdef RANDOMSWAP
	candidate = parallel_tempering_decide_swap_random(sinmod, n_beta, 1);
#else
	candidate = parallel_tempering_decide_swap_now(sinmod, n_beta);
	/*candidate = parallel_tempering_decide_swap_nonrandom(sinmod, n_beta, n_swap, iter);*/
#endif
	if (candidate != -1) {
		/* wait for threads to reach iteration */
		parallel_tempering_do_swap(sinmod, n_beta, candidate);
		inc_swapcount(sinmod[candidate]);
	}
	reset_to_best(sinmod, n_beta);
}

