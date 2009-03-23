

#include <signal.h>
#include "mcmc.h"
#include "gsl_helper.h"
#include "simplesin.h"
#include "debug.h"
#include <gsl/gsl_sf.h>

#define DUMP_PROB_INTERVAL 2000
#define MAX_ITERATIONS 10000

double sigma;
int run;

void ctrl_c_handler(int signal);

void calc_model(mcmc * m) {
	unsigned int i;
	double x;
	double y;
	double param0 = gsl_vector_get(m->params, 0);
	double param1 = gsl_vector_get(m->params, 1);
	double param2 = gsl_vector_get(m->params, 2);

	for(i = 0; i < m->x_dat->size; i++) {
		x = gsl_vector_get(m->x_dat, i);
		y = param0 * gsl_sf_sin(2.0 * M_PI * param1 * x + param2);
		gsl_vector_set(m->model, i, y);
	}
}

void calc_prob(mcmc * m) {
	gsl_vector * diff = dup_vector(m->model);
	require(gsl_vector_sub(diff, m->y_dat));
	m->prob = calc_vector_squaresum(diff) / (2 * sigma*sigma);
}

void analyse(parallel_tempering_mcmc ** sinmod, int n_beta);

void simplesin(const char * filename) {
	int n_beta;
	int n_par;
	double beta_0;
	double delta_beta;
	int i;
	const char ** params_descr;
	parallel_tempering_mcmc ** sinmod;

	/* TODO: load these from config */
	sigma = 0.5;
	n_beta = 12;
	beta_0 = 0.001;

	delta_beta = (1.0 - beta_0) / (n_beta - 1);
	sinmod = (parallel_tempering_mcmc**)calloc(n_beta, sizeof(parallel_tempering_mcmc*));

	printf("Initializing parallel tempering for %d chains\n", n_beta);
	for(i = 0; i < n_beta; i++) {
		printf("\tChain %d - beta = %f\n", i+1, 1.0 - i * delta_beta);

		/* That is kind of stupid and could be optimized. not critical though. */
		sinmod[i]->m = mcmc_load(filename);

		set_beta(sinmod[i], 1.0 - i * delta_beta);
	}
	params_descr = get_params_descr(sinmod[0]->m);
	n_par = get_n_par(sinmod[0]->m);

	/* here was code for removing *_results.dat files.
	 * this should be done externally before calling the program. */

	debug("Initializing models");
	for(i = 0; i < n_beta; i++) {
		calc_model(sinmod[i]->m);
		calc_prob(sinmod[i]->m);
	}
	debug("Starting markov chain calibration");

	markov_chain_calibrate(sinmod[0]->m, 10000, 0.5, 20000, 0.85, DEFAULT_ADJUST_STEP);

	debug("Setting startingpoint for the calibration of all hotter distribution to the")
	debug("best parameter values of the (beta=1)-distribution")
	for(i = 0 + 1; i < n_beta; i++) {
		set_params(sinmod[i]->m, get_params_best(sinmod[i]->m));
		calc_model(sinmod[i]->m);
		calc_prob(sinmod[i]->m);
		sinmod[i]->m->prob *= sinmod[i]->beta;

		markov_chain_calibrate(sinmod[i]->m, 10000, 0.5, 20000, 0.85, DEFAULT_ADJUST_STEP);
	}
	signal(SIGINT, ctrl_c_handler);
	analyse(sinmod, n_beta);
}

#define randomu() gsl_rng_uniform(get_random(sinmod[0]->m))

void parallel_tempering_swap(parallel_tempering_mcmc ** sinmod, int n_beta, int n_swap) {
	double swap_probability;
	double a_beta, b_beta;
	double a_prob, b_prob;
	double r, c;
	int a, b;
	gsl_vector * temp;
	debug("checking if we do a swap");
	swap_probability = randomu();
	if (swap_probability < 1.0 / n_swap) {
		a = (n_beta - 1) * randomu();
		b = (a + 1) % n_beta;
		assert(a >= 0 && a < n_beta);
		assert(b >= 0 && b < n_beta);
		a_prob = get_prob(sinmod[a]->m);
		b_prob = get_prob(sinmod[b]->m);
		a_beta = get_beta(sinmod[a]);
		b_beta = get_beta(sinmod[a]);
		r = a_beta * b_prob / b_beta + b_beta * a_prob / a_beta - (a_prob
				+ b_prob);
		dump_i("swapping", a);
		dump_i("with", b);
		dump_d("likelihood: ", r);
		c = gsl_sf_log(randomu());
		if (r > c) {
			dump_d("really swapping", c);
			temp = get_params(sinmod[a]->m);
			set_params(sinmod[a]->m, get_params(sinmod[b]->m));
			set_params(sinmod[b]->m, temp);
			/* TODO: update best? */
		}else{
			dump_d("not swapping", c);
		}
	}
}

void analyse(parallel_tempering_mcmc ** sinmod, int n_beta) {
	int i;
	int n_swap = 30;
	unsigned long iter = sinmod[0]->m->n_iter;
	run = 1;
	debug("starting the analysis");

	while (run && iter < MAX_ITERATIONS) {
		/* TODO: maybe we can do 100 operations off these in threads using OpenMP ? */
		for (i = 0; i < n_beta; i++) {
			dump_i("one markov-chain step for ", i);
			markov_chain_step(sinmod[i]->m);
		}
		/* TODO: what happens if a different one (!=0) finds a really good solution? */
		mcmc_check_best(sinmod[0]->m);
		mcmc_append_current_parameters(sinmod[0]->m, iter);

		if (iter % DUMP_PROB_INTERVAL == DUMP_PROB_INTERVAL - 1) {
			debug("dumping distribution");
			dump_ul("iteration", iter);
			dump_ul("acceptance rate: accepts", get_params_accepts_sum(sinmod[0]->m));
			dump_ul("acceptance rate: rejects", get_params_rejects_sum(sinmod[0]->m));
			mcmc_dump_probabilities(sinmod[0]->m, DUMP_PROB_INTERVAL);
			dump(sinmod[0]->m);
		}

		parallel_tempering_swap(sinmod, n_beta, n_swap);
		iter++;
	}
}
#undef randumu

/* TODO: add ctrl-c handler */
void ctrl_c_handler(int signal) {
	printf("received Ctrl-C (%d). Stopping ... (please be patient)\n", signal);
	run = 0;
}

int main(void) {
	simplesin("simplesin/input");
	return 0;
}

