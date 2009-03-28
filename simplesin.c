

#include <signal.h>
#include "mcmc.h"
#include "gsl_helper.h"
#include "simplesin.h"
#include "debug.h"
#include <gsl/gsl_sf.h>

#define DUMP_PROB_INTERVAL  1000
#define PRINT_PROB_INTERVAL 1000
#define MAX_ITERATIONS 1000000

double sigma;
int run;

void ctrl_c_handler(int signal);

double apply_formula(mcmc * m, unsigned int i, double param0, double param1, double param2) {
	double x = gsl_vector_get(m->x_dat, i);
	double y = param0 * gsl_sf_sin(2.0 * M_PI * param1 * x + param2);
	gsl_vector_set(m->model, i, y);
	return y;
}

void calc_model(mcmc * m, const gsl_vector * old_values) {
	unsigned int i;
	double param0 = gsl_vector_get(m->params, 0);
	double param1 = gsl_vector_get(m->params, 1);
	double param2 = gsl_vector_get(m->params, 2);
	double y;
	double square_sum = 0;

	if (old_values != NULL && calc_same(old_values, m->params) == 1)
		return;
	/*dump_v("recalculating model for parameter values", m->params);*/
	for(i = 0; i < m->x_dat->size; i++) {
		y = apply_formula(m, i, param0, param1, param2) - gsl_vector_get(m->y_dat, i);
		square_sum += y*y;
	}
	m->prob = square_sum / (-2 * sigma*sigma);
	/*debug("model done");*/
}
void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void)i;
	(void)old_value;

	calc_model(m, NULL);
}

/*
void calc_prob(mcmc * m) {
	gsl_vector * diff = dup_vector(m->model);
	require(gsl_vector_sub(diff, m->y_dat));
	m->prob = - calc_vector_squaresum(diff) / (2 * sigma*sigma);
	gsl_vector_free(diff);
}*/

void analyse(parallel_tempering_mcmc ** sinmod, int n_beta);

void simplesin(const char * filename) {
	int n_beta;
	int n_par;
	double beta_0;
	double delta_beta;
	int i;
	const char ** params_descr;
	long burn_in_iterations;
	double rat_limit;
	long iter_limit;
	double mul;
	parallel_tempering_mcmc ** sinmod;

	/* TODO: load these from config */
	sigma = 0.5;
	n_beta = 12 / 3;
	beta_0 = 0.001;
	burn_in_iterations = 10000;
	rat_limit = 0.5;
	iter_limit = 20000;
	mul = 0.85;

	delta_beta = (1.0 - beta_0) / (n_beta - 1);
	sinmod = (parallel_tempering_mcmc**)calloc(n_beta, sizeof(parallel_tempering_mcmc*));
	assert(sinmod != NULL);

	printf("Initializing parallel tempering for %d chains\n", n_beta);
	for(i = 0; i < n_beta; i++) {
		printf("\tChain %d - beta = %f\n", i, 1.0 - i * delta_beta);
		sinmod[i] = (parallel_tempering_mcmc*)malloc(sizeof(parallel_tempering_mcmc));
		assert(sinmod[i] != NULL);
		/* That is kind of stupid (doublicate execution) and could be optimized.
		 * not critical though. */
		sinmod[i]->m = mcmc_load(filename);
		assert(sinmod[i]->m->y_dat->size == sinmod[i]->m->x_dat->size);
		set_beta(sinmod[i], 1.0 - i * delta_beta);
		mcmc_check(sinmod[i]->m);
	}
	params_descr = get_params_descr(sinmod[0]->m);
	n_par = get_n_par(sinmod[0]->m);

	/* here was code for removing *_results.dat files.
	 * this should be done externally before calling the program. */

	printf("Initializing models\n");
	for(i = 0; i < n_beta; i++) {
		calc_model(sinmod[i]->m, NULL);
		mcmc_check(sinmod[i]->m);
	}
	printf("Starting markov chain calibration\n");
	wait();
	markov_chain_calibrate(sinmod[0]->m, burn_in_iterations, rat_limit, iter_limit, mul, DEFAULT_ADJUST_STEP);

	printf("Setting startingpoint for the calibration of all hotter distribution to \n");
	printf("  the best parameter values of the (beta=1)-distribution\n");
	wait();
	for(i = 0 + 1; i < n_beta; i++) {
		printf("\tCalibrating chain %d\n", i);
		set_params(sinmod[i]->m, dup_vector(get_params_best(sinmod[0]->m)));
		calc_model(sinmod[i]->m, NULL);
		sinmod[i]->m->prob *= get_beta(sinmod[i]);

		markov_chain_calibrate(sinmod[i]->m, burn_in_iterations, rat_limit, iter_limit, mul, DEFAULT_ADJUST_STEP);
	}
	printf("all chains calibrated.\n");
	wait();
	signal(SIGINT, ctrl_c_handler);
	analyse(sinmod, n_beta);
	for(i = 0; i < n_beta; i++) {
		sinmod[i]->m = mcmc_free(sinmod[i]->m);
		free(sinmod[i]);
	}
	free(sinmod);
}

void parallel_tempering_swap(parallel_tempering_mcmc ** sinmod, int n_beta, int n_swap) {
	double swap_probability;
	double a_beta, b_beta;
	double a_prob, b_prob;
	double r, c;
	int a, b;
	gsl_vector * temp;
	assert(n_beta > 0);
	if(n_beta == 1)
		return;
	IFVERBOSE
		debug("checking if we do a swap");
	swap_probability = get_next_urandom(sinmod[0]->m);
	if (swap_probability < 1.0 / n_swap) {
		a = (n_beta - 1) * get_next_urandom(sinmod[0]->m);
		b = (a + 1) % n_beta;
		assert(a >= 0 && a < n_beta);
		assert(b >= 0 && b < n_beta);
		mcmc_check(sinmod[a]->m);
		mcmc_check(sinmod[b]->m);
		a_prob = get_prob(sinmod[a]->m);
		b_prob = get_prob(sinmod[b]->m);
		a_beta = get_beta(sinmod[a]);
		b_beta = get_beta(sinmod[b]);
		r = a_beta * b_prob / b_beta + b_beta * a_prob / a_beta - (a_prob
				+ b_prob);
		printf("swapping %d with %d with probability %f\n", a, b, r);
		c = get_next_alog_urandom(sinmod[0]->m);
		if (r > c) {
			dump_d("we are really swapping", c);
			temp = dup_vector(get_params(sinmod[a]->m));
			set_params(sinmod[a]->m, dup_vector(get_params(sinmod[b]->m)));
			set_params(sinmod[b]->m, temp);
			mcmc_check(sinmod[a]->m);
			mcmc_check(sinmod[b]->m);
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
	printf("starting the analysis\n");
	wait();

	while (run && iter < MAX_ITERATIONS) {
		/* TODO: maybe we can do 100 operations off these in threads using OpenMP ? */
		for (i = 0; i < n_beta; i++) {
			/*dump_i("one markov-chain step for ", i);*/
			markov_chain_step(sinmod[i]->m, 0);
		}
		/* TODO: what happens if a different one (!=0) finds a really good solution? */
		mcmc_check_best(sinmod[0]->m);
		mcmc_append_current_parameters(sinmod[0]->m);

		if (iter % PRINT_PROB_INTERVAL == PRINT_PROB_INTERVAL - 1) {
			IFDEBUG; else {
				printf("iteration: %lu, a/r: %lu/%lu v:%.5e|%.5e|%.5e \r", iter,
						get_params_accepts_sum(sinmod[0]->m),
						get_params_rejects_sum(sinmod[0]->m),
						gsl_vector_get(get_params(sinmod[0]->m),0),
						gsl_vector_get(get_params(sinmod[0]->m),1),
						gsl_vector_get(get_params(sinmod[0]->m),2));
				fflush(stdout);
			}
			debug("dumping distribution");
			dump_ul("iteration", iter);
			dump_ul("acceptance rate: accepts", get_params_accepts_sum(sinmod[0]->m));
			dump_ul("acceptance rate: rejects", get_params_rejects_sum(sinmod[0]->m));
			dump(sinmod[0]->m);
		}
		if (iter % DUMP_PROB_INTERVAL == DUMP_PROB_INTERVAL - 1) {
			mcmc_dump_probabilities(sinmod[0]->m, -DUMP_PROB_INTERVAL*2);
		}

		parallel_tempering_swap(sinmod, n_beta, n_swap);
		iter++;
	}
}

void ctrl_c_handler(int signal) {
	printf("received Ctrl-C (%d). Stopping ... (please be patient)\n", signal);
	run = 0;
}

int main(void) {
	simplesin("simplesin/input");
	return 0;
}

