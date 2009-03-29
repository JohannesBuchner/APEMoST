

#include <signal.h>
#include <gsl/gsl_sf.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include "simplesin.h"
#include "debug.h"

#define DUMP_PROB_LENGTH  -1
#define PRINT_PROB_INTERVAL 1000
#ifdef BENCHMARK
#define MAX_ITERATIONS   40000
#else
#define MAX_ITERATIONS 1000000
#endif

double sigma;
int run;
int dumpflag;

void ctrl_c_handler(int signal);
void sigusr_handler(int signal);

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
	for (i = 0; i < m->x_dat->size; i++) {
		y = apply_formula(m, i, param0, param1, param2) - gsl_vector_get(
				m->y_dat, i);
		square_sum += y * y;
	}
	set_prob(m, get_beta(m)* square_sum / (-2 * sigma*sigma));
	/*debug("model done");*/
}

void calc_model_for(mcmc * m, const unsigned int i, const double old_value) {
	(void)i;
	(void)old_value;

	calc_model(m, NULL);
}

void analyse(mcmc ** sinmod, int n_beta);

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
	mcmc ** sinmod;

	/* TODO: load these from config */
	sigma = 0.5;
	n_beta = 12 / 3;
	beta_0 = 0.001;
	burn_in_iterations = 10000;
	rat_limit = 0.5;
	iter_limit = 20000;
	mul = 0.85;

	delta_beta = (1.0 - beta_0) / (n_beta - 1);
	sinmod = (mcmc**)calloc(n_beta, sizeof(mcmc*));
	assert(sinmod != NULL);

	printf("Initializing parallel tempering for %d chains\n", n_beta);
	for(i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, 1.0 - i * delta_beta);
		/* That is kind of stupid (duplicate execution) and could be optimized.
		 * not critical though. */
		sinmod[i] = mcmc_load(filename);
		mcmc_check(sinmod[i]);
		sinmod[i]->additional_data = malloc(sizeof(parallel_tempering_mcmc));
		set_beta(sinmod[i], 1.0 - i * delta_beta);
		/*gsl_vector_scale(sinmod[i]->params_step, 2000);*/
		printf("\tsteps: "); dump_vector(get_steps(sinmod[i]));
		mcmc_check(sinmod[i]);
	}
	params_descr = get_params_descr(sinmod[0]);
	n_par = get_n_par(sinmod[0]);

	/* here was code for removing *_results.dat files.
	 * this should be done externally before calling the program. */

	printf("Initializing models\n");
	for(i = 0; i < n_beta; i++) {
		calc_model(sinmod[i], NULL);
		mcmc_check(sinmod[i]);
	}
	printf("Starting markov chain calibration\n");
	wait();
	markov_chain_calibrate(sinmod[0], burn_in_iterations, rat_limit, iter_limit, mul, DEFAULT_ADJUST_STEP);

	printf("Setting startingpoint for the calibration of all hotter distribution to \n");
	printf("  the best parameter values of the (beta=1)-distribution\n");
	wait();
	for(i = 0 + 1; i < n_beta; i++) {
		printf("\tCalibrating chain %d\n", i);
		set_params(sinmod[i], dup_vector(get_params_best(sinmod[0])));
		calc_model(sinmod[i], NULL);
		sinmod[i]->prob *= get_beta(sinmod[i]);

		markov_chain_calibrate(sinmod[i], burn_in_iterations, rat_limit, iter_limit, mul, DEFAULT_ADJUST_STEP);
	}
	printf("all chains calibrated.\n");
	for(i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, 1.0 - i * delta_beta);
		printf("\tsteps: "); dump_vector(get_steps(sinmod[i]));
	}
	wait();

	signal(SIGINT, ctrl_c_handler);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);

	analyse(sinmod, n_beta);
	for(i = 0; i < n_beta; i++) {
		sinmod[i] = mcmc_free(sinmod[i]);
		free(sinmod[i]);
	}
	free(sinmod);
}

void parallel_tempering_swap(mcmc ** sinmod, int n_beta, int n_swap) {
	double swap_probability;
	double a_beta, b_beta;
	double a_prob, b_prob;
	double r, c;
	int a, b;
	gsl_vector * temp;
	assert(n_beta > 0);
	if(n_beta == 1)
		return;
	(void) n_swap;
	IFVERBOSE
		debug("checking if we do a swap");
	swap_probability = get_next_urandom(sinmod[0]);
	if (swap_probability < 1.0 / 10000) {
		/* reset chain */
		a = (int)(n_beta*1000 * get_next_urandom(sinmod[0])) % n_beta;
		set_params(sinmod[a], dup_vector(get_params_best(sinmod[a])));
	}else if (swap_probability < 1.0 / n_swap) {
		a = (int)(n_beta*1000 * get_next_urandom(sinmod[0])) % n_beta;
		b = (a + 1) % n_beta;
		assert(a >= 0 && a < n_beta);
		assert(b >= 0 && b < n_beta);
		mcmc_check(sinmod[a]);
		mcmc_check(sinmod[b]);
		a_prob = get_prob(sinmod[a]);
		b_prob = get_prob(sinmod[b]);
		a_beta = get_beta(sinmod[a]);
		b_beta = get_beta(sinmod[b]);
		r = a_beta * b_prob / b_beta + b_beta * a_prob / a_beta - (a_prob
				+ b_prob);
		c = get_next_alog_urandom(sinmod[0]);
		if (r > c) {
			IFDEBUG printf("swapping %d with %d with probability %f\n", a, b, r);
			dump_d("we are really swapping", c);
			temp = dup_vector(get_params(sinmod[a]));
			set_params(sinmod[a], dup_vector(get_params(sinmod[b])));
			set_params(sinmod[b], temp);

			temp = dup_vector(get_params_best(sinmod[a]));
			set_params_best(sinmod[a], dup_vector(get_params_best(sinmod[b])));
			set_params_best(sinmod[b], temp);

			mcmc_check(sinmod[a]);
			mcmc_check(sinmod[b]);
			/* TODO: update best? */


		}else{
			IFVERBOSE dump_d("not swapping", c);
		}
	}
}

#include <time.h>

int get_duration() {
	static clock_t stored = 0;
	clock_t new = stored;
	stored = clock();
	return stored - new;
}

void print_current_positions(mcmc ** sinmod, int n_beta){
	int i;
	printf("printing chain parameters: \n");
	for(i = 0; i < n_beta; i++) {
		printf("\tchain %d: current: ", i); dump_vector(get_params(sinmod[i]));
		printf("\tchain %d: best: ", i);    dump_vector(get_params_best(sinmod[i]));
	}
	fflush(stdout);
}

void analyse(mcmc ** sinmod, int n_beta) {
	int i;
	int n_swap = 30;
	unsigned long iter = sinmod[0]->n_iter;
	/*int subiter;*/
	get_duration();
	run = 1;
	dumpflag = 0;
	printf("starting the analysis\n");
	wait();

	while (run && iter < MAX_ITERATIONS) {
		/* TODO: maybe we can do 100 operations off these in threads using OpenMP ? */
		for (i = 0; i < n_beta; i++) {
			/*dump_i("one markov-chain step for ", i);*/
			markov_chain_step(sinmod[i], 0);
		}
		mcmc_check_best(sinmod[0]);
		mcmc_append_current_parameters(sinmod[0]);
		/* TODO: what happens if a different one (!=0) finds a really good solution? */

		iter++;
		parallel_tempering_swap(sinmod, n_beta, n_swap);

		if (iter % PRINT_PROB_INTERVAL == 0) {
			if (dumpflag) {
				mcmc_dump_probabilities(sinmod[0], DUMP_PROB_LENGTH);
				print_current_positions(sinmod, n_beta);
				dumpflag = 0;
			}
			IFDEBUG {
				debug("dumping distribution");
				dump_ul("iteration", iter);
				dump_ul("acceptance rate: accepts", get_params_accepts_sum(sinmod[0]));
				dump_ul("acceptance rate: rejects", get_params_rejects_sum(sinmod[0]));
				dump(sinmod[0]);
			} else {
				printf(
						"iteration: %lu, a/r: %lu/%lu v:%.5e|%.5e|%.5e [%d/%lu ticks]\r",
						iter, get_params_accepts_sum(sinmod[0]),
						get_params_rejects_sum(sinmod[0]), gsl_vector_get(
								get_params(sinmod[0]), 0), gsl_vector_get(
								get_params(sinmod[0]), 1), gsl_vector_get(
								get_params(sinmod[0]), 2), get_duration(),
						CLOCKS_PER_SEC);
				fflush(stdout);
			}
		}
	}
}

void ctrl_c_handler(int signalnr) {
	printf("\nreceived Ctrl-C (%d). Stopping ... (please be patient)\n\n", signalnr);
	run = 0;
}
void sigusr_handler(int signalnr) {
	printf("\nreceived SIGUSR (%d). Will dump at next opportunity.\n\n", signalnr);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);
	dumpflag = 1;
}

int main(void) {
	simplesin("simplesin/input");
	return 0;
}

