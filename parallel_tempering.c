#include <signal.h>
#include <gsl/gsl_sf.h>
#include <omp.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include "parallel_tempering.h"
#include "debug.h"

#ifdef BENCHMARK
#define MAX_ITERATIONS  100000
#else
#define MAX_ITERATIONS 1000000
#endif

int run;
int dumpflag;

void set_beta(mcmc * m, double newbeta) {
	((parallel_tempering_mcmc *) m->additional_data)->beta = newbeta;
}
double get_beta(mcmc * m) {
	return ((parallel_tempering_mcmc *) m->additional_data)->beta;
}

void ctrl_c_handler(int signalnr) {
	printf("\nreceived Ctrl-C (%d). Stopping ... (please be patient)\n\n",
			signalnr);
	run = 0;
}
void sigusr_handler(int signalnr) {
	printf("\nreceived SIGUSR (%d). Will dump at next opportunity.\n\n",
			signalnr);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);
	dumpflag = 1;
}

#include <time.h>

int get_duration() {
	static clock_t stored = 0;
	clock_t new = stored;
	stored = clock();
	return stored - new;
}

void print_current_positions(mcmc ** sinmod, int n_beta) {
	int i;
	printf("printing chain parameters: \n");
	for (i = 0; i < n_beta; i++) {
		printf("\tchain %d: current: ", i);
		dump_vectorln(get_params(sinmod[i]));
		printf("\tchain %d: best: ", i);
		dump_vectorln(get_params_best(sinmod[i]));
	}
	fflush(stdout);
}

void report(mcmc ** sinmod, int n_beta, int swapcount) {
	print_current_positions(sinmod, n_beta);
	printf("swapped %d times.", swapcount);
	mcmc_dump_probabilities(sinmod[0], -1);
}

int run;
int dumpflag;

void analyse(mcmc ** sinmod, int n_beta);

void parallel_tempering(const char * filename, int n_beta, double beta_0,
		unsigned long burn_in_iterations, double rat_limit,
		unsigned long iter_limit, double mul) {
	int n_par;
	double delta_beta;
	int i;
	const char ** params_descr;
	mcmc ** sinmod;

	/* TODO: load these from config */
	delta_beta = (1.0 - beta_0) / (n_beta - 1);
	sinmod = (mcmc**) calloc(n_beta, sizeof(mcmc*));
	assert(sinmod != NULL);

	printf("Initializing parallel tempering for %d chains\n", n_beta);
	#pragma omp parallel for
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, 1.0 - i * delta_beta);
		/* That is kind of stupid (duplicate execution) and could be optimized.
		 * not critical though. */
		sinmod[i] = mcmc_load(filename);
		mcmc_check(sinmod[i]);
		sinmod[i]->additional_data = malloc(sizeof(parallel_tempering_mcmc));
		set_beta(sinmod[i], 1.0 - i * delta_beta);
		printf("\tsteps: ");
		dump_vectorln(get_steps(sinmod[i]));
		mcmc_check(sinmod[i]);
	}
	params_descr = get_params_descr(sinmod[0]);
	n_par = get_n_par(sinmod[0]);

	/* here was code for removing *_results.dat files.
	 * this should be done externally before calling the program. */

	printf("Initializing models\n");
	for (i = 0; i < n_beta; i++) {
		calc_model(sinmod[i], NULL);
		mcmc_check(sinmod[i]);
	}
	printf("Starting markov chain calibration\n");
	wait();
	markov_chain_calibrate(sinmod[0], burn_in_iterations, rat_limit,
			iter_limit, mul, DEFAULT_ADJUST_STEP);

	printf(
			"Setting startingpoint for the calibration of all hotter distribution to \n");
	printf("  the best parameter values of the (beta=1)-distribution\n");
	wait();

	/* this could be parallelized */
	#pragma omp parallel for
	for (i = 0 + 1; i < n_beta; i++) {
		printf("\tCalibrating chain %d\n", i);
		set_params(sinmod[i], dup_vector(get_params_best(sinmod[0])));
		calc_model(sinmod[i], NULL);

		markov_chain_calibrate(sinmod[i], burn_in_iterations, rat_limit,
				iter_limit, mul, DEFAULT_ADJUST_STEP);
	}
	printf("all chains calibrated.\n");
	for (i = 0; i < n_beta; i++) {
		printf("\tChain %2d - beta = %f ", i, 1.0 - i * delta_beta);
		printf("\tsteps: ");
		dump_vectorln(get_steps(sinmod[i]));
	}
	wait();

	signal(SIGINT, ctrl_c_handler);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);

	analyse(sinmod, n_beta);
	for (i = 0; i < n_beta; i++) {
		sinmod[i] = mcmc_free(sinmod[i]);
		free(sinmod[i]);
	}
	free(sinmod);
}

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
	r = a_beta * b_prob / b_beta + b_beta * a_prob / a_beta - (a_prob
			+ b_prob);
	c = get_next_alog_urandom(a);
	if (r > c) {
		return 1;
	} else {
		return 0;
	}
}

int parallel_tempering_decide_swap_random(mcmc ** sinmod, int n_beta, int n_swap) {
	double swap_probability;
	int a, b;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	IFVERBOSE
		debug("checking if we do a swap");
	swap_probability = get_next_urandom(sinmod[0]);
	if (swap_probability < 1.0 / n_swap) {
		a = (int) (n_beta * 1000 * get_next_urandom(sinmod[0])) % n_beta;
		b = (a + 1) % n_beta;
		if(check_swap_probability(sinmod[a], sinmod[b]) == 1)
			return a;
	}
	return -1;
}
int parallel_tempering_decide_swap_nonrandom(mcmc ** sinmod, int n_beta, int n_swap, int iter) {
	int a, b;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	if (iter % n_swap == 0) { 
		a = iter/n_swap % n_beta;
		b = (a + 1) % n_beta;
		if(check_swap_probability(sinmod[a], sinmod[b]) == 1)
			return a;
	}
	return -1;
}
int parallel_tempering_decide_swap_now(mcmc ** sinmod, int n_beta, int n_swap, int iter) {
	int a, b;
	assert(n_beta > 0);
	if (n_beta == 1)
		return -1;
	a = iter/n_swap % n_beta;
	b = (a + 1) % n_beta;
	if(check_swap_probability(sinmod[a], sinmod[b]) == 1)
		return a;
	return -1;
}

void parallel_tempering_do_swap(mcmc ** sinmod, int n_beta, int a) {
	double r;
	int b;
	gsl_vector * temp;
	b = (a + 1) % n_beta;
	IFDEBUG
		printf("swapping %d with %d\n", a, b);
	temp = dup_vector(get_params(sinmod[a]));
	set_params(sinmod[a], dup_vector(get_params(sinmod[b])));
	set_params(sinmod[b], temp);

	r = get_prob_best(sinmod[a]);
	if( r > get_prob_best(sinmod[b])) {
		set_prob_best(sinmod[b], r);
		set_params_best(sinmod[b], dup_vector(get_params_best(sinmod[a])));
	}else{
		r = get_prob_best(sinmod[b]);
		set_prob_best(sinmod[a], r);
		set_params_best(sinmod[a], dup_vector(get_params_best(sinmod[b])));
	}

	mcmc_check(sinmod[a]);
	mcmc_check(sinmod[b]);
}

void analyse(mcmc ** sinmod, int n_beta) {
	int i;
	int n_swap = 30;
	unsigned long iter = sinmod[0]->n_iter;
	int candidate;
	int subiter;
	int swapcount = 0;
	get_duration();
	run = 1;
	dumpflag = 0;
	printf("starting the analysis\n");
	wait();

	while (run && iter < MAX_ITERATIONS) {
		/* TODO: maybe we can do 100 operations off these in threads using OpenMP ? */
		#pragma omp parallel for
		for (i = 0; i < n_beta; i++) {
			for (subiter = 0; subiter < n_swap; subiter++) {
				markov_chain_step(sinmod[i], 0);
				mcmc_check_best(sinmod[i]);
				if (i == 0) {
					mcmc_append_current_parameters(sinmod[0]);
				}
			}
		}
		/* do this in master-thread */
		iter += n_swap;
		#ifdef RANDOMSWAP		
		candidate = parallel_tempering_decide_swap_random(sinmod, n_beta, 1);
		#else
		candidate = parallel_tempering_decide_swap_now(sinmod, n_beta, n_swap, iter);
		/*candidate = parallel_tempering_decide_swap_nonrandom(sinmod, n_beta, n_swap, iter);*/
		#endif
		if(candidate != -1) {
			/* wait for threads to reach iteration */
			parallel_tempering_do_swap(sinmod, n_beta, candidate);
			swapcount++;
		}
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
				printf("iteration: %lu, a/r: %lu/%lu v:", iter,
						get_params_accepts_sum(sinmod[0]),
						get_params_rejects_sum(sinmod[0]));
				dump_vector(get_params(sinmod[0]));
				printf(" [%d/%lu ticks]\r", get_duration(), CLOCKS_PER_SEC);
				fflush(stdout);
			}
		}
	}
	report(sinmod, n_beta, swapcount);
}

