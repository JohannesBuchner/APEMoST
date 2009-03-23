/**
 * Public methods of the class
 */

#ifndef MCMC
#define MCMC

#include <stdio.h>
#include <stdlib.h>

#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif

#include "mcmc_struct.h"

/**
 * create and initialize a mcmc class using the configuration given in
 * @param filename
 */
mcmc * mcmc_load(const char * filename);
/**
 * frees the memory used by the class
 */
void mcmc_free(mcmc * m);

/**
 * adds the current parameter values to params_distr as nth iteration
 * was: add_values
 * @param m
 * @param n_iter
 */
void mcmc_append_current_parameters(mcmc * m, int n_iter);

void mcmc_dump_model(mcmc * m);
void mcmc_dump_y_dat(mcmc * m);

/**
 * write probability/distribution (params_distr) out to files.
 *
 * The files are named after the parameter names, with .prob.dump appended.
 * You can use those files to build a histogram of where the algorithm has been
 * (and how often) to find out
 *
 * @param m
 * @param n_values use the n last iterations. if negative, all are used.
 */
void mcmc_dump_probabilities(mcmc * m, int n_values);

/**
 * check if a new best value has been found
 */
void mcmc_check_best(mcmc * m);


/**
 * create/calibrate the markov-chain
 */
void markov_chain_calibrate(mcmc * m, double rat_limit, unsigned int burn_in_iterations,
		unsigned int iter_limit, double mul);
/**
 * take a step using the markov-chain
 */
void markov_chain_step(mcmc * m);




/* getter + setter */
const char ** get_params_descr(mcmc * m);
long get_params_accepts(mcmc * m);
long get_params_rejects(mcmc * m);
long get_params_accepts_for(mcmc * m, int i);
long get_params_rejects_for(mcmc * m, int i);
gsl_vector * get_params(mcmc * m);
gsl_vector * get_params_best(mcmc * m);
gsl_vector * get_x(mcmc * m);
gsl_vector * get_y(mcmc * m);
int get_n_par(mcmc * m);
gsl_rng * get_random(mcmc * m);
double get_prob(mcmc * m);
double get_prob_best(mcmc * m);

void set_prob(mcmc * m, double new_prob);
void set_prob_best(mcmc * m, double new_prob_best);
void set_minmax_for(mcmc * m, double new_min, double new_max, int i);
void set_model(mcmc * m, gsl_vector * new_model);
void set_n_par(mcmc * m, int new_n_par);
void set_params_best(mcmc * m, gsl_vector * new_params_best);
void set_params_for(mcmc * m, double new_param, int i);
void set_params(mcmc * m, gsl_vector * new_params);
void set_params_descr_all(mcmc * m, const char ** new_par_descr);
void set_params_descr_for(mcmc * m, const char * new_par_descr, int i);
void set_random(mcmc * m, gsl_rng * newrandom);
void set_prob(mcmc * m, double new_prob);
void set_x(mcmc * m, gsl_vector * new_x);
void set_x_copy(mcmc * m, gsl_vector * new_x);
void set_y(mcmc * m, gsl_vector * new_y);
void set_y_copy(mcmc * m, gsl_vector * new_y);
void set_steps_for(mcmc * m, double new_steps, int i);
void set_steps_all(mcmc * m, double * new_steps);

#endif

