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

/* getter + setter */
long get_params_accepts(mcmc * m);
long get_params_rejects(mcmc * m);
long get_params_accepts_for(mcmc * m, int i);
long get_params_rejects_for(mcmc * m, int i);
double get_params_ar_for(mcmc * m, int i);
void set_params_ar(mcmc * m, gsl_vector ** new_params_ar);
void set_params_ar_for(mcmc * m, gsl_vector * new_params_ar, int i);
void set_prob_best(mcmc * m, double new_prob_best);
void set_minmax_for(mcmc * m, double new_min, double new_max, int i);
void set_model(mcmc * m, gsl_vector * new_model);
void set_n_par(mcmc * m, int new_n_par);
void set_params_best(mcmc * m, gsl_vector * new_params_best);
void set_params_for(mcmc * m, double new_param, int i);
void set_params(mcmc * m, gsl_vector * new_params);
void set_params_descr_all(mcmc * m, char ** new_par_descr);
void set_params_descr_for(mcmc * m, char * new_par_descr, int i);
void set_seed(mcmc * m, gsl_vector * new_seed);
void set_probability(mcmc * m, double new_prob);
void set_x(mcmc * m, gsl_vector * new_x);
void set_x_copy(mcmc * m, gsl_vector * new_x);
void set_y(mcmc * m, gsl_vector * new_y);
void set_y_copy(mcmc * m, gsl_vector * new_y);
void set_steps_for(mcmc * m, double new_steps, int i);
void set_steps_all(mcmc * m, double * new_steps);



#endif

