#ifndef MCMC_GETTERSETTER_H_
#define MCMC_GETTERSETTER_H_

#include "mcmc.h"

/* getter + setter */
const char ** get_params_descr(const mcmc * m);
long get_params_accepts_sum(const mcmc * m);
long get_params_rejects_sum(const mcmc * m);
gsl_vector * get_accept_rate(const mcmc * m);
long get_params_accepts_for(const mcmc * m, int i);
long get_params_rejects_for(const mcmc * m, int i);
gsl_vector * get_params(const mcmc * m);
gsl_vector * get_params_best(const mcmc * m);
const gsl_vector * get_x(const mcmc * m);
const gsl_vector * get_y(const mcmc * m);
unsigned int get_n_par(const mcmc * m);
gsl_rng * get_random(const mcmc * m);
double get_prob(const mcmc * m);
double get_prob_best(const mcmc * m);
gsl_vector * get_steps(const mcmc * m);

void set_prob(mcmc * m, const double new_prob);
void set_prob_best(mcmc * m, const double new_prob_best);
void set_minmax_for(mcmc * m, const double new_min, const double new_max,
		const unsigned int i);
void set_model(mcmc * m, gsl_vector * new_model);
void set_n_par(mcmc * m, const int new_n_par);
void set_params_best(mcmc * m, gsl_vector * new_params_best);
void set_params_for(mcmc * m, const double new_param, const unsigned int i);
void set_params(mcmc * m, gsl_vector * new_params);
void set_params_descr_all(mcmc * m, const char ** new_par_descr);
void set_params_descr_for(mcmc * m, const char * new_par_descr,
		const unsigned int i);
void set_random(mcmc * m, gsl_rng * newrandom);
void set_prob(mcmc * m, const double new_prob);
void set_x(mcmc * m, const gsl_vector * new_x);
void set_x_copy(mcmc * m, const gsl_vector * new_x);
void set_y(mcmc * m, const gsl_vector * new_y);
void set_y_copy(mcmc * m, const gsl_vector * new_y);
void set_steps_for(mcmc * m, const double new_steps, int i);
void set_steps_all(mcmc * m, const double * new_steps);
void set_params_accepts_for(mcmc * m, const long new_params_accept,
		const unsigned int i);
void set_params_rejects_for(mcmc * m, const long new_params_reject,
		const unsigned int i);

void inc_params_accepts_for(mcmc * m, const unsigned int i);
void inc_params_rejects_for(mcmc * m, const unsigned int i);
void inc_params_accepts(mcmc * m);
void inc_params_rejects(mcmc * m);

/**
 * get next random number (uniformly distributed between 0 and 1)
 */
double get_next_urandom(const mcmc * m);

/**
 * get next random number (logarithmic uniformly distributed between -inf and 0)
 */
double get_next_alog_urandom(const mcmc * m);

#endif /* MCMC_GETTERSETTER_H_ */