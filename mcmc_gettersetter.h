#ifndef MCMC_GETTERSETTER_H_
#define MCMC_GETTERSETTER_H_

#include "mcmc.h"

/* getter + setter */
const char ** get_params_descr(mcmc * m);
long get_params_accepts_sum(mcmc * m);
long get_params_rejects_sum(mcmc * m);
gsl_vector * get_accept_rate(mcmc * m);
gsl_vector * get_reject_rate(mcmc * m);
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
gsl_vector * get_steps(mcmc * m);

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
void set_params_accepts_for(mcmc * m, long new_params_accept, int i);
void set_params_rejects_for(mcmc * m, long new_params_reject, int i);

#endif /* MCMC_GETTERSETTER_H_ */
