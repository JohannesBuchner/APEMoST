#ifndef MCMC
#define MCMC

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>


typedef struct {
	/* number of parameters */
	unsigned int n_par; 
	/* number of accepted steps for MCMC (after calibration) */
	unsigned long accept;
	/* number of rejected steps for MCMC (after calibration) */
	unsigned long reject;
	/* probability of the most recently evaluated parameter values */
	double prob;
	/* probability of best parameter values yets */
	double prob_best;
	/* pointer to 1D-array as seed for random number generators */
	gsl_vector * seed;
	/* pointer to 1D-array containing the parameter values */
	gsl_vector * params; /* size = n_par */
	/* pointer to 1D-array containing the best parameters yet */
	gsl_vector * params_best;  /* size = n_par */
	/* pointer to 2D-array containing the resulting values for each parameter */
	gsl_vector ** params_distr; /* size = iter, n_par */
	/* pointer to 1D-array containing description of parameters in string 
	 * format */
	char ** params_descr; /* size = n_par */
	/* pointer to 1D-array containing the number of accepted/rejected steps for 
	 * individual parameters*/
	long * params_accepts; /* size = n_par */
	/* pointer to 1D-array containing the number of accepted/rejected steps for 
	 * individual parameters*/
	long * params_rejects; /* size = n_par */
	/* pointer to 1D-array containing the current step width for individual 
	 * parameters */
	gsl_vector * params_step; /* size = n_par; set by calibration */
	/* pointer to 1D-array containing the lower limit for each parameter */
	gsl_vector * params_min; /* size = n_par */
	/* pointer to 1D-array containing the upper limit for each parameter */
	gsl_vector * params_max; /* size = n_par */
	/* pointer to 1D-array containing the abscissa of the observations (e.g., 
	 * date, frequency, wavelength) */
	gsl_vector * x_dat; /* size = x-size */
	/* pointer to 1D-array containing the ordinate of the observations (e.g., 
	 * intensity, power, radial vel.) */
	gsl_vector * y_dat; /* size = x-size */
	/* pointer to 1D-array containing the model values corresponding to the 
	 * observed abscissa values */
	gsl_vector * model; /* size = x-size; synthetic y_dat value */
	
	/* number of iterations for which space is allocated (params_distr) */
	unsigned long size;
	/* number of iterations calculated */
	unsigned long iter;
} mcmc;

double get_random_number();
mcmc * mcmc_init(unsigned int n_pars);
mcmc * mcmc_load(const char * filename);
void mcmc_free(mcmc * m);
void add_values(mcmc * m, int n_iter);
void write2files(mcmc * m);
long get_params_accepts(mcmc * m);
long get_params_rejects(mcmc * m);
long get_params_accepts_for(mcmc * m, int i);
long get_params_rejects_for(mcmc * m, int i);
double get_params_ar_for(mcmc * m, int i);
gsl_histogram * get_hist(mcmc * m, int i, int nbins);
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
gsl_vector ** copy_gsl_vector_array(gsl_vector ** arr, const gsl_vector ** src, size_t n);
void set_steps_all(mcmc * m, double * new_steps);



#endif

