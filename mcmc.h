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
	int n_par; 
	/* number of accepted steps for MCMC (after calibration) */
	long accept;
	/* number of rejected steps for MCMC (after calibration) */
	long reject;
	/* probability of the most recently evaluated parameter values */
	double prob;
	/* probability of best parameter values yets */
	double prob_best;
	/* pointer to 1D-array as seed for random number generators */
	gsl_vector * seed;
	/* pointer to 1D-array containing the parameter values */
	gsl_vector * params;
	/* pointer to 1D-array containing the best parameters yet */
	gsl_vector * params_best;
	/* pointer to 2D-array containing the resulting values for each parameter */
	gsl_vector ** params_distr;
	/* pointer to 2D-array which acts as a buffer for storing parameter values 
	 * for quickly */
	gsl_vector ** params_distr_buf;
	/* pointer to 1D-array containing description of parameters in string 
	 * format */
	char ** params_descr;
	/* pointer to 1D-array containing the current value of the prior probability
	 * for each parameter */
	gsl_vector * params_priors;
	/* pointer to 2D-array containing the number of accepted/rejected steps for 
	 * individual parameters*/
	gsl_vector ** params_ar;
	/* pointer to 1D-array containing the current step width for individual 
	 * parameters */
	gsl_vector * params_step;
	/* pointer to 2D-array containing the lower and upper limits for each 
	 * parameter */
	gsl_vector ** params_minmax;
	/* pointer to 1D-array containing the abscissa of the observations (e.g., 
	 * date, frequency, wavelength) */
	gsl_vector * x_dat; 
	/* pointer to 1D-array containing the ordinate of the observations (e.g., 
	 * intensity, power, radial vel.) */
	gsl_vector * y_dat;
	/* pointer to 1D-array containing the model values corresponding to the 
	 * observed abscissa values */
	gsl_vector * model;
} mcmc_data;




#endif

