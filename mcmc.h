#ifndef MCMC
#define MCMC

#include <stdio.h>
#include <stdlib.h>

#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>

/**
 * The main class of operation.
 */
typedef struct {
	/** number of parameters */
	unsigned int n_par;
	/** number of accepted steps for MCMC (after calibration) */
	unsigned long accept;
	/** number of rejected steps for MCMC (after calibration) */
	unsigned long reject;
	/** probability of the most recently evaluated parameter values */
	double prob;
	/** probability of best parameter values yet */
	double prob_best;
	/**
	 * seed for random number generators
	 */
	gsl_vector * seed;
	/**
	 * current parameters
	 * size = n_par
	 */
	gsl_vector * params;
	/**
	 * best parameters yet
	 * size = n_par
	 */
	gsl_vector * params_best;
	/**
	 * pointer to 2D-array containing the resulting values for each parameter.
	 * for each iteration, and parameter, the values of the calculation are
	 * kept.
	 * size = iter, n_par
	 */
	gsl_vector ** params_distr;
	/**
	 * descriptions of parameters
	 * size = n_par
	 */
	char ** params_descr;
	/**
	 * number of accepted steps for individual parameters
	 * size = n_par
	 */
	long * params_accepts;
	/**
	 * number of rejected steps for individual parameters
	 * size = n_par
	 */
	long * params_rejects;
	/**
	 * current step widths for individual parameters
	 * size = n_par; set by calibration
	 */
	gsl_vector * params_step;
	/**
	 * lower limits for each parameter
	 * size = n_par
	 */
	gsl_vector * params_min;
	/**
	 * upper limits for each parameter
	 * size = n_par
	 */
	gsl_vector * params_max;
	/**
	 * pointer to 1D-array containing the abscissa of the observations (e.g.,
	 * date, frequency, wavelength)
	 * size = x-size
	 */
	gsl_vector * x_dat;
	/**
	 * pointer to 1D-array containing the ordinate of the observations (e.g.,
	 * intensity, power, radial vel.).
	 * size = x-size.
	 */
	gsl_vector * y_dat;
	/** pointer to 1D-array containing the model values corresponding to the
	 * observed abscissa values.
	 * This is a synthetic (calculated) y_dat value for the current parameters.
	 * size = x-size
	 */
	gsl_vector * model;

	/** number of iterations for which space is allocated (params_distr) */
	unsigned long size;
	/** number of iterations calculated */
	unsigned long iter;
} mcmc;

/**
 * \private
 */
double get_random_number();
/**
 * create class
 * \private
 * @param n_pars parameters
 */
mcmc * mcmc_init(unsigned int n_pars);
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
 * prepare the calculation of (next) iteration, i.e., allocate space
 * @param iter number of iteration
 */
void mcmc_prepare_iteration(mcmc * m, unsigned long iter);

/**
 * was: add_values
 */
void mcmc_append_current_parameters(mcmc * m, int n_iter);

void mcmc_dump_model(mcmc * m);
void mcmc_dump_y_dat(mcmc * m);

/**
 * @see calc_hist
 */
gsl_histogram * get_hist(mcmc * m, int i, int nbins);

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
gsl_vector ** copy_gsl_vector_array(gsl_vector ** arr, const gsl_vector ** src, size_t n);
void set_steps_all(mcmc * m, double * new_steps);



#endif

