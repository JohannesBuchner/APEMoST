#ifndef MCMC_STRUCT_H_
#define MCMC_STRUCT_H_

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_rng.h>

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
	 * random number generator
	 */
	gsl_rng * random;
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
	const char ** params_descr;
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
	const gsl_vector * x_dat;
	/**
	 * pointer to 1D-array containing the ordinate of the observations (e.g.,
	 * intensity, power, radial vel.).
	 * size = x-size.
	 */
	const gsl_vector * y_dat;
	/** pointer to 1D-array containing the model values corresponding to the
	 * observed abscissa values.
	 * This is a synthetic (calculated) y_dat value for the current parameters.
	 * size = x-size
	 */
	gsl_vector * model;

	/** number of iterations for which space is allocated (params_distr) */
	unsigned long size;
	/** number of iterations calculated */
	unsigned long n_iter;
} mcmc;

#endif /* MCMC_STRUCT_H_ */
