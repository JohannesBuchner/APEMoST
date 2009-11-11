/*
    APEMoST - Automated Parameter Estimation and Model Selection Toolkit
    Copyright (C) 2009  Johannes Buchner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef MCMC_STRUCT_H_
#define MCMC_STRUCT_H_

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
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
	/** explicit prior, so it can be subtracted from prob */
	double prior;
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
	 * files where visited nodes are written to.
	 */
	FILE ** files;
	/**
	 * descriptions of parameters
	 * size = n_par
	 */
	const char ** params_descr;
	/**
	 * number of accepted steps for individual parameters
	 * size = n_par
	 */
	unsigned long * params_accepts;
	/**
	 * number of rejected steps for individual parameters
	 * size = n_par
	 */
	unsigned long * params_rejects;
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
	 * arbitrary sized array containing the observation as found in the
	 * file "data"
	 *
	 * column 0 is the x-data
	 * column 1 is the y-data
	 * etc.
	 */
	const gsl_matrix * data;

	/** number of iterations calculated */
	unsigned long n_iter;

	/** any extensions can be stored here **/
	void * additional_data;
} mcmc;

#endif /* MCMC_STRUCT_H_ */
