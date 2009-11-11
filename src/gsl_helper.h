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

#ifndef MCMC_GSL_HELPER
#define MCMC_GSL_HELPER

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif

/**
 * sums the values
 */
double calc_vector_sum(const gsl_vector * v);

/**
 * sums the squared values
 */
double calc_vector_squaresum(const gsl_vector * v);
/**
 * returns a duplicate.
 *
 * The caller has to free the returned vector.
 */
gsl_vector * dup_vector(const gsl_vector * v);

/**
 * normalizes the vector, i.e. the values are scaled so that the sum of
 * all values is 1
 *
 * The caller has to free the returned vector.
 */
gsl_vector * calc_normalized(const gsl_vector * v);

/**
 * @return 1 if vectors contain the same entries
 */
int calc_same(const gsl_vector * a, const gsl_vector * b);

/**
 * a = max(a, b)
 */
void max_vector(gsl_vector * a, const gsl_vector * b);
/**
 * a = min(a, b)
 */
void min_vector(gsl_vector * a, const gsl_vector * b);




/**
 * sorts all vectors by the entries in the first vector.
 *
 * selection sort
 */
void sort(gsl_vector ** vs, unsigned int nvectors, unsigned int vector_size);

/**
 * given the first index of the matrix, iterate through the second to find
 * the smallest entry
 */
double min_column(const gsl_matrix * m, const unsigned int i);

/**
 * given the second index of the matrix, iterate through the first to find
 * the smallest entry
 */
double min_row(const gsl_matrix * m, const unsigned int i);

/**
 * given the first index of the matrix, iterate through the second to find
 * the largest entry
 */
double max_column(const gsl_matrix * m, const unsigned int i);

/**
 * given the second index of the matrix, iterate through the first to find
 * the largest entry
 */
double max_row(const gsl_matrix * m, const unsigned int i);

/**
 * n-dimensional weighted square deviation
 */
double calc_deviation(const gsl_matrix * x, const gsl_vector * y,
		const gsl_vector * k, const double d, const gsl_vector * weights);

/**
 * n-dimensional weighted linear regression
 * returns k and d.
 */
gsl_vector * linreg_n(const gsl_matrix * x, const gsl_vector * y, double * d,
		const gsl_vector * weights);

#endif
