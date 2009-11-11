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

/**
 * Here are the "private" methods of the class and helper functions
 */
#ifndef MCMC_INTERNAL_H_
#define MCMC_INTERNAL_H_

#include "mcmc.h"
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_sf.h>

/**
 * create class
 * \private
 * @param n_pars parameters
 */
mcmc * mcmc_init(const unsigned int n_pars);

/**
 * count the lines (\n) in the file
 * @param filename
 */
unsigned int countlines(const char * filename);

/**
 * a modulo operator for double values
 */
/*double mod_double(const double x, const double div);*/
#define mod_double(x, div) 	((x) < 0 ? \
	(x) - (div) * (int) ((x) / (div) - 1) : \
	(x) - (div) * (int) ((x) / (div)))

/**
 * the value with positive sign.
 */
/*double abs_double(const double x);*/
#define abs_double(x) 	((x) < 0 ? -(x) : (x))

#endif /* MCMC_INTERNAL_H_ */
