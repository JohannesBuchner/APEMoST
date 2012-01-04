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

#ifndef PARALLEL_TEMPERING_CONFIG_H_
#define PARALLEL_TEMPERING_CONFIG_H_

#include "mcmc.h"
#include "parallel_tempering_beta.h"

#define CALIBRATION_FILE "calibration_results"

void write_params_file(mcmc * m);

void write_calibration_summary(mcmc ** chains, unsigned int n_chains);

mcmc ** setup_chains();

void read_calibration_file(mcmc ** chains, unsigned int n_chains);

void write_calibrations_file(mcmc ** chains, const unsigned int n_chains);

#endif /* PARALLEL_TEMPERING_CONFIG_H_ */
