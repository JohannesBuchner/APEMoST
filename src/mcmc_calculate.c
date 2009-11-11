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

#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"

/*#define m->n_par 5*/

void mcmc_append_current_parameters(mcmc * m) {
	mcmc_dump_current(m);
	m->n_iter++;
}

void mcmc_check_best(mcmc * m) {
	if (m->prob > m->prob_best) {
		dump_v("found a better solution", m->params);
		m->prob_best = m->prob;
		set_params_best(m, m->params);
	}
}

