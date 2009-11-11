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

#ifndef TEMPERING_INTERACTION_H_
#define TEMPERING_INTERACTION_H_

#include "mcmc.h"

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * this shouldn't make any difference at the moment ...
 */
#define RANDOMSWAP
#endif

/**
 * does swapping chains and other mixes.
 *
 * @param chains
 * @param iter number of iterations that passed. Is a multiple of n_swap.
 * @param n_beta size of chains
 */
void tempering_interaction(mcmc ** chains, unsigned int n_beta,
		unsigned long iter);

#endif /* TEMPERING_INTERACTION_H_ */
