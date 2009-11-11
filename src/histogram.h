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

#ifndef HISTOGRAM
#define HISTOGRAM

/**
 * create a inclusive histogram where min and max are inside the value range.
 */
gsl_histogram * create_hist(int nbins, double min, double max);

/**
 * calculate a histogram
 * @param v vector to look at
 * @param nbins number of bins to use for the histogram
 */
gsl_histogram * calc_hist(const gsl_vector * v, int nbins);

/**
 * append the input of a file with n columns to the according histograms.
 */
void append_to_hists(gsl_histogram ** hists, unsigned int n,
		const char * filename);

/**
 * find the smallest and largest values of a file for each column, then set
 * the min/max vectors
 */
void find_min_max(char * filename, gsl_vector * min, gsl_vector * max);

/**
 * find the smallest and largest values of a file for each column, then update
 * the min/max vectors
 */
void update_min_max(char * filename, gsl_vector * min, gsl_vector * max);

#endif
