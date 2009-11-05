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
