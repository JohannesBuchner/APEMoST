#ifndef MCMC
#define MCMC

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_histogram.h>
#include <gsl/gsl_vector.h>

gsl_histogram * get_hist(gsl_vector ** param_distr, int index, int nbins);


#endif

