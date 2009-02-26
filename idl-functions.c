

#include "mcmc.h"

/**
 * PTR_NEW: 
 *   creates a pointer to a heap variable
 *   argument is copied and a pointer is returned.
 **
 * DBLARR: 
 *   creates a double/float array/vectro
 *   argument: dimension, up to 8
 **
 * (*var)[a,b]
 *   dereferences pointer, accesses index a of var, and of that accesses index b
 * (*var)[a,*]
 *   full vector of 2D array (*var), var being a pointer to the 2D-array
 ** 
 * ?INDGEN:
 *   generates a array where each entry has the value of its subscript.
 *   examples: findgen (float), dindgen (double)
 * 
 **
 * Histograms:
 *   count stuff that falls into bins.
 *   return the number of elements in the bin for each bin.
 *   usual arguments: nbins, min, max
 */

#define IFDEBUG if(0)
#define DEBUG(str)     IFDEBUG printf(str "\n");
#define DEBUGD(str, v) IFDEBUG printf(str ": %f\n", v);

gsl_histogram * get_hist(gsl_vector ** param_distr, int index, int nbins) {
	double max;
	double min;
	unsigned int i;
	double binwidth;
	double sum = 0;
	double v;
	gsl_histogram * h;
	gsl_vector * dat = param_distr[index-1];
	
	gsl_vector_minmax (dat, &min, &max);
	binwidth = (max - min)/nbins;
	DEBUGD("min", min);
	DEBUGD("max", max);
	
	DEBUG("allocating the histogram");
	h = gsl_histogram_alloc (dat->size);
	DEBUG("setting range");
	gsl_histogram_set_ranges_uniform (h, min, max);
	
	/* with out the following, the max element doesn't fall in the last bin */
	h->range[h->n] += 1; 
	
	DEBUG("summing up");
	for(i=0; i<dat->size; i++){ 
		v = gsl_vector_get (dat, i);
		sum += v;
		gsl_histogram_increment (h, v);
	}
	DEBUG("scaling");
	/* double gsl_histogram_sum (const gsl_histogram * h) */
	gsl_histogram_scale (h, 1/sum);
	DEBUG("done");
	return h;
}



