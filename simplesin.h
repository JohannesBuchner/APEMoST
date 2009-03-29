#ifndef SIMPLESIN_H_
#define SIMPLESIN_H_
#include "mcmc.h"

typedef struct {
	/**
	 * likelihood of acceptance (TODO: yes?)
	 */
	double beta;

} parallel_tempering_mcmc;

void set_beta(mcmc * m, double newbeta) {
	((parallel_tempering_mcmc *) m->additional_data)->beta = newbeta;
}
double get_beta(mcmc * m) {
	return ((parallel_tempering_mcmc *) m->additional_data)->beta;
}

#endif /* SIMPLESIN_H_ */
