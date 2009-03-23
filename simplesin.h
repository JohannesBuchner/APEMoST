#ifndef SIMPLESIN_H_
#define SIMPLESIN_H_
#include "mcmc.h"

typedef struct {
	mcmc * m;
	/**
	 * likelihood of acceptance (TODO: yes?)
	 */
	double beta;

} parallel_tempering_mcmc;

void set_beta(parallel_tempering_mcmc * m, double newbeta) {
	m->beta = newbeta;
}
double get_beta(parallel_tempering_mcmc * m) {
	return m->beta;
}


#endif /* SIMPLESIN_H_ */
