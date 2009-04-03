#ifndef TEMPERING_INTERACTION_H_
#define TEMPERING_INTERACTION_H_

#include "mcmc.h"

void tempering_interaction(mcmc ** sinmod, unsigned int n_beta,
		unsigned long iter);

#endif /* TEMPERING_INTERACTION_H_ */
