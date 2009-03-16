#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include "debug.h"

#define free(p) { dump_p("about to free", (void*)p); free(p); }

double get_random_number() {
	gsl_rng * random;
	double v;
	
	gsl_rng_env_setup();
	random = gsl_rng_alloc(gsl_rng_default);
	v = gsl_rng_uniform(random);
	gsl_rng_free(random);
	
	return v;
}

#define ALLOCATION_CHUNKS 1
/**
 * for allocation, we don't want to call alloc too often, rather grow in steps
 */
unsigned long get_new_size(unsigned long new_iter) {
	return ((new_iter+1)/ALLOCATION_CHUNKS + 1)*ALLOCATION_CHUNKS;
}

/**
 * make space available as set in m->size
 */
void resize(mcmc * m, unsigned long new_size) {
	unsigned long orig_size = m->size;
	unsigned long i;
	/* TODO: we can not allocate more than the int-space
	 * if we have more iterations than that, we need to move our data to the
	 * disk.
	 */
	dump_ul("resizing params_distr to", new_size);
	if(new_size < orig_size) {
		debug("shrinking -> freeing vectors");
		for(i = orig_size;i > new_size; i--) {
			dump_ul("freeing vector", i - 1);
			gsl_vector_free(m->params_distr[i - 1]);
		}
	}
	m->params_distr = realloc(m->params_distr, new_size);
	
	if(new_size > orig_size) {
		debug("growing -> allocating vectors");
		for(i = orig_size; i < new_size; i++) {
			dump_ul("allocating vector", i);
			m->params_distr[i] = gsl_vector_alloc(m->n_par);
		}
	}
	m->size = new_size;
	debug("done resizing");
}

void prepare_iter(mcmc * m, unsigned long iter) {
	unsigned long new_size = get_new_size(iter);
	if(m->size != new_size) {
		resize(m, new_size);
	}
}

void init_seed(mcmc * m) {
	/* TODO: find right type and initialization for random generator */
	/* TODO: fill arrays with initial values */
	get_random_number();
	m->seed = NULL; 
}

mcmc * mcmc_init(unsigned int n_pars) {
	mcmc * m;
	debug("allocating mcmc struct");
	m = (mcmc*)malloc(sizeof(mcmc));
	assert(m != NULL);
	m->n_par = n_pars;
	m->accept = 0;
	m->reject = 0;
	m->prob = -1e+10;
	m->prob_best = -1e+10;
	
	init_seed(m);
	
	m->params = gsl_vector_alloc(m->n_par);
	assert(m->params != NULL);
	m->params_best = gsl_vector_alloc(m->n_par);
	assert(m->params_best != NULL);
	m->params_distr = NULL;
	prepare_iter(m, 0);
	assert(m->params_distr != NULL);
	
	m->params_accepts = (long*)calloc(m->n_par, sizeof(long));
	assert(m->params_accepts != NULL);
	m->params_rejects = (long*)calloc(m->n_par, sizeof(long));
	assert(m->params_rejects != NULL);
	m->params_step   = gsl_vector_calloc(m->n_par);
	assert(m->params_step != NULL);
	m->params_min    = gsl_vector_calloc(m->n_par);
	assert(m->params_min != NULL);
	m->params_max    = gsl_vector_calloc(m->n_par);
	assert(m->params_max != NULL);

	m->params_descr = NULL;
	m->x_dat = NULL;
	m->y_dat = NULL;
	m->model = NULL;
	return m;
}

void mcmc_free(mcmc * m) {
	m->seed = NULL; /* TODO */
	debug("freeing params");
	gsl_vector_free(m->params);
	debug("freeing params_best");
	gsl_vector_free(m->params_best);
	debug("freeing params_distr");
	resize(m, 0);
	debug("freeing params_descr");
	free(m->params_descr);
	debug("freeing accepts/rejects");
	free(m->params_accepts);
	free(m->params_rejects);
	debug("freeing step/min/max");
	gsl_vector_free(m->params_step);
	gsl_vector_free(m->params_min);
	gsl_vector_free(m->params_max);
	free(m->x_dat);
	free(m->y_dat);
	free(m->model);
	free(m);
}

void mcmc_load_data(mcmc * m, char * filename) {
	(void)filename;
	m->params_descr = NULL;
	m->x_dat = NULL;
	m->y_dat = NULL;
	m->model = NULL;
}

/* TODO: cleanup */

/* TODO: on setting the parameter length, preallocate all 2D arrays. It is a 
 *       slight coding issue when the 2D-arrays remain NULL-pointers.
 */

/* TODO: change from i-1 addressing to i */

/* TODO: all setters seem to need copy functions, check if we really need those
 *       This can probably be done by the caller (using gsl_vector_alloc() and
 *       gsl_vector_memcpy())
 */

/* TODO: check if we need getters. so far they look straight-forward, so 
 *       I'd let people access the struct directly.
 */

extern void calc_model(mcmc * m);
extern void calc_model_for(mcmc * m, int i);

extern double calc_prob(mcmc * m);

void add_values(mcmc * m, int n_iter) {
	/* TODO */
	(void)m;
	n_iter += 0;
}

void write2files(mcmc * m) {
	(void)m;
	/* TODO */
}

void setup(mcmc * m, const char * filename) {
	(void)m;
	(void)filename;
	/* TODO */
}

long get_params_accepts(mcmc * m) {
	unsigned int i;
	long sum = 0;
	for (i = 0; i < m->n_par; i++) {
		sum += m->params_accepts[i];
	}
	return sum;
}
long get_params_rejects(mcmc * m) {
	unsigned int i;
	long sum = 0;
	for (i = 0; i < m->n_par; i++) {
		sum += m->params_rejects[i];
	}
	return sum;
}
long get_params_accepts_for(mcmc * m, int i) {
	return m->params_accepts[i];
}
long get_params_rejects_for(mcmc * m, int i) {
	return m->params_rejects[i];
}

gsl_histogram * get_hist(mcmc * m, int index, int nbins) {
	return calc_hist(m->params_distr, index, nbins);
}


void set_params_accepts_for(mcmc * m, long new_params_accept, int i) {
	m->params_accepts[i] = new_params_accept;
}
void set_params_rejects_for(mcmc * m, long new_params_reject, int i) {
	m->params_rejects[i] = new_params_reject;
}

void set_prob_best(mcmc * m, double new_prob_best) {
	m->prob_best = new_prob_best;
}

void set_minmax_for(mcmc * m, double new_min, double new_max, int i) {
	gsl_vector_set(m->params_min, i, new_min);
	gsl_vector_set(m->params_max, i, new_max);
}

void set_steps_for(mcmc * m, double new_step, int i) {
	gsl_vector_set(m->params_step, i, new_step);
}

void set_model(mcmc * m, gsl_vector * new_model) {
	m->model = new_model;
}

void set_n_par(mcmc * m, int new_n_par) {
	m->n_par = new_n_par;
}

void set_params_best(mcmc * m, gsl_vector * new_params_best) {
	m->params_best = new_params_best;
}

void set_params_for(mcmc * m, double new_param, int i) {
	gsl_vector_set(m->params, i, new_param);
}

void set_params(mcmc * m, gsl_vector * new_params) {
	m->params = new_params;
}

void set_params_descr_all(mcmc * m, char ** new_par_descr) {
	m->params_descr = new_par_descr;
}

void set_params_descr_for(mcmc * m, char * new_par_descr, int i) {
	m->params_descr[i-1] = new_par_descr;
}

void set_seed(mcmc * m, gsl_vector * new_seed) {
	m->seed = new_seed;
}

void set_probability(mcmc * m, double new_prob) {
	m->prob = new_prob;
}

void set_x(mcmc * m, gsl_vector * new_x) {
	/*if(m->x_dat != NULL)
		gsl_vector_free(m->x_dat);*/
	m->x_dat = new_x;
}

/* TODO: check if we need that
void set_x_copy(mcmc * m, gsl_vector * new_x) {
	gsl_vector_memcpy(m->x_dat, new_x);
}*/

void set_y(mcmc * m, gsl_vector * new_y) {
	/*if(m->x_dat != NULL)
		gsl_vector_free(m->x_dat);*/
	m->y_dat = new_y;
}

/* TODO: check if we need that
void set_y_copy(mcmc * m, gsl_vector * new_y) {
	gsl_vector_memcpy(m->y_dat, new_y);
}*/

void free_gsl_vector_array(gsl_vector ** arr) {
	int i = 0;
	if(arr != NULL) {
		while(arr[i]!=NULL)
			gsl_vector_free(arr[i]);
	}
}

/*
gsl_vector * alloc_gsl_vector_array(unsigned int size) {
	int i = 0;
	arr = (gsl_vector *)calloc(size, sizeof(gsl_vector *));
	assert(arr != NULL);
	for(i = 0; i < size; i++) {
		arr[i] = gsl_vector_alloc(n);
		gsl_vector_memcpy(arr[i], src[i]);
	}
	return arr;
}
*/
/*
gsl_vector ** copy_gsl_vector_array(gsl_vector ** arr, const gsl_vector ** src, size_t n) {
	int i = 0;
	assert(src != NULL);
	if(arr == NULL) {
		arr = (gsl_vector *)calloc(n, sizeof(gsl_vector *));
		assert(arr != NULL);
		for(i = 0; i < n; i++) {
			arr[i] = gsl_vector_alloc(n);
			gsl_vector_memcpy(arr[i], src[i]);
		}
		return arr;
	}else{
		for(i = 0; i < n; i++) {
			gsl_vector_memcpy(arr[i], src[i]);
		}
	}
}*/

void set_steps_all(mcmc * m, double * new_steps) {
	unsigned int i;
	for(i = 1; i < m->n_par + 1; i++) {
		set_steps_for(m, new_steps[i+1], i);
	}
}

