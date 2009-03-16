#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "gsl_helper.h"
#include <gsl/gsl_rng.h>
#include "debug.h"

#define free(p) { IFSEGV dump_p("about to free", (void*)p); free(p); }
#define NCOLUMNS 5
#define MAX_LINE_LENGTH 256

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
	IFSEGV debug("allocating mcmc struct");
	m = (mcmc*)malloc(sizeof(mcmc));
	assert(m != NULL);
	m->iter = 0;
	m->size = 0;
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

	m->params_descr = (char**)calloc(m->n_par, sizeof(char*));;
	m->x_dat = NULL;
	m->y_dat = NULL;
	m->model = NULL;
	IFSEGV debug("allocating mcmc struct done");
	return m;
}

void mcmc_free(mcmc * m) {
	unsigned int i;
	
	m->seed = NULL; /* TODO */
	IFSEGV debug("freeing params");
	gsl_vector_free(m->params);
	IFSEGV debug("freeing params_best");
	gsl_vector_free(m->params_best);
	IFSEGV debug("freeing params_distr");
	resize(m, 0);
	
	IFSEGV debug("freeing params_descr");
	for (i = 0; i < m->n_par; i++) {
		free(m->params_descr[i]);
	}
	free(m->params_descr);
	
	IFSEGV debug("freeing accepts/rejects");
	free(m->params_accepts);
	free(m->params_rejects);
	IFSEGV debug("freeing step/min/max");
	gsl_vector_free(m->params_step);
	gsl_vector_free(m->params_min);
	gsl_vector_free(m->params_max);
	free(m->x_dat);
	free(m->y_dat);
	free(m->model);
	free(m);
}
#undef SKIP_DEBUG
#define SKIP_DEBUG 0

FILE * openfile(const char * filename) {
	FILE * input = fopen(filename, "r");
	if(input == NULL) {
		fprintf(stderr, "error opening file %s\n", filename);
		perror("file could not be opened");
		exit(1);
	}
	return input;
}

unsigned int countlines(const char * filename) {
	int nlines = 0;
	int c;
	FILE * input = openfile(filename);
	while(1) {
		c = fgetc(input);
		if(c == '\n') {
			nlines++;
		}
		if(c == EOF)
			break;
	}
	fclose(input);
	return nlines;
}

int strnlen(char * s, int maxlen) {
	int i;
	for(i=0;i<maxlen && s[i]!=0;i++);
	return i;
}

/**
 * returns 0 on success.
 */
int load_parameter(mcmc * m, FILE * input, int i) {
	int col = 0;
	double start;
	double min;
	double max;
	double step;
	char * descr = (char*)calloc(MAX_LINE_LENGTH, sizeof(char));
	dump_i("parsing line", i);
	
	col = fscanf(input, "%lf\t%lf\t%lf\t%s\t%lf\n", &start, &min, &max, descr, &step);
	if(col != 5) {
		fprintf(stderr, "only %d fields matched.\n", col);
		return 1;
	}
	if(step < 0 || step > 1) {
		fprintf(stderr, "start (column 1) must be between 0 and 1. currently: %f.\n", start);
		return 1;
	}
	if(!(descr != NULL && 
		strnlen(descr, MAX_LINE_LENGTH) > 0 && 
		strnlen(descr, MAX_LINE_LENGTH) < MAX_LINE_LENGTH
	)) {
		fprintf(stderr, "description invalid: %s\n", descr);
		return 1;
	}
	debug("setting values");
	gsl_vector_set(m->params, i, start);
	gsl_vector_set(m->params_min,  i, min);
	gsl_vector_set(m->params_max,  i, max);
	m->params_descr[i] = descr;
	gsl_vector_set(m->params_step, i, min + start*(max-min));
	debug("setting values done.");
	return 0;
}

int load_datapoint(mcmc * m, FILE * input, int i) {
	double x;
	double y;
	int col;
	
	col = fscanf(input, "%lf%lf", &x, &y);
	if(col != 2) {
		fprintf(stderr, "only %d fields matched.\n", col);
		return 1;
	}
	gsl_vector_set(m->x_dat, i, x);
	gsl_vector_set(m->y_dat, i, y);
	return 0;
}

void load_data(mcmc * m, const char * filename) {
	FILE * input;
	int i = 0;
	int npoints = countlines(filename);
	dump_i("lines", npoints);
		
	m->x_dat = gsl_vector_alloc(npoints);
	m->y_dat = gsl_vector_alloc(npoints);
	m->model = gsl_vector_alloc(npoints);
	
	input = openfile(filename);
	for (i = 0; i < npoints; i++) {
		if(load_datapoint(m, input, i) != 0) {
			fprintf(stderr, "Line %d of %s is of incorrect format.", 
				i + 1, filename);
			exit(1);
		}
	}
	dump_i("loaded data points", npoints);
	
	
}

char *strdup(const char *s) {
	char *buf = calloc(strlen(s)+1, sizeof(char));
	if (buf != NULL)
		strcpy(buf, s);
	return buf;
}

mcmc * mcmc_load(const char * filename) {
	mcmc * m;
	FILE * input;
	char datafilename[MAX_LINE_LENGTH];
	char datafilepath[MAX_LINE_LENGTH];
	char * datadir;
	int currentline = 0;
	const int pretext = 1;
	
	int nlines;
	nlines = countlines(filename);
	dump_i("number of lines", nlines);
	
	m = mcmc_init(nlines - pretext);
	input = openfile(filename);
	fscanf(input, "%s\n", datafilename);
	currentline++;
	
	while(currentline<nlines) {
		if(load_parameter(m, input, currentline - pretext) != 0) {
			fprintf(stderr, "Line %d of %s is of incorrect format.", 
				currentline + 1, filename);
			exit(1);
		}
		currentline++;
	}
	fclose(input);
	
	datadir = dirname(strdup(filename));
	sprintf(datafilepath, "%s/%s", datadir, datafilename);
	dump_s("looking for data in file", datafilepath);
	
	load_data(m, datafilepath);
	return m;
}

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
	m->params_descr[i] = new_par_descr;
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

