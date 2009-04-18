#include "climber.c"

gsl_vector ** coefficients;

gsl_vector ** generate_coefficients(unsigned int dimensions) {
	unsigned int i;
	unsigned int order;
	gsl_vector ** c = (gsl_vector**) calloc(dimensions, sizeof(gsl_vector*));

	for (i = 0; i < dimensions; i++) {
		order = gsl_rng_uniform(get_rng_instance()) * 10 + 1;
		assert(order <= 10);
		c[i] = get_random_uniform_vector(order);
		gsl_vector_add_constant(c[i], -0.5);
		gsl_vector_scale(c[i], 20);
		printf("  coeff %d: ", i);
		dump_vectorln(c[i]);
	}
	return c;
}

double my_polynome(gsl_vector * x) {
	double v = 1;
	double w;
	unsigned int i, j;
	for (i = 0; i < x->size; i++) {
		w = 0;
		for (j = 0; j < coefficients[i]->size; j++) {
			w += pow(gsl_vector_get(x, i), j) * gsl_vector_get(coefficients[i],
					j);
		}
		v *= w;
	}
	return v;
}

double f(gsl_vector * x) {
	return my_polynome(x);
}
void free_coefficients(unsigned int ndim) {
	unsigned int i;
	for (i = 0; i < ndim; i++) {
		gsl_vector_free(coefficients[i]);
	}
	free(coefficients);
}

int main(int argc, char ** argv) {
	unsigned int i;
	long count = 0;
	unsigned int limit;
	unsigned int ndim;
	double exactness;
	gsl_vector * start;
	if (argc == 4) {
		limit = atoi(argv[1]);
		ndim = atoi(argv[2]);
		exactness = 1.0 / atoi(argv[3]);
	} else {
		printf("%s: SYNOPSIS: <number of runs> <ndim> <exactness>\n"
			"\n"
			"\n\texactness\tinverse. so if you want 0.001, give 1000."
			"\n", argv[0]);
		exit(1);
	}
	setup_rng();
	for (i = 0; i < limit; i++) {
		coefficients = generate_coefficients(ndim);
		start = get_random_uniform_vector(ndim);
		count += find_local_maximum(ndim, exactness, start);
		free_coefficients(ndim);
		gsl_vector_free(start);
	}
	gsl_rng_free(get_rng_instance());
	printf("%lu steps\n", count);
	return 0;
}
