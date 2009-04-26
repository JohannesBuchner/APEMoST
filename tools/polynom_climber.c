#include "climber.c"

gsl_vector ** coefficients = NULL;

gsl_vector ** generate_coefficients(unsigned int dimensions) {
	unsigned int i;
	unsigned int order;
	gsl_vector ** c = (gsl_vector**) calloc(dimensions, sizeof(gsl_vector*));

	for (i = 0; i < dimensions; i++) {
		order = gsl_rng_uniform(get_rng_instance()) * 10 + 1;
		assert(order <= 10);
		c[i] = get_random_uniform_vector(order);
		gsl_vector_add_constant(c[i], -0.5);
		gsl_vector_scale(c[i], 2000);
		IFDEBUG {
			printf("  coeff %d: ", i);
			dump_vectorln(c[i]);
		}
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

char * get_var_name(int n) {
	char * name;
	if (n == 0)
		return "x";
	if (n == 1)
		return "y";
	if (n == 2)
		return "z";
	name = (char *) calloc(100, sizeof(char));
	sprintf(name, "x%d", n);
	return name;
}

void print_coefficients(unsigned int ndim, gsl_vector ** c, gsl_vector * start,gsl_vector * end,
		long count) {
	unsigned int i, j;

	printf("new worst case (%li): start: ", count);
	dump_vectorln(start);
	for (i = 0; i < ndim; i++) {
		printf("  coeff %d: ", i);
		for (j = 0; j < c[i]->size - 1; j++) {
			printf("%f*%s**%d + ", gsl_vector_get(c[i], j), get_var_name(i), j);
		}
		printf("%f*%s**%d\n", gsl_vector_get(c[i], c[i]->size - 1),
				get_var_name(i), j);
	}
	printf("found optimum at: ");
	dump_vectorln(end);
	printf("\n");
}

int main(int argc, char ** argv) {
	unsigned int i;
	long sum = 0;
	long count = 0;
	unsigned int limit;
	unsigned int ndim;
	double exactness;
	gsl_vector * start;
	gsl_vector * end;
	if (argc == 4) {
		exactness = atof(argv[1]);
		limit = atoi(argv[2]);
		ndim = atoi(argv[3]);
	} else {
		printf("%s: SYNOPSIS: <exactness> <number of runs> <ndim>\n"
			"\n"
			"\texactness: \thow detailled should we refine the search\n"
			"\n", argv[0]);
		exit(1);
	}
	setup_rng();
	for (i = 0; i < limit; i++) {
		if (i % 100 == 0) {
			if (coefficients != NULL)
				free_coefficients(ndim);
			coefficients = generate_coefficients(ndim);
		}
		start = get_random_uniform_vector(ndim);
		end = dup_vector(start);
		count = find_local_maximum(ndim, exactness, end);
#ifdef WORSTCASE
		if (count > sum) {
			sum = count;
			print_coefficients(ndim, coefficients, start, end, count);
		}
#else
		sum += count;
#endif
		gsl_vector_free(start);
	}
	gsl_rng_free(get_rng_instance());
	printf("%lu steps\n", sum);
	return 0;
}
