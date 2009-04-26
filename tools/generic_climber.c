#include "climber.c"

gsl_vector * min;
gsl_vector * max;
gsl_vector_int * round;

void print_real_vector(gsl_vector * x) {
	double realval;
	unsigned int i;

	for (i = 0; i < x->size; i++) {
		realval = gsl_vector_get(x, i) * (gsl_vector_get(max, i)
				- gsl_vector_get(min, i)) + gsl_vector_get(min, i);
		if (gsl_vector_int_get(round, i) == 1) {
			printf("%d", (int) realval);
		} else {
			printf("%f", realval);
		}
		if (i != x->size - 1) {
			printf("\t");
		} else {
			printf("\n");
		}
	}
	fflush(NULL);
}

double f_read_from_stdin(gsl_vector * x) {
	double value;
	print_real_vector(x);

	if (scanf("%lf", &value) != 1) {
		fprintf(stderr, "user aborted.\n");
		exit(1);
	}
	fprintf(stderr, "%f\n", value);
	return value;
}

double f(gsl_vector * x) {
#ifndef NOCACHE
	return f_cached(x, f_read_from_stdin);
#else
	return f_read_from_stdin(x);
#endif
}

void usage(char * progname) {
	printf("%s: SYNOPSIS: <exactness> parameters...\n"
		"\n"
		"\texactness: \thow detailled should we refine the search\n"
		"\tparameters:\teach parameter is a quadrupel of type, min, max, start\n"
		"\t\ttype: i for integer, d for double\n"
		"\t\tmin : lower bound of parameter values\n"
		"\t\tmax : upperbound of parameter values\n"
		"\t\tstart : starting point (0..1)\n"
		"\n"
		"example: %s 0.001 i -10 100 0 d 12.3 54 20\n"
		"\n", progname, progname);
	exit(1);

}

int main(int argc, char ** argv) {
	unsigned int i;
	unsigned int ndim;
	double exactness;
	gsl_vector * start;
	if (argc < 2 + 4 * 1) {
		usage(argv[0]);
		return 1;
	}
	exactness = atof(argv[1]);
	if ((argc - 2) % 4) {
		fprintf(stderr, "%s: wrong number of arguments\n", argv[0]);
		usage(argv[0]);
	}
	ndim = (argc - 2) / 4;
	round = gsl_vector_int_alloc(ndim);
	min = gsl_vector_alloc(ndim);
	max = gsl_vector_alloc(ndim);
	start = gsl_vector_alloc(ndim);

	for (i = 0; i < ndim; i++) {
		if (strcmp(argv[2 + i * 4], "i") == 0)
			gsl_vector_int_set(round, i, 1);
		else if (strcmp(argv[2 + i * 4], "d") == 0)
			gsl_vector_int_set(round, i, 0);
		else
			usage(argv[0]);
		gsl_vector_set(min, i, atof(argv[2 + i * 4 + 1]));
		gsl_vector_set(max, i, atof(argv[2 + i * 4 + 2]));
		gsl_vector_set(start, i, (atof(argv[2 + i * 4 + 3]) - gsl_vector_get(
				min, i)) / (gsl_vector_get(max, i) - gsl_vector_get(min, i)));
	}
	setup_rng();

	find_local_maximum(ndim, exactness, start);
	fprintf(stderr, "--- best: ---");
	print_real_vector(start);

	return 0;
}
