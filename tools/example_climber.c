#include "climber.c"


double f(gsl_vector * x) {
	return -(gsl_vector_get(x, 0)*gsl_vector_get(x, 0) + gsl_vector_get(x, 1)*gsl_vector_get(x, 1));
}

int main(int argc, char ** argv) {
	long count = 0;
	double exactness;
	gsl_vector * start;
	if (argc == 2) {
		exactness = atof(argv[1]);
		dump_d("using exactness", exactness);
	} else {
		printf("%s: SYNOPSIS: <exactness>\n"
			"\n"
			"\texactness: \thow detailled should we refine the search\n"
			"\n", argv[0]);
		exit(1);
	}
	setup_rng();
	start = get_random_uniform_vector(2);
	gsl_vector_add_constant(start, -0.5);
	gsl_vector_scale(start, 2);
	count = find_local_maximum(2, exactness, start);
	printf("%lu steps\n", count);
	gsl_vector_free(start);
	gsl_rng_free(get_rng_instance());
	return 0;
}
