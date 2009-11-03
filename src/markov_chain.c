#include <string.h>
#include <stdio.h>
#include <libgen.h>

#include "mcmc.h"
#include "mcmc_internal.h"
#include "debug.h"
#include "gsl_helper.h"
#include <gsl/gsl_sf.h>

#define BETWEEN(x, min, max) ( (x) >= (min) && (x) <= (max) )
#define MAX(a, b) ((a) > (b) ? a : b)
#define MIN(a, b) ((a) < (b) ? a : b)

static void restart_from_best(mcmc * m) {
	set_params(m, dup_vector(get_params_best(m)));
	set_prob(m, get_prob_best(m));
}

void burn_in(mcmc * m, const unsigned int burn_in_iterations) {
	unsigned long iter;
	unsigned long subiter;

	gsl_vector * original_steps = get_steps(m);
	m->params_step = dup_vector(m->params_max);
	gsl_vector_sub(m->params_step, m->params_min);
	gsl_vector_scale(m->params_step, 0.1);

	debug("Beginning calibration of MCMC ...");
	debug("Starting burn-in ...");
	mcmc_check(m);
	for (iter = 0; iter < burn_in_iterations / 2; iter++) {
		for (subiter = 0; subiter < 200; subiter++) {
			markov_chain_step(m);
		}
		iter += subiter;
		dump_ul("\tBurn-in Iteration", iter);
		IFVERBOSE {
			dump_v("stepwidth", get_steps(m));
			dump_v("params", get_params(m));
		}
		mcmc_check_best(m);
	}
	debug("Re-initializing burn-in ...");
	restart_from_best(m);
	gsl_vector_scale(m->params_step, 0.5);
	for (; iter < burn_in_iterations; iter++) {
		for (subiter = 0; subiter < 200; subiter++) {
			markov_chain_step(m);
		}
		iter += subiter;
		dump_ul("\tBurn-in Iteration", iter);
		IFVERBOSE {
			dump_v("stepwidth", get_steps(m));
			dump_v("params", get_params(m));
		}
		mcmc_check_best(m);
	}
	debug("Burn-in done, adjusting steps ...");
	gsl_vector_memcpy(get_steps(m), original_steps);
	gsl_vector_free(original_steps);
	mcmc_check(m);
	debug("Burn-in done.");

}

void clear_bit(char * bitfield, unsigned int i) {
	bitfield[i / 8] &= ~(1 << (i % 8));
}
void set_bit(char * bitfield, unsigned int i) {
	bitfield[i / 8] |= (1 << (i % 8));
}
int get_bit(char * bitfield, unsigned int i) {
	return bitfield[i / 8] & (1 << (i % 8));
}

#ifndef ACCURACY_DEVIATION_FACTOR
/**
 * How good should the acceptance rate be calculated in dependence of
 * deviation from the desired value?
 * accuracy = factor * deviation
 */
#define ACCURACY_DEVIATION_FACTOR 0.25
#endif
/**
 * Get acceptance rate.
 * The closer the acceptance rate is to the desired acceptance rate, the more
 * accurately will it be assessed.
 *
 * @param min_accuracy you can request a upper limit on the accuracy, e.g. 1%.
 * 		any calculation will have at most the accuracy of 1% then. Otherwise put 0 here.
 * @param max_accuracy you can request a lower limit on the accuracy, e.g. 0.1%
 * 		any calculation will have at least the accuracy of 0.1% then. Otherwise put 1 here.
 *
 * @param acceptance_rate here the a/r gets stored
 * @param accuracy here the accuracy gets stored
 *
 * @return iterations used
 */
unsigned int assess_acceptance_rate(mcmc * m, unsigned int param,
		double desired_acceptance_rate, double min_accuracy,
		double max_accuracy, double * acceptance_rate, double * accuracy) {
	unsigned int i = 0;
	unsigned int j;
	unsigned int n = 40;
	unsigned int accepts = 0;
	double stdev = 0;
	unsigned int maxdev = 0;
	double accept_rate;
	double required_accuracy = min_accuracy;
	char * acceptslog = NULL;
	gsl_vector * start_params = dup_vector(get_params(m));
	double start_prob = get_prob(m);

	reset_accept_rejects(m);

	while (1) {
		IFVERBOSE
			printf("calculating %d steps.\n", n);
		acceptslog = (char*) realloc(acceptslog, n * sizeof(char));
		assert(acceptslog != NULL);

		for (; i < n; i++) {
			accepts = get_params_accepts_for(m, param);
			/*restart_from_best(m);*/
			markov_chain_step_for(m, param);
			mcmc_check_best(m);
			set_prob(m, start_prob);
			set_params(m, dup_vector(start_params));
			if (accepts == get_params_accepts_for(m, param)) {
				/* had a reject -> set bit to 0 */
				clear_bit(acceptslog, i);
			} else {
				/* had a accept -> set bit to 1 */
				set_bit(acceptslog, i);
			}
		}
		accept_rate = accepts / (double) n;
		IFVERBOSE
			printf("accept rate: %f (%d/%d)\n", accept_rate, accepts, n);

		/* get max deviation */
		accepts = 0;
		stdev = 0;
		maxdev = 0 + 1;
		for (j = 0; j < n; j++) {
			if (get_bit(acceptslog, j) != 0) {
				accepts++;
			}
			stdev += pow(accepts - accept_rate * j, 2);
			if (abs(accepts - accept_rate * j) > maxdev) {
				maxdev = abs(accepts - accept_rate * j);
			}
		}
		stdev = sqrt(stdev / n) * 2;

		/*
		 * if we are way off, we don't need to be that accurate.
		 * if we are close, we want to be more accurate
		 * 30% could also be 10% to be more cautious
		 */
		required_accuracy = abs_double(accept_rate - desired_acceptance_rate)
				* ACCURACY_DEVIATION_FACTOR;
		if (required_accuracy < 0.005)
			required_accuracy = 0.005;
		if (required_accuracy < min_accuracy) {
			required_accuracy = min_accuracy;
		}
		if (required_accuracy > max_accuracy) {
			required_accuracy = max_accuracy;
		}

		/*
		 * we assume we have a deviation of maxdev at the end.
		 * how many values do we need to get below required_accuracy
		 */
		*acceptance_rate = accept_rate;
		*accuracy = maxdev / 1. / n;
		IFVERBOSE
			printf("accuracy wanted: %f, got: %f\n", required_accuracy,
					*accuracy);

		if (*accuracy <= required_accuracy) {
			break;
		}
		/*
		 * we need (maxdev / min_accuracy) values to achieve min_accuracy
		 */
		assert(maxdev / required_accuracy >= n);
		n = ((unsigned int) ((maxdev / 1. / required_accuracy) / 8) + 1) * 8;
	}
	return n;
}

void markov_chain_calibrate_linear_regression(mcmc * m,
		double desired_acceptance_rate, const double max_ar_deviation,
		const unsigned int iter_limit, gsl_matrix * all_stepwidths,
		gsl_matrix * all_acceptance_rates, gsl_matrix * all_accuracies) {
	unsigned int i;
	unsigned int j;
	unsigned int l;
	unsigned int n;
	double xy_bar;
	double xx_bar;
	double x_bar;
	double y_bar;
	double weight_sum;
	double k;
	double d;
	double sigma;
	unsigned int n_par = get_n_par(m);
	double calc_accuracy = max_ar_deviation * 2;
	double nextpoint;
	unsigned int iter = 0;
	double current_acceptance_rate;
	double accuracy;
	double min_stepwidth;
	double max_stepwidth;
	double weight;
	gsl_vector_int * enough_points_discovered = gsl_vector_int_alloc(n_par);
	FILE * progress_plot_file = fopen("calibration_progress.data", "w");
	gsl_vector_int_set_all(enough_points_discovered, 0);

	assert(progress_plot_file != NULL);
	for (j = 0; j < all_stepwidths->size2 && iter < iter_limit; j++) {
		for (i = 0; i < n_par; i++) {
			/*if (gsl_matrix_get(all_acceptance_rates, i, j) == 0)
			 gsl_matrix_set(all_acceptance_rates, i, j, -1);*/
			if (gsl_matrix_get(all_acceptance_rates, i, j) < 0)
				continue;

			fprintf(progress_plot_file, "%d\t%d\t%f\t%f\t%f\n", i + 1, iter,
					gsl_matrix_get(all_stepwidths, i, j), gsl_matrix_get(
							all_acceptance_rates, i, j), gsl_matrix_get(
							all_accuracies, i, j));
		}
	}
	debug("calibrating using linear regression");
	for (l = 0; l < all_stepwidths->size2 && iter < iter_limit; l++) {
		for (i = 0; i < n_par; i++) {
			if (gsl_vector_int_get(enough_points_discovered, i) != 0) {
				continue;
			}
			/* lets make a linear regression of it, to obtain k and d */
			xy_bar = 0;
			xx_bar = 0;
			x_bar = 0;
			y_bar = 0;
			weight_sum = 0;
			n = 0;
			min_stepwidth = 1;
			max_stepwidth = 0;
			for (j = 0; j < all_stepwidths->size2; j++) {
				if (gsl_matrix_get(all_acceptance_rates, i, j) < 0)
					continue;

				if (min_stepwidth > gsl_matrix_get(all_stepwidths, i, j))
					min_stepwidth = gsl_matrix_get(all_stepwidths, i, j);
				if (max_stepwidth < gsl_matrix_get(all_stepwidths, i, j))
					max_stepwidth = gsl_matrix_get(all_stepwidths, i, j);
				n++;

				weight = 1. / MIN(abs_double(gsl_matrix_get(
										all_acceptance_rates, i, j)
								- desired_acceptance_rate), 1);

				x_bar += gsl_matrix_get(all_stepwidths, i, j);
				y_bar += gsl_matrix_get(all_acceptance_rates, i, j) * weight;
				weight_sum += weight;
			}
			x_bar /= n;
			y_bar = y_bar / weight_sum;
			printf("  xbar = %f, ybar = %f, n=%d\n", x_bar, y_bar, n);

			for (j = 0; j < all_stepwidths->size2; j++) {
				if (gsl_matrix_get(all_acceptance_rates, i, j) < 0)
					continue;

				weight = 1. / MIN(abs_double(gsl_matrix_get(
										all_acceptance_rates, i, j)
								- desired_acceptance_rate), 1);

				printf("  %f | %f | weight = %f\n", gsl_matrix_get(
						all_stepwidths, i, j), gsl_matrix_get(
						all_acceptance_rates, i, j), weight);

				xy_bar += (gsl_matrix_get(all_stepwidths, i, j) - x_bar)
						* (gsl_matrix_get(all_acceptance_rates, i, j) - y_bar)
						* weight;
				xx_bar += pow(gsl_matrix_get(all_stepwidths, i, j) - x_bar, 2);
			}
			IFDEBUG
				printf("  %d points, xy = %f, xx = %f, weightsum = %f\n", n,
						xy_bar, xx_bar, weight_sum);

			k = xy_bar * n / weight_sum / xx_bar;
			d = y_bar - k * x_bar;

			/* also calculate sigma */
			sigma = 0;
			for (j = 0; j < all_stepwidths->size2; j++) {
				if (gsl_matrix_get(all_acceptance_rates, i, j) < 0)
					continue;
				weight = 1. / MIN(abs_double(gsl_matrix_get(
										all_acceptance_rates, i, j)
								- desired_acceptance_rate), 1);
				sigma += weight * pow(k * gsl_matrix_get(all_stepwidths, i, j)
						+ d - gsl_matrix_get(all_acceptance_rates, i, j), 2);
			}
			sigma = sqrt(sigma / weight_sum);

			IFDEBUG
				printf("%d: k = %g; d = %g; sigma = %g (from %d points)\n", i,
						k, d, sigma, n);

			/* with k and d, we can suggest a new point */
			nextpoint = (desired_acceptance_rate - d) / k;
			if (nextpoint > max_stepwidth + 0.2 * (max_stepwidth
					- min_stepwidth))
				nextpoint = max_stepwidth + 0.2 * (max_stepwidth
						- min_stepwidth);
			if (nextpoint < min_stepwidth - 0.2 * (max_stepwidth
					- min_stepwidth))
				nextpoint = min_stepwidth - 0.2 * (max_stepwidth
						- min_stepwidth);
			if (nextpoint > 1)
				nextpoint = 1;
			if (nextpoint < 0)
				nextpoint = min_stepwidth * 0.1;

			printf("%d: next stepwidth: %f\n", i, nextpoint);

			set_steps_for_normalized(m, nextpoint, i);

			/*
			 * now: sigma may be very high. So if we already have a point within
			 * the noise box:  max(calc_accuracy/3, sigma) / k
			 * it is not worth evaluating further. The noise is simply too high
			 * to get a better value around here. We would be test around forever.
			 */
			if (sigma < calc_accuracy / 3)
				sigma = calc_accuracy / 3;
			else
				sigma = sigma / 3;

			printf("%d: next stepwidth: %f +- %f\n", i, nextpoint,
					abs_double(sigma / k));

			if (n > 4) {
				for (j = 0; j < all_stepwidths->size2; j++) {
					if (gsl_matrix_get(all_acceptance_rates, i, j) < 0)
						continue;

					if (k < 0
							&& abs_double(gsl_matrix_get(all_stepwidths, i, j) - nextpoint)
									< abs_double(sigma / k)) {

						/* we already tested where we would go. */
						printf("%d: best stepwidth possible reached.\n", i);

						gsl_vector_int_set(enough_points_discovered, i, 1);
						break;
					}
				}
			}

			n = all_stepwidths->size2;
			for (j = 0; j < all_stepwidths->size2; j++) {
				if (gsl_matrix_get(all_acceptance_rates, i, j) < 0) {
					if (n == all_stepwidths->size2)
						n = j;
					continue;
				}
			}
			if (n == all_stepwidths->size2) {
				printf("%d: no space to store new points. should be enough.\n",
						i);
				gsl_vector_int_set(enough_points_discovered, i, 1);
			}

			if (gsl_vector_int_get(enough_points_discovered, i) == 0) {
				/* we evaluate the point, and add it at n */
				iter += assess_acceptance_rate(m, i, desired_acceptance_rate,
						sigma, sigma * 3, &current_acceptance_rate, &accuracy);

				gsl_matrix_set(all_stepwidths, i, n, nextpoint);
				gsl_matrix_set(all_acceptance_rates, i, n,
						current_acceptance_rate);
				gsl_matrix_set(all_accuracies, i, n, accuracy);
				IFDEBUG
					printf("%d: a/r: %f (+-%f)\n", i, current_acceptance_rate,
							accuracy);

				fprintf(progress_plot_file, "%d\t%d\t%f\t%f\t%f\t%f\n", i + 1,
						iter, gsl_matrix_get(all_stepwidths, i, n),
						gsl_matrix_get(all_acceptance_rates, i, n),
						gsl_matrix_get(all_accuracies, i, n), k);
				fflush(progress_plot_file);
			}
		}
	}
	fclose(progress_plot_file);
}

void markov_chain_calibrate_quadratic(mcmc * m,
		const double desired_acceptance_rate, const double max_ar_deviation,
		const unsigned int iter_limit, double mul, const double adjust_step) {

	unsigned int n_par = get_n_par(m);
	unsigned int i;
	unsigned int j;
	unsigned int k;
	unsigned int iter = 0;
	double current_acceptance_rate;
	double accuracy;
	double denominator;
	double nextpoint;
	double calc_accuracy = 0.01;
	double root;
	unsigned int n = 0;
	double a;
	double b;
	double c;
	gsl_matrix * stepwidths = gsl_matrix_alloc(n_par, 3);
	gsl_matrix * acceptance_rates = gsl_matrix_alloc(n_par, 3);
	gsl_vector_int * points_discovered = gsl_vector_int_alloc(n_par);
	gsl_matrix * all_stepwidths = gsl_matrix_alloc(n_par, 30);
	gsl_matrix * all_acceptance_rates = gsl_matrix_alloc(n_par, 30);
	gsl_matrix * all_accuracies = gsl_matrix_alloc(n_par, 30);

	gsl_matrix_set_all(acceptance_rates, 1);
	gsl_vector_int_set_all(points_discovered, 0);
	gsl_matrix_set_all(all_acceptance_rates, -1);

	mul = iter_limit + adjust_step; /* avoiding unused */

	nextpoint = 1;
	for (j = 0; nextpoint != -1; j++) {
		printf(" ==== ITERATION %d ==== \n", j);
		nextpoint = -1;
		for (i = 0; i < n_par; i++) {
			if (gsl_vector_int_get(points_discovered, i) == 0) {
				/* 1) discover first point */
				gsl_matrix_set(stepwidths, i, 0, get_steps_for_normalized(m, i));
				gsl_vector_int_set(points_discovered, i, 1);
				IFDEBUG
					printf("%d: point 0 (given): %f\n", i, gsl_matrix_get(
							stepwidths, i, 0));
			}
			if (gsl_vector_int_get(points_discovered, i) == 1) {

				iter += assess_acceptance_rate(m, i, desired_acceptance_rate,
						0, calc_accuracy, &current_acceptance_rate, &accuracy);

				gsl_matrix_set(acceptance_rates, i, 0, current_acceptance_rate);
				gsl_matrix_set(all_acceptance_rates, i, n,
						current_acceptance_rate);
				gsl_matrix_set(all_accuracies, i, n, accuracy);
				gsl_matrix_set(all_stepwidths, i, n, gsl_matrix_get(stepwidths,
						i, 0));
				n++;
				IFDEBUG
					printf("%d: a/r 0: %f (+-%f)\n", i, gsl_matrix_get(
							acceptance_rates, i, 0), accuracy);

				/*
				 * we have the first point.
				 * suggest a new one. We do not know much, so just make a suggestion
				 * We *do* know though, in which direction we should move:
				 */
				gsl_matrix_set(stepwidths, i, 1, gsl_matrix_get(stepwidths, i,
						0) * (1 + 5 * (current_acceptance_rate
						- desired_acceptance_rate)));

				if (gsl_matrix_get(stepwidths, i, 1) < 0)
					gsl_matrix_set(stepwidths, i, 1, 0.05 * gsl_matrix_get(
							stepwidths, i, 0));
				if (gsl_matrix_get(stepwidths, i, 1) > 1)
					gsl_matrix_set(stepwidths, i, 1, 1);
				IFDEBUG
					printf("%d: point 1 (guess): %f\n", i, gsl_matrix_get(
							stepwidths, i, 1));
				gsl_vector_int_set(points_discovered, i, 2);
			}
			if (gsl_vector_int_get(points_discovered, i) == 2) {
				/* 2) discover second point. */
				set_steps_for_normalized(m, gsl_matrix_get(stepwidths, i, 1), i);

				iter += assess_acceptance_rate(m, i, desired_acceptance_rate,
						0, calc_accuracy, &current_acceptance_rate, &accuracy);
				gsl_matrix_set(acceptance_rates, i, 1, current_acceptance_rate);
				gsl_matrix_set(all_acceptance_rates, i, n,
						current_acceptance_rate);
				gsl_matrix_set(all_accuracies, i, n, accuracy);
				gsl_matrix_set(all_stepwidths, i, n, gsl_matrix_get(stepwidths,
						i, 1));
				n++;
				IFDEBUG
					printf("%d: a/r 1: %f (+-%f)\n", i, gsl_matrix_get(
							acceptance_rates, i, 1), accuracy);

				/* 3) we have two points, lets make a line
				 * we want the third point to be at desired_acceptance_rate
				 * a * x^2 + b * x + c = acceptance_rate
				 */

				/* assume linear */
				a = 0;

				/* linear gradient */
				b = (gsl_matrix_get(acceptance_rates, i, 0) - gsl_matrix_get(
						acceptance_rates, i, 1)) / (gsl_matrix_get(stepwidths,
						i, 0) - gsl_matrix_get(stepwidths, i, 1));

				/* offset */
				c = gsl_matrix_get(acceptance_rates, i, 1) - b
						* gsl_matrix_get(stepwidths, i, 1);

				IFDEBUG
					printf(" a = %f; b = %f; c = %f (line %f*x+%f)\n", a, b, c,
							b, c);

				/* solve to desired_acceptance_rate */
				gsl_matrix_set(stepwidths, i, 2, (desired_acceptance_rate - c)
						/ b);
				if (gsl_matrix_get(stepwidths, i, 2) < 0)
					gsl_matrix_set(stepwidths, i, 2, 0.05 * gsl_matrix_get(
							stepwidths, i, 1));
				if (gsl_matrix_get(stepwidths, i, 2) > 1)
					gsl_matrix_set(stepwidths, i, 2, 1);
				IFDEBUG
					printf("%d: point 2 (line): %f\n", i, gsl_matrix_get(
							stepwidths, i, 2));
				gsl_vector_int_set(points_discovered, i, 3);
			}
			if (BETWEEN(gsl_vector_int_get(points_discovered, i), 3, 99)) {
				/* 3) discover third point. */
				set_steps_for_normalized(m, gsl_matrix_get(stepwidths, i, 2), i);

				iter += assess_acceptance_rate(m, i, desired_acceptance_rate,
						0, calc_accuracy, &current_acceptance_rate, &accuracy);
				gsl_matrix_set(acceptance_rates, i, 2, current_acceptance_rate);
				gsl_matrix_set(all_acceptance_rates, i, n,
						current_acceptance_rate);
				gsl_matrix_set(all_accuracies, i, n, accuracy);
				gsl_matrix_set(all_stepwidths, i, n, gsl_matrix_get(stepwidths,
						i, 2));
				n++;
				IFDEBUG
					printf("%d: a/r %d: %f (+-%f)\n", i, gsl_vector_int_get(
							points_discovered, i) - 1, gsl_matrix_get(
							acceptance_rates, i, 2), accuracy);

				if (j > 1) {
					/**
					 * break condition:
					 */
					printf("%d: a/r currently between [%f..%f] \n", i,
							min_column(acceptance_rates, 2), max_column(
									acceptance_rates, 2));
					if (max_column(acceptance_rates, 2)
							< desired_acceptance_rate + max_ar_deviation * 4
							&& min_column(acceptance_rates, 2)
									> desired_acceptance_rate
											- max_ar_deviation * 4) {
						printf("%d: a/r sufficient.\n", i);
						break;
					}
					if (j > 100)
						break;
				}

				/*
				 * We can now seek a quadratic polynome with the three points
				 * a * x^2 + b * x + c = acceptance_rate
				 *
				 * an algebra system tells us
				 *
				 * a=-(v_1*(x_3-x_2)-v_2*x_3+v_3*x_2+(v_2-v_3)*x_1)/(x_1*(x_3^2-x_2^2)-x_2*x_3^2+x_2^2*x_3+x_1^2*(x_2-x_3))
				 * b=(v_1*(x_3^2-x_2^2)-v_2*x_3^2+v_3*x_2^2+(v_2-v_3)*x_1^2)/(x_1*(x_3^2-x_2^2)-x_2*x_3^2+x_2^2*x_3+x_1^2*(x_2-x_3))
				 * c=(v_1*(x_2^2*x_3-x_2*x_3^2)+x_1*(v_2*x_3^2-v_3*x_2^2)+x_1^2*(v_3*x_2-v_2*x_3))/(x_1*(x_3^2-x_2^2)-x_2*x_3^2+x_2^2*x_3+x_1^2*(x_2-x_3))
				 *
				 * which boils down to the following ugly lines
				 * cross your fingers that you don't have to debug this
				 * (they are correct)
				 */
				denominator = gsl_matrix_get(stepwidths, i, 0) * (pow(
						gsl_matrix_get(stepwidths, i, 2), 2) - pow(
						gsl_matrix_get(stepwidths, i, 1), 2)) - gsl_matrix_get(
						stepwidths, i, 1) * pow(
						gsl_matrix_get(stepwidths, i, 2), 2) + pow(
						gsl_matrix_get(stepwidths, i, 1), 2) * gsl_matrix_get(
						stepwidths, i, 2) + pow(
						gsl_matrix_get(stepwidths, i, 0), 2) * (gsl_matrix_get(
						stepwidths, i, 1) - gsl_matrix_get(stepwidths, i, 2));

				a = -(gsl_matrix_get(acceptance_rates, i, 0) * (gsl_matrix_get(
						stepwidths, i, 2) - gsl_matrix_get(stepwidths, i, 1))
						- gsl_matrix_get(acceptance_rates, i, 1)
								* gsl_matrix_get(stepwidths, i, 2)
						+ gsl_matrix_get(acceptance_rates, i, 2)
								* gsl_matrix_get(stepwidths, i, 1)
						+ (gsl_matrix_get(acceptance_rates, i, 1)
								- gsl_matrix_get(acceptance_rates, i, 2))
								* gsl_matrix_get(stepwidths, i, 0))
						/ denominator;

				b = (gsl_matrix_get(acceptance_rates, i, 0) * (pow(
						gsl_matrix_get(stepwidths, i, 2), 2) - pow(
						gsl_matrix_get(stepwidths, i, 1), 2)) - gsl_matrix_get(
						acceptance_rates, i, 1) * pow(gsl_matrix_get(
						stepwidths, i, 2), 2) + gsl_matrix_get(
						acceptance_rates, i, 2) * pow(gsl_matrix_get(
						stepwidths, i, 1), 2) + (gsl_matrix_get(
						acceptance_rates, i, 1) - gsl_matrix_get(
						acceptance_rates, i, 2)) * pow(gsl_matrix_get(
						stepwidths, i, 0), 2)) / denominator;

				c = (gsl_matrix_get(acceptance_rates, i, 0) * (pow(
						gsl_matrix_get(stepwidths, i, 1), 2) * gsl_matrix_get(
						stepwidths, i, 2) - gsl_matrix_get(stepwidths, i, 1)
						* pow(gsl_matrix_get(stepwidths, i, 2), 2))
						+ gsl_matrix_get(stepwidths, i, 0) * (gsl_matrix_get(
								acceptance_rates, i, 1) * pow(gsl_matrix_get(
								stepwidths, i, 2), 2) - gsl_matrix_get(
								acceptance_rates, i, 2) * pow(gsl_matrix_get(
								stepwidths, i, 1), 2)) + pow(gsl_matrix_get(
						stepwidths, i, 0), 2) * (gsl_matrix_get(
						acceptance_rates, i, 2) * gsl_matrix_get(stepwidths, i,
						1) - gsl_matrix_get(acceptance_rates, i, 1)
						* gsl_matrix_get(stepwidths, i, 2))) / denominator;

				IFDEBUG
					printf(
							" a = %f; b = %f; c = %f (polynomial %f*x**2+%f*x+%f )\n",
							a, b, c, a, b, c);
				/*
				 * sigh.
				 * ok, assuming that is right, we can get the next point by solving
				 * a * x^2 + b * x + c = desired_acceptance_rate
				 * to x (quadratic equation)
				 *
				 * algebra system tells us:
				 * x=(+-sqrt(4*a*desired_acceptance_rate-4*a*c+b^2)-b)/(2*a)
				 */
				assert(4*a*(desired_acceptance_rate - c) + b*b >= 0);

				root = sqrt(4* a * (desired_acceptance_rate - c) + pow(b, 2));
				IFDEBUG
					printf(" polynomial solutions: %f, %f\n", (root - b) / (2
							* a), (-root - b) / (2* a ));
				/* select the solution we want */

				/* prefer solution between previous points */
				IFDEBUG
					printf(
							" should be in [%f..%f]; or in [0..1] or in [0..]\n",
							min_column(stepwidths, i),
							max_column(stepwidths, i));

				if (BETWEEN((root - b) / (2* a ), min_column(stepwidths, i),
						max_column(stepwidths, i))) {
					nextpoint = (root - b) / (2* a );
				} else if (BETWEEN((-root - b) / (2* a ), min_column(stepwidths, i),
						max_column(stepwidths, i))) {
					nextpoint = (-root - b) / (2* a );
				} else {
					/* prefer solution in [0..1] */
					/*
					 * TODO: add bounds: if we have prior knowledge of a definite
					 * high or low bound, jump somewhere between.
					 *
					 * e.g. we know a point that has a acceptance rate 5% too
					 * high, so don't jump beyond that
					 * or we know a point that has a acceptance rate 5% too low
					 * so don't jump below that.
					 */
					if (BETWEEN((root - b) / (2* a ), -0.02, 1.1)) {
						nextpoint = (root - b) / (2* a );
					} else if (BETWEEN((-root - b) / (2* a ), -0.02, 1.1)) {
						nextpoint = (-root - b) / (2* a );
					} else {
						/* ok, this may suggest a high stepwidth (>1). */
						if ((root - b) / (2* a ) >= 0) {
							nextpoint = (root - b) / (2* a );
						} else if ((-root - b) / (2* a ) >= 0) {
							nextpoint = (-root - b) / (2* a );
						} else {
							fprintf(stderr,
									"ERROR: found no solution where stepwidth > 0\n");
							exit(1);
						}
						IFDEBUG
							printf("WARNING, landed at high stepwidth\n");
					}
				}
				IFDEBUG
					printf("%d: point %d (polynomial): %f\n", i,
							gsl_vector_int_get(points_discovered, i), nextpoint);
				if (gsl_vector_int_get(points_discovered, i) > 4
						&& BETWEEN((nextpoint - min_column(stepwidths, i))/
								(max_column(stepwidths, i) - min_column(stepwidths, i)),
								-0.1, 1.1)) {
					IFDEBUG
						printf(
								"  we are extrapolating. dropping out of quad mode\n");
					gsl_vector_int_set(points_discovered, i,
							gsl_vector_int_get(points_discovered, i) + 100);
				}
				if (gsl_vector_int_get(points_discovered, i) > 2) {
					if ((b + 2 * a * MIN(nextpoint, min_column(stepwidths, i))
							> 0) || (b + 2 * a
							* MAX(nextpoint, max_column(stepwidths, i)) > 0)) {
						IFDEBUG
							printf(
									"polynome not monotone. dropping out of quad mode\n");
						gsl_vector_int_set(points_discovered, i,
								gsl_vector_int_get(points_discovered, i) + 100);
					}
				}

				/* ok, we got our next point. but to repeat, we want it in this
				 * order: (x1, x2, x3). currently we have (x0, x1, x2)
				 */
				gsl_matrix_set(stepwidths, i, 0, gsl_matrix_get(stepwidths, i,
						1));
				gsl_matrix_set(stepwidths, i, 1, gsl_matrix_get(stepwidths, i,
						2));
				gsl_matrix_set(stepwidths, i, 2, nextpoint);
				if (gsl_matrix_get(stepwidths, i, 2) < 0)
					gsl_matrix_set(stepwidths, i, 2, 0.1 * gsl_matrix_get(
							stepwidths, i, 1));
				if (gsl_matrix_get(stepwidths, i, 2) > 1)
					gsl_matrix_set(stepwidths, i, 2, 1);
				nextpoint = gsl_matrix_get(stepwidths, i, 2);
				gsl_matrix_set(acceptance_rates, i, 0, gsl_matrix_get(
						acceptance_rates, i, 1));
				gsl_matrix_set(acceptance_rates, i, 1, gsl_matrix_get(
						acceptance_rates, i, 2));
				gsl_matrix_set(acceptance_rates, i, 2, 0 /* unknown yet */);
				gsl_vector_int_set(points_discovered, i, gsl_vector_int_get(
						points_discovered, i) + 1);

				if (gsl_vector_int_get(points_discovered, i) > 100) {
					/*
					 * dropping out now.
					 *
					 * we have obtained n data points, but some may be too far away
					 * for a linear relationship.
					 * We have a last point on our hands, that has not been
					 * discovered yet. [at nextpoint]
					 *
					 * We also have a estimate of the gradient! it is
					 *   k = 2*a*stepwidth + b
					 *
					 * at nextpoint, the tangent is
					 *   (2*a*nextpoint + b)*stepwidth + c
					 *
					 * Lets look for the furthest datapoints that are within
					 * the box
					 *   stepwidth = nextpoint +- desired_accuracy / k * factor
					 */
					printf(
							"deleting points more than %f outside\n",
							abs_double(1. / (2 * a * nextpoint + b) * calc_accuracy * 10));
					j = 0;
					for (k = 0; gsl_matrix_get(all_acceptance_rates, i, k) >= 0; k++) {
						if (abs_double(gsl_matrix_get(all_stepwidths, i, k) - nextpoint)
								> abs_double(1. / (2 * a * nextpoint + b) * calc_accuracy
										* 10)) {
							IFDEBUG
								printf("  deleting stepwidth %f - a/r %f\n",
										gsl_matrix_get(all_stepwidths, i, k),
										gsl_matrix_get(all_acceptance_rates, i,
												k));
							gsl_matrix_set(all_acceptance_rates, i, k, -1);
						} else {
							j++;
						}
					}
					if (j < 3) {
						/* add border points of this rectangle */
						set_steps_for_normalized(m, nextpoint, i);

						iter += assess_acceptance_rate(m, i,
								desired_acceptance_rate, 0, calc_accuracy,
								&current_acceptance_rate, &accuracy);
						gsl_matrix_set(acceptance_rates, i, 1,
								current_acceptance_rate);
						gsl_matrix_set(all_acceptance_rates, i, n,
								current_acceptance_rate);
						gsl_matrix_set(all_accuracies, i, n, accuracy);
						gsl_matrix_set(all_stepwidths, i, n,
								get_steps_for_normalized(m, i));
						n++;
						IFDEBUG
							printf("%d: a/r %d: %f (+-%f)\n", i, n,
									gsl_matrix_get(acceptance_rates, i, 1),
									accuracy);
						j++;
					}
					if (j < 3) {
						/* add border points of this rectangle */
						set_steps_for_normalized(m, nextpoint - calc_accuracy
								* 3 / (2 * a * nextpoint + b), i);

						iter += assess_acceptance_rate(m, i,
								desired_acceptance_rate, 0, calc_accuracy,
								&current_acceptance_rate, &accuracy);
						gsl_matrix_set(acceptance_rates, i, 1,
								current_acceptance_rate);
						gsl_matrix_set(all_acceptance_rates, i, n,
								current_acceptance_rate);
						gsl_matrix_set(all_accuracies, i, n, accuracy);
						gsl_matrix_set(all_stepwidths, i, n,
								get_steps_for_normalized(m, i));
						n++;
						IFDEBUG
							printf("%d: a/r %d: %f (+-%f)\n", i, n,
									gsl_matrix_get(acceptance_rates, i, 1),
									accuracy);

						set_steps_for_normalized(m, nextpoint + calc_accuracy
								* 3 / (2 * a * nextpoint + b), i);

						iter += assess_acceptance_rate(m, i,
								desired_acceptance_rate, 0, calc_accuracy,
								&current_acceptance_rate, &accuracy);
						gsl_matrix_set(acceptance_rates, i, 1,
								current_acceptance_rate);
						gsl_matrix_set(all_acceptance_rates, i, n,
								current_acceptance_rate);
						gsl_matrix_set(all_accuracies, i, n, accuracy);
						gsl_matrix_set(all_stepwidths, i, n,
								get_steps_for_normalized(m, i));
						n++;
						IFDEBUG
							printf("%d: a/r %d: %f (+-%f)\n", i, n,
									gsl_matrix_get(acceptance_rates, i, 1),
									accuracy);

					}
				}
			}
		}
	}
	markov_chain_calibrate_linear_regression(m, desired_acceptance_rate,
			max_ar_deviation, iter_limit - iter, all_stepwidths,
			all_acceptance_rates, all_accuracies);
	gsl_matrix_free(stepwidths);
	gsl_matrix_free(acceptance_rates);
	gsl_vector_int_free(points_discovered);
	gsl_matrix_free(all_stepwidths);
	gsl_matrix_free(all_acceptance_rates);
	gsl_matrix_free(all_accuracies);
}

#ifndef MAX_ACCURACY_IMPROVEMENT
#define MAX_ACCURACY_IMPROVEMENT 2.8
#endif

#ifndef SCALE_LIN_WORST
#define SCALE_LIN_WORST 5
#endif
#ifndef SCALE_MIN
#define SCALE_MIN 0.4
#endif

void markov_chain_calibrate_alt(mcmc * m, double desired_acceptance_rate,
		const double max_ar_deviation, const unsigned int iter_limit,
		double mul, const double adjust_step) {

	unsigned int i;
	unsigned int j;
	double current_acceptance_rate;
	double accuracy;
	unsigned int n_par = get_n_par(m);
	double scale = 1.2;
	double movedirection;
	double move;
	double max_deviation;
	double worst_accuracy = 0;
	double worst_accuracy_previous = 0;
	double best_worst_accuracy = 1;
	unsigned int iter = 0;
	FILE * progress_plot_file = fopen("calibration_progress.data", "w");
	gsl_vector * accuracies = gsl_vector_alloc(n_par);
	gsl_vector_set_all(accuracies, 0);

	mul = iter_limit + adjust_step; /* avoiding unused */

	/*
	 * we use assess_acceptance_rate for a n-dim point in the stepwidth-space
	 * we want to minimize | acceptance_rate - rat_limit|.
	 * if acceptance_rate >> rat_limit, we should increase the stepwidth
	 * if acceptance_rate << rat_limit, we should decrease the stepwidth
	 */

	while (1) {
		max_deviation = 0;
		/* assess current situation */
		for (j = 0; j < 1; j++) {
			printf("calculating for up to %f accuracy\n",
					worst_accuracy_previous / MAX_ACCURACY_IMPROVEMENT);
			worst_accuracy = 0;
			for (i = 0; i < n_par; i++) {
				if (gsl_vector_get(accuracies, i) < 0.1
						* worst_accuracy_previous) {
					continue;
				}

				/*
				 * the idea is to reuse the accuracy of the worst parameter in the
				 * previous round. So in this round, we only want an accuracy
				 * improvement of 3 times that. Why? This is a n-dimensional
				 * minimization problem and the parameters are not independent.
				 * So it wouldn't help us if we determined the value of one
				 * dimension extremely accurately, but be far away in another dim.
				 */
				iter += assess_acceptance_rate(m, i, desired_acceptance_rate,
						worst_accuracy_previous / MAX_ACCURACY_IMPROVEMENT,
						1 /* no restriction */, &current_acceptance_rate,
						&accuracy);
				printf("%d: a/r: %f (+-%f); desired: %f; steps: %f\n", i,
						current_acceptance_rate, accuracy,
						desired_acceptance_rate, get_steps_for_normalized(m, i));

				fprintf(progress_plot_file, "%d\t%d\t%f\t%f\t%f\n", i + 1,
						iter, get_steps_for_normalized(m, i),
						current_acceptance_rate, accuracy);

				/* keep track of worst performer */
				/*if (worst_accuracy < accuracy) {*/
				worst_accuracy += accuracy;
				/*}*/
				gsl_vector_set(accuracies, i, accuracy);

				movedirection = current_acceptance_rate
						- desired_acceptance_rate;
				/*
				 * reduce scale if we have already settled down before.
				 * We don't want the steps jumping around wildly
				 */
				scale = best_worst_accuracy * SCALE_LIN_WORST + SCALE_MIN;
				assert(scale > 0);
				move = movedirection * scale;
				if (move < -1)
					move = -0.9;
				if (max_deviation < abs_double(movedirection)) {
					max_deviation = abs_double(movedirection);
				}
				/* 10% too high => increase steps by 10% */
				/* 10% too low  => decrease steps by 10% */

				set_steps_for(m, get_steps_for(m, i) * (1 + move), i);

				printf("%d: new steps: %f\n", i, get_steps_for_normalized(m, i));
			}
			if (iter > iter_limit * n_par) {
				fprintf(stderr, "calibration failed: iteration limit reached\n");
				exit(1);
			}
			worst_accuracy_previous = worst_accuracy / n_par;
			if (worst_accuracy_previous < best_worst_accuracy) {
				best_worst_accuracy = worst_accuracy;
			}
		}
		printf("max deviation: %f; ", max_deviation);
		dump_v("current values", get_params(m));

		if (max_deviation < max_ar_deviation && worst_accuracy
				< max_ar_deviation * 2) {
			printf("small deviation: %f; quitting\n", max_deviation);
			break;
		}
	}
	fclose(progress_plot_file);

}

void markov_chain_calibrate_orig(mcmc * m, double rat_limit,
		const double max_rat_deviation, const unsigned int iter_limit,
		double mul, const double adjust_step) {
	/* we aim a acceptance rate between 20 and 30% */
	unsigned int i;

	int reached_perfection = 0;
	gsl_vector * accept_rate = NULL;
	double delta_reject_accept_t;

	unsigned long iter = 0;
	unsigned long subiter;
	int nchecks_without_rescaling = 0;
	int rescaled;

	/* we want a acceptance rate for each parameter */
	rat_limit = pow(rat_limit, 1.0 / get_n_par(m));

	gsl_vector_scale(m->params_step, adjust_step);
	debug("Calibrating step widths ...");
	reset_accept_rejects(m);

	while (1) {
		for (i = 0; i < get_n_par(m); i++) {
			markov_chain_step_for(m, i);
			mcmc_check_best(m);
		}
		iter++;
		if (iter % ITER_READJUST == 0) {
			accept_rate = get_accept_rate(m);

			dump_ul("------------------------------------------------ iteration", iter);
			dump_v("params", get_params(m));
			dump_v("acceptance rate: ", accept_rate);
			dump_v("steps", get_steps(m));

			rescaled = 0;
			for (i = 0; i < get_n_par(m); i++) {
				IFDEBUG
					printf(
							"\t\tneeded acceptance rate: <%f, >%f; got %f for %i",
							rat_limit + 0.05, rat_limit - 0.05, gsl_vector_get(
									accept_rate, i), i);
				if (gsl_vector_get(accept_rate, i) > rat_limit + 0.05) {
					set_steps_for(m, get_steps_for(m, i) / mul, i);
					IFDEBUG
						printf("\t scaling up   ^");
					if (rescaled == 0)
						rescaled = -1;
					if (get_steps_for_normalized(m, i) > 1) {
						printf(
								"\nWARNING: step width of %s is quite big! %d times the param space\n",
								get_params_descr(m)[i], (int) (gsl_vector_get(
										m->params_step, i) / (gsl_vector_get(
										m->params_max, i) - gsl_vector_get(
										m->params_min, i))));
						printf(
								"\nWARNING: This can mean the parameter is independent.\n");
						printf("\n SETTING PARAMETER STEP TO PARAMETER RANGE\n");
						set_steps_for_normalized(m, 1, i);
						if (rescaled == -1)
							rescaled = 0;
					}
					if (get_steps_for_normalized(m, i) > 10000) {
						fprintf(
								stderr,
								"calibration failed: step width of %s became too large.\n",
								get_params_descr(m)[i]);
						exit(1);
					}
					if (rescaled == -1)
						rescaled = 1;
				}
				if (gsl_vector_get(accept_rate, i) < rat_limit - 0.05) {
					set_steps_for(m, get_steps_for(m, i) * mul, i);
					IFDEBUG
						printf("\t scaling down v");
					if (get_steps_for_normalized(m, i) < 10E-10) {
						printf(
								"\nWARNING: step width of %s is quite small! %e times the param space\n",
								get_params_descr(m)[i],
								get_steps_for_normalized(m, i));
					}
					rescaled = 1;
				}
				IFDEBUG
					printf("\n");
				assert(gsl_vector_min(get_steps(m)) > 0);
			}
			if (rescaled == 0)
				nchecks_without_rescaling++;
			else dump_v("steps", m->params_step);
			restart_from_best(m);
			reset_accept_rejects(m);
			for (subiter = 0; subiter < ITER_READJUST; subiter++) {
				markov_chain_step(m);
				mcmc_check_best(m);
			}
			gsl_vector_free(accept_rate);
			accept_rate = get_accept_rate(m);
			dump_v("New overall accept rate after reset", accept_rate);
			gsl_vector_free(accept_rate);
			delta_reject_accept_t = get_accept_rate_global(m)
					- TARGET_ACCEPTANCE_RATE;
			dump_d("Compared to desired rate", delta_reject_accept_t);
			if (abs_double(delta_reject_accept_t) < max_rat_deviation) {
				reached_perfection = 1;
				debug("calibration reached the desired acceptance rate");
				printf("\n %d steps without rescaling \n",
						nchecks_without_rescaling);
			} else {
				reached_perfection = 0;
				if (delta_reject_accept_t < 0) {
					rat_limit /= 0.99;
				} else {
					rat_limit *= 0.99;
				}
			}
			if (nchecks_without_rescaling >= NO_RESCALING_LIMIT
					&& reached_perfection == 1 && rescaled == 0) {
				debug("quitting calibration because we did not need to rescale for several times");
				break;
			}
			if (iter > iter_limit) {
				fprintf(stderr,
						"calibration failed: limit of %d iterations reached.",
						iter_limit);
				exit(1);
			}
		}
	}
	reset_accept_rejects(m);
	debug("calibration of markov-chain done.");
}

void markov_chain_calibrate(mcmc * m, const unsigned int burn_in_iterations,
		double desired_acceptance_rate, const double max_ar_deviation,
		const unsigned int iter_limit, double mul, const double adjust_step) {

	burn_in(m, burn_in_iterations);

	dump_d("desired acceptance rate per parameter", desired_acceptance_rate);

#ifdef CALIBRATE_QUADRATIC
	markov_chain_calibrate_quadratic
#else
#ifdef CALIBRATE_ALTERNATE
	markov_chain_calibrate_alt
#else
	markov_chain_calibrate_orig
#endif
#endif
	(m, desired_acceptance_rate, max_ar_deviation, iter_limit, mul, adjust_step);
}

void do_step_for(mcmc * m, const unsigned int i) {
	const double step = gsl_vector_get(m->params_step, i);
	const double old_value = gsl_vector_get(m->params, i);
	double new_value = old_value + get_next_gauss_random(m, step);
	const double max = gsl_vector_get(m->params_max, i);
	const double min = gsl_vector_get(m->params_min, i);
	/* dump_d("Jumping from", old_value); */

#if CIRCULAR_PARAMS == 0
	while (new_value > max || new_value < min) {
		new_value = old_value + get_next_gauss_random(m, step);
		IFVERBOSE
			printf("Value borders reached; looking for new starting point"
				" for %d \n", i);
	}
#else
	unsigned int j = 0;
	/** circular parameters **/
	unsigned int parameters[] = {CIRCULAR_PARAMS, 0};

	if (new_value > max || new_value < min) {
		while (1) {
			if (parameters[j] == 0) {
				/* non-circular parameter */
				do {
					new_value = old_value + get_next_gauss_random(m, step);
				}while (new_value > max || new_value < min);
				break;
			}
			if (parameters[j] == i + 1) {
				/* circular parameter */
				new_value = min + mod_double(new_value - min, max - min);
				break;
			}
			j++;
		}
	}
#endif
	assert(new_value <= max);
	assert(new_value >= min);
	/* dump_d("To", new_value); */
	set_params_for(m, new_value, i);
}

static void do_step(mcmc * m) {
	unsigned int i;
	for (i = 0; i < get_n_par(m); i++) {
		do_step_for(m, i);
	}
}

/**
 * @returns 1 if accept, 0 if rejecting
 */
static int check_accept(mcmc * m, const double prob_old) {
	double prob_new = get_prob(m);
	double prob_still_accept;

	/* shortcut */
	if (prob_new == prob_old) {
		return 1;
	}
	IFVERBOSE
		dump_v("suggesting parameter", get_params(m));

	if (prob_new > prob_old) {
		IFVERBOSE
			dump_d("accepting improvement of", prob_new - prob_old);
		return 1;
	} else {
		prob_still_accept = get_next_alog_urandom(m);
		if (prob_still_accept < (prob_new - prob_old)) {
			IFVERBOSE {
				dump_d("accepting probability", prob_still_accept);
				dump_d("accepting worsening of", prob_new - prob_old);
			}
			return 1;
		} else {
			IFVERBOSE
				dump_d("rejecting worsening of", prob_new - prob_old);
			return 0;
		}
	}
}

static void revert(mcmc * m, const double prob_old) {
	set_prob(m, prob_old);
}

void markov_chain_step_for(mcmc * m, const unsigned int index) {
	double prob_old = get_prob(m);
	double old_value = gsl_vector_get(m->params, index);

	mcmc_check(m);
	do_step_for(m, index);

	calc_model_for(m, index, old_value);

	if (check_accept(m, prob_old) == 1) {
		inc_params_accepts_for(m, index);
	} else {
		revert(m, prob_old);
		set_params_for(m, old_value, index);
		inc_params_rejects_for(m, index);
	}
}

#ifndef MINIMAL_STEPWIDTH
#define MINIMAL_STEPWIDTH 0.0000001
#endif
#ifndef MAXIMAL_STEPWIDTH
#define MAXIMAL_STEPWIDTH 1000000
#endif

void rmw_adapt_stepwidth(mcmc * m, const double prob_old) {
	unsigned int i;
	double step;
	double min;
	double scale;
	double max;
	double alpha = exp(get_prob(m) - prob_old);
	if (alpha > 1)
		alpha = 1;
	for (i = 0; i < get_n_par(m); i++) {
		scale = (gsl_vector_get(m->params_max, i) - gsl_vector_get(
				m->params_min, i));
		min = MINIMAL_STEPWIDTH * scale;
		max = MAXIMAL_STEPWIDTH * scale;

		step = gsl_vector_get(get_steps(m), i);
		step += get_next_uniform_random(m) / sqrt(m->n_iter) * (alpha
				- TARGET_ACCEPTANCE_RATE) * scale;
		;
		if (step < min)
			step = min;
		if (step > max)
			step = max;
		gsl_vector_set(get_steps(m), i, step);
	}
}

void markov_chain_step(mcmc * m) {
	double prob_old = get_prob(m);
	gsl_vector * old_values = dup_vector(m->params);

	mcmc_check(m);
	do_step(m);

	calc_model(m, old_values);

	if (check_accept(m, prob_old) == 1) {
		inc_params_accepts(m);
		gsl_vector_free(old_values);
	} else {
		revert(m, prob_old);
		set_params(m, old_values);
		inc_params_rejects(m);
	}
}
