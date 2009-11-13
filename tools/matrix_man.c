/*
    APEMoST - Automated Parameter Estimation and Model Selection Toolkit
    Copyright (C) 2009  Johannes Buchner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
 
/**
 * matrix_man allows manipulation of double values tables, such as
 * subtracting, adding, multiplying columns and rows in a automatable way.
 */

#include <stdio.h>
#include <stdlib.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>

#include <string.h>
#include <ctype.h>

#define SEP '\t'
#define OUTPUTFORMAT "%e\t"

char * outputformat = "%e\t";

#ifdef DEBUG
#define IFDEBUG if(1)
#else
#define IFDEBUG if(0)
#endif

void usage(char * progname) {
	printf("%s: SYNOPSIS: <command>\n"
		"\n"
		"matrix_man allows manipulation of double values tables\n"
		"It reads csv or tab-seperated values (double) from stdin,\n"
		"performs the commands and prints the result to stdout.\n"
		"\ncurrently supported commands\n", progname);
	printf(""
		"\trm <i>           remove column i\n"
		"\tnew <i>          add a new column at i and fill it with zeroes\n"
		"\tswap <i> <o>  swap the value in column i+o with i\n"
		"\tinvm <i>         set i to i^(-1)\n"
		"\tinva <i>         set i to -i\n"
		"\tdiff <i>         set i to the difference to the previous line\n");
	printf(""
		"\tadd <i> <j>      add the value in column j to i. \n"
		"\t                 i and j can only specify one column each\n"
		"\tadd_rel <i> <o>  add the value in column i+o to i\n"
		"\taddc <i> <v>     add the value v to i\n"
		"\tmul <i> <j>      multiply the value in column j onto the value in i\n"
		"\t                 i and j can only specify one column each\n"
		"\tmulc <i> <v>     multiply i with the value v\n"
		"\tmul_rel <i> <o>  multiply the value in column i+o onto i\n");
	printf("\n"
		"For any column number a column selection can be given:\n"
		"\ti%%2=0,i>3,i<2,i=2,i/4\n"
		"The character ',' acts as a boolean 'and'. '/' means unequal. \n"
		"above, o is always a offset number\n"
		"\n");
	exit(1);

}

/**
 * @returns -1 on error
 */
int qualifies_parse(int i, char * qualification_in) {
	int j = 0;
	int k = 0;
	int last;
	int len = strlen(qualification_in);
	char * qualification;
	int tokens[5] = { -1, -1, -1, -1, -1 };
	int result = -1;

	qualification = (char*) calloc(len + 1, sizeof(char));
	if (qualification == NULL) {
		perror("could not allocate buffer");
		exit(-1);
	}
	strcpy(qualification, qualification_in);

	while (qualification[k] != ',' && k < len) {
		if (qualification[k] == '%' || qualification[k] == '<'
				|| qualification[k] == '>' || qualification[k] == '='
				|| qualification[k] == '/') {
			tokens[j] = k;
			j++;
		}
		k++;
	}
	if (k == len)
		last = 1;
	else {
		last = 0;
		qualification[k] = 0;
	}
	IFDEBUG
		printf("DEBUG: qual: %s, i: %d\n", qualification, i);
	IFDEBUG
		printf("DEBUG: qual: len: %d, mylen: %d, ntokens: %d\n", len, k, j);
	if (j == 1 && qualification[0] == 'i') {
		IFDEBUG
			printf("DEBUG: qual: relation\n");
		if (qualification[tokens[0]] == '<') {
			if (i + 1 < atoi(&qualification[tokens[0] + 1])) {
				result = 1;
			} else {
				result = 0;
			}
		} else if (qualification[tokens[0]] == '>') {
			if (i + 1 > atoi(&qualification[tokens[0] + 1])) {
				result = 1;
			} else {
				result = 0;
			}
		} else if (qualification[tokens[0]] == '=') {
			if (i + 1 == atoi(&qualification[tokens[0] + 1])) {
				result = 1;
			} else {
				result = 0;
			}
		} else if (qualification[tokens[0]] == '/') {
			if (i + 1 != atoi(&qualification[tokens[0] + 1])) {
				result = 1;
			} else {
				result = 0;
			}
		}
	} else if (j == 2 && qualification[0] == 'i' && qualification[tokens[0]]
			== '%' && qualification[tokens[1]] == '=') {
		IFDEBUG
			printf("DEBUG: qual: modulo\n");
		qualification[tokens[0]] = 0;
		qualification[tokens[1]] = 0;
		if ((i + 1) % atoi(&qualification[tokens[0] + 1]) == atoi(
				&qualification[tokens[1] + 1])) {
			result = 1;
		} else {
			result = 0;
		}
	}
	free(qualification);
	IFDEBUG
		printf("DEBUG: qual: result: %d\n", result);
	if (result == -1)
		return result;
	if (result == 0 || last == 1) {
		return result;
	} else {
		return qualifies_parse(i, &qualification_in[k + 1]);
	}
}

int qualifies(int i, char * qualification) {
	int result = qualifies_parse(i, qualification);
	if (result == -1) {
		fprintf(stderr, "error on qualification string\n");
		exit(2);
	} else {
		return result;
	}
}

void command_rm(int count, gsl_vector * current, char * qualification) {
	int i;
	for (i = 0; i < count; i++) {
		if (!qualifies(i, qualification))
			printf(outputformat, gsl_vector_get(current, i));
	}
	printf("\n");
	fflush(stdout);
}
void command_inva(int count, gsl_vector * current, char * qualification) {
	int i;
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification))
			printf(outputformat, -gsl_vector_get(current, i));
		else
			printf(outputformat, gsl_vector_get(current, i));
	}
	printf("\n");
	fflush(stdout);
}
void command_invm(int count, gsl_vector * current, char * qualification) {
	int i;
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification))
			printf(outputformat, 1.0 / gsl_vector_get(current, i));
		else
			printf(outputformat, gsl_vector_get(current, i));
	}
	printf("\n");
	fflush(stdout);
}
void command_diff(int count, gsl_vector * prev, gsl_vector * current,
		char * qualification) {
	int i;
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification))
			printf(outputformat, gsl_vector_get(current, i) - gsl_vector_get(
					prev, i));
		else
			printf(outputformat, gsl_vector_get(current, i));
	}
	printf("\n");
	fflush(stdout);
}

void command_new(int count, gsl_vector * current, char * qualification) {
	int offset = 0;
	int i;
	for (i = 0; i < count; i++) {
		if (qualifies(i + offset, qualification)) {
			printf(outputformat, 0.0);
			offset++;
			i--;
		} else {
			printf(outputformat, gsl_vector_get(current, i));
		}
	}
	if (qualifies(i, qualification))
		printf(outputformat, 0.0);

	printf("\n");
	fflush(stdout);
}
void command_add(int count, gsl_vector * current, char * qualification_source,
		char * qualification_target) {
	int source = -1;
	int target = -1;
	int i;
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification_source) && source == -1) {
			source = i;
		}
		if (qualifies(i, qualification_target) && target == -1) {
			target = i;
		}
	}
	if (source == -1) {
		fprintf(stderr, "no source match\n");
		exit(-1);
	}
	if (target == -1) {
		fprintf(stderr, "no target match\n");
		exit(-1);
	}
	for (i = 0; i < count; i++) {
		if (i == target)
			printf(outputformat, gsl_vector_get(current, i) + gsl_vector_get(
					current, source));
		else
			printf(outputformat, gsl_vector_get(current, i));
	}
	printf("\n");
	fflush(stdout);
}
void command_add_rel(int count, gsl_vector * current,
		char * qualification_source, char * offset_str) {
	int i;
	int offset = atoi(offset_str);
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification_source)) {
			printf(outputformat, gsl_vector_get(current, i) + gsl_vector_get(
					current, i + offset));
		} else {
			printf(outputformat, gsl_vector_get(current, i));
		}
	}
	printf("\n");
	fflush(stdout);
}

void command_addc(int count, gsl_vector * current, char * qualification_source,
		char * value_str) {
	int i;
	double value;
	sscanf(value_str, "%lf", &value);
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification_source)) {
			printf(outputformat, gsl_vector_get(current, i) + value);
		} else {
			printf(outputformat, gsl_vector_get(current, i));
		}
	}
	printf("\n");
	fflush(stdout);
}
void command_mul(int count, gsl_vector * current, char * qualification_source,
		char * qualification_target) {
	int source = -1;
	int target = -1;
	int i;
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification_source) && source == -1) {
			source = i;
		}
		if (qualifies(i, qualification_target) && target == -1) {
			target = i;
		}
	}
	if (source == -1) {
		fprintf(stderr, "no source match\n");
		exit(-1);
	}
	if (target == -1) {
		fprintf(stderr, "no target match\n");
		exit(-1);
	}
	for (i = 0; i < count; i++) {
		if (i == target)
			printf(outputformat, gsl_vector_get(current, i) * gsl_vector_get(
					current, source));
		else
			printf(outputformat, gsl_vector_get(current, i));
	}
	printf("\n");
	fflush(stdout);
}
void command_mul_rel(int count, gsl_vector * current,
		char * qualification_source, char * offset_str) {
	int i;
	int offset = atoi(offset_str);
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification_source)) {
			printf(outputformat, gsl_vector_get(current, i) * gsl_vector_get(
					current, i + offset));
		} else {
			printf(outputformat, gsl_vector_get(current, i));
		}
	}
	printf("\n");
	fflush(stdout);
}
void command_mulc(int count, gsl_vector * current, char * qualification_source,
		char * value_str) {
	int i;
	double value;
	sscanf(value_str, "%lf", &value);
	for (i = 0; i < count; i++) {
		if (qualifies(i, qualification_source)) {
			printf(outputformat, gsl_vector_get(current, i) * value);
		} else {
			printf(outputformat, gsl_vector_get(current, i));
		}
	}
	printf("\n");
	fflush(stdout);
}

void row_man(int argc, char ** argv) {
	gsl_vector * prev = NULL;
	gsl_vector * current = NULL;
	FILE * input = stdin;
	char buf[10000];
	int count;
	int i;
	int col;
	int at_whitespace;
	int line = 0;
	char * command = argv[1];
	double v;

	while (!feof(input)) {
		line++;
		current = gsl_vector_alloc(300);
		if (fgets(buf, 10000, input) == NULL || feof(input))
			break;
		count = 0;
		i = 0;
		at_whitespace = 1;
		while (buf[i] != 0) {
			if (isspace(buf[i])) {
				if (at_whitespace == 0)
					count++;
				at_whitespace = 1;
			} else {
				if (at_whitespace == 1) {
					col = sscanf(&buf[i], "%lf", &v);
					IFDEBUG
						printf("field %d: %f\n", count, v);
					gsl_vector_set(current, count, v);

					if (col != 1) {
						fprintf(stderr,
								"field could not be read: %d, line %d\n", i
										+ count, line);
						exit(1);
					}
				}
				at_whitespace = 0;
			}
			i++;
		}
		if (at_whitespace == 0)
			count++;

		IFDEBUG
			printf("command: %s; nargs = %d\n", command, argc);

		if (0 == strcmp(command, "rm") && argc == 3)
			command_rm(count, current, argv[2]);
		else if (0 == strcmp(command, "new") && argc == 3)
			command_new(count, current, argv[2]);
		else if (0 == strcmp(command, "add") && argc == 4)
			command_add(count, current, argv[2], argv[3]);
		else if (0 == strcmp(command, "mul") && argc == 4)
			command_mul(count, current, argv[2], argv[3]);
		else if (0 == strcmp(command, "add_rel") && argc == 4)
			command_add_rel(count, current, argv[2], argv[3]);
		else if (0 == strcmp(command, "mul_rel") && argc == 4)
			command_mul_rel(count, current, argv[2], argv[3]);
		else if (0 == strcmp(command, "addc") && argc == 4)
			command_addc(count, current, argv[2], argv[3]);
		else if (0 == strcmp(command, "mulc") && argc == 4)
			command_mulc(count, current, argv[2], argv[3]);
		else if (0 == strcmp(command, "invm") && argc == 3)
			command_invm(count, current, argv[2]);
		else if (0 == strcmp(command, "inva") && argc == 3)
			command_inva(count, current, argv[2]);
		else if (0 == strcmp(command, "diff") && argc == 3) {
			if (prev != NULL) {
				command_diff(count, prev, current, argv[2]);
				gsl_vector_free(prev);
			}
			prev = current;
			current = NULL;
		} else {
			fprintf(stderr, "Unknown command\n");
			usage(argv[0]);
		}
		if (current != NULL)
			gsl_vector_free(current);
		fflush(stdout);
	}

}

int main(int argc, char ** argv) {
	if (argc < 3) {
		fprintf(stderr, "Wrong number of arguments\n");
		usage(argv[0]);
	} else {
		if (strcmp(argv[1], "-f") == 0) {
			outputformat = argv[2];

			argc -= 2;
			argv += 2;
		}

		row_man(argc, argv);
	}

	return 0;
}

