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

#include "utils.h"
#include <ctype.h>

FILE * openfile(const char * filename) {
	FILE * input = fopen(filename, "r");
	if (input == NULL) {
		fprintf(stderr, "error opening file %s\n", filename);
		perror("file could not be opened");
		exit(1);
	}
	return input;
}

unsigned int countlines(const char * filename) {
	int nlines = 0;
	int c;
	int r;
	FILE * input = openfile(filename);
	while (1) {
		c = fgetc(input);
		if (c == '\n') {
			nlines++;
		}
		if (c == EOF)
			break;
	}
	r = fclose(input);
	assert (r == 0);
	return nlines;
}

unsigned int get_column_count(const char * filename) {
	FILE * input;
	char buf[10000];
	unsigned int count = 0;
	int i = 0;
	int at_whitespace = 1;

	input = openfile(filename);
	if (fgets(buf, 10000, input) == NULL) {
		fprintf(stderr, "error: file %s is empty!", filename);
		exit(1);
	}
	while (buf[i] != 0) {
		if (isspace(buf[i])) {
			if (at_whitespace == 0)
				count++;
			at_whitespace = 1;
		} else {
			at_whitespace = 0;
		}
		i++;
	}
	if (at_whitespace == 0)
		count++;
	return count;
}
