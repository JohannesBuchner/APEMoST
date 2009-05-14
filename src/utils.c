
#include "utils.h"

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
