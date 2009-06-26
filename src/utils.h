#ifndef MCMC_UTILS
#define MCMC_UTILS

#include <stdio.h>
#include <stdlib.h>

#ifdef NOASSERT
#define assert(cond)
#else
#include <assert.h>
#endif

/**
 * open the file or die
 */
FILE * openfile(const char * filename);

/**
 * a wc -l on the file.
 */
unsigned int countlines(const char * filename);

/**
 * how many columns does this file have? (looks at first line only)
 */
unsigned int get_column_count(const char * filename);
#endif
