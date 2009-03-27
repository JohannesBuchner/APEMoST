/**
 * Enabling the garbage collector
 */

#ifndef MEMORY_H_
#define MEMORY_H_

#include "debug.h"

#define FREEMSG(x) IFSEGV dump_p("about to free", (void*)x);

#define WITH_GARBAGE_COLLECTOR

#ifdef WITH_GARBAGE_COLLECTOR

#include <gc.h>

#define mem_malloc(x) GC_malloc(x)
#define mem_calloc(n,x) GC_malloc((n)*(x))
#define mem_realloc(p,x) GC_realloc((p),(x))
#define mem_free(x) { FREEMSG(x); (x) = NULL; }

#else

#define mem_malloc(x) malloc(x)
#define mem_calloc(n,x) calloc(n, x)
#define mem_realloc(p,x) realloc(p,x)
#define mem_free(x) { FREEMSG(x); free(x); (x) = NULL; }

#endif

/*#define gsl_vector_free(v)*/

#endif /* MEMORY_H_ */
