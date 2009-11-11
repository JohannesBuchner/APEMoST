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
 * Enabling the garbage collector
 */

#ifndef MEMORY_H_
#define MEMORY_H_

#include "debug.h"

#define FREEMSG(x) IFSEGV dump_p("about to free", (void*)x);

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * The garbage collector (boehmgc) is enabled by default, 
 * but is not required. This disables it.
 * 
 * You may also want to remove the -lgc flag.
 **/
#define WITHOUT_GARBAGE_COLLECTOR
#endif

#ifdef WITHOUT_GARBAGE_COLLECTOR

#define mem_malloc(x) malloc(x)
#define mem_calloc(n,x) calloc(n, x)
#define mem_realloc(p,x) realloc(p,x)
#define mem_free(x) { FREEMSG(x); free((void*)x); }

#else

#include <gc/gc.h>

#define mem_malloc(x) GC_malloc(x)
#define mem_calloc(n,x) GC_malloc((n)*(x))
#define mem_realloc(p,x) GC_realloc((p),(x))
#define mem_free(x) { FREEMSG(x); (x) = NULL; }

#endif

#ifdef NOFREE
#define gsl_vector_free(v)
#endif

#endif /* MEMORY_H_ */
