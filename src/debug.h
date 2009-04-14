#ifndef DEBUGHELPER
#define DEBUGHELPER

#include <stdio.h>

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Turns on debug output.
 *
 * One could think that turning this off would improve performance,
 * test have shown the additional call time is not significant.
 */
#define DEBUG
#endif

#ifdef DEBUG
#define IFDEBUG if(1)
#else
#define IFDEBUG if(0)
#endif

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * Turns on debug output for allocation and free()
 * to track segfaults.
 *
 * You can also use libduma or valgrind:
 *
 * <h3>DUMA</h3>
 * LD_PRELOAD=libduma.so.0.0.0 yourprogram
 * <h3>valgrind</h3>
 * valgrind --tool=memcheck --leak-check=full yourprogram
 */
#define SEGV
#endif

#ifdef SEGV
#define IFSEGV if(1)
#else
#define IFSEGV if(0)
#endif

#ifdef WAIT
#define IFWAIT if(1)
#else
#define IFWAIT if(0)
#endif

#ifdef __NEVER_SET_FOR_DOCUMENTATION_ONLY
/**
 * additionally to DEBUG, print more.
 */
#define VERBOSE
#endif

#ifdef VERBOSE
#define IFVERBOSE if(1)
#else
#define IFVERBOSE if(0)
#endif

#include "memory.h"
#include "mcmc.h"

/* http://www.decompile.com/cpp/faq/file_and_line_error_string.htm */
#define STRINGIFY(x) #x
#define TOSTRING(x) STRINGIFY(x)
#define AT __FILE__ ":" TOSTRING(__LINE__)

#define debug(str)        IFDEBUG { printf("\tDEBUG[%s]: %s\n", AT, str); fflush(NULL); }
#define dump_i(str, var)  IFDEBUG { printf("\tDEBUG[%s]: %s: %i\n", AT, str, var); fflush(NULL); }
#define dump_ui(str, var) IFDEBUG { printf("\tDEBUG[%s]: %s: %u\n", AT, str, var); fflush(NULL); }
#define dump_d(str, var)  IFDEBUG { printf("\tDEBUG[%s]: %s: %f\n", AT, str, var); fflush(NULL); }
#define dump_l(str, var)  IFDEBUG { printf("\tDEBUG[%s]: %s: %l\n", AT, str, var); fflush(NULL); }
#define dump_ul(str, var) IFDEBUG { printf("\tDEBUG[%s]: %s: %lu\n", AT, str, var); fflush(NULL); }
#define dump_size(str, var) IFDEBUG { printf("\tDEBUG[%s]: %s: %lu\n", AT, str, (unsigned long)var); fflush(NULL); }
#define dump_i_s(str, index, var) \
                          IFDEBUG { printf("\tDEBUG[%s]: %s[%i]: %s\n", AT, str, index, var); fflush(NULL); }
#define dump_s(str, var)  IFDEBUG { printf("\tDEBUG[%s]: %s: %s\n", AT, str, var); fflush(NULL); }
#define dump_p(str, var)  IFDEBUG { printf("\tDEBUG[%s]: %s: %p\n", AT, str, var); fflush(NULL); }
#define dump_v(str, v)    IFDEBUG { printf("\tDEBUG[%s]: %s: ", AT, str); dump_vectorln(v); fflush(NULL); }
#define dump_m(str, m)    IFDEBUG { printf("\tDEBUG[%s]: %s\n", AT, str); dump(m); fflush(NULL); }

/**
 * Wait for keyboard input
 */
#define wait() IFWAIT { printf("if it is ok to continue, press return. "); fflush(NULL); fgetc(stdin); }

void dump(mcmc * m);
void dump_vector(gsl_vector * v);
void dump_vectorln(gsl_vector * v);

#define require(x) (x) /* there is a gsl handler so we don't need that */
#ifndef require
void require(const int returncode);
#endif

#endif
