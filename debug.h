#ifndef DEBUGHELPER
#define DEBUGHELPER

#include <stdio.h>
#ifdef DEBUG
#define IFDEBUG if(1)
#else
#define IFDEBUG if(0)
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
#define dump_v(str, v)    IFDEBUG { printf("\tDEBUG[%s]: %s: ", AT, str); dump_vector(v); fflush(NULL); }
#define dump_m(str, m)    IFDEBUG { printf("\tDEBUG[%s]: %s\n", AT, str); dump(m); fflush(NULL); }

#define wait() IFWAIT { printf("if it is ok to continue, press return. "); fflush(NULL); fgetc(stdin); }

void dump(mcmc * m);
void dump_vector(gsl_vector * v);
void require(const int returncode);

#endif

