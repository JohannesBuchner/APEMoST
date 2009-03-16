#ifndef DEBUGHELPER
#define DEBUGHELPER

#include <stdio.h>
#include "mcmc.h"

#ifdef DEBUG
#define IFDEBUG if(1)
#else
#define IFDEBUG if(0)
#endif

#define debug(str)        IFDEBUG printf("DEBUG: %s\n", str);
#define dump_i(str, var)  IFDEBUG printf("DEBUG: %s: %i\n", str, var);
#define dump_ui(str, var) IFDEBUG printf("DEBUG: %s: %u\n", str, var);
#define dump_d(str, var)  IFDEBUG printf("DEBUG: %s: %f\n", str, var);
#define dump_l(str, var)  IFDEBUG printf("DEBUG: %s: %l\n", str, var);
#define dump_ul(str, var) IFDEBUG printf("DEBUG: %s: %lu\n", str, var);
#define dump_i_s(str, index, var) IFDEBUG printf("DEBUG: %s[%i]: %s\n", str, index, var);
#define dump_s(str, var)  IFDEBUG printf("DEBUG: %s: %s\n", str, var);
#define dump_p(str, var)  IFDEBUG printf("DEBUG: %s: %p\n", str, var);
#define dump_m(str, m)    IFDEBUG dump(m);

void dump(mcmc * m); 

#endif

