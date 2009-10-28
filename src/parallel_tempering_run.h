#ifndef _PARALLEL_TEMPERING_RUN
#define _PARALLEL_TEMPERING_RUN

extern int run;
extern int dumpflag;

int get_duration();
void register_signal_handlers();
long unsigned int get_ticks_per_second();

#endif
