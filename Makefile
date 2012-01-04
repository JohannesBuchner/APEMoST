## Makefile for APEMoST
## 
## You can call it like this if you have segfaults
##  $ LD_PRELOAD=libduma.so.0.0.0 CCFLAGS="-DDEBUG -DSEGV" make tests
## 
## You might want to compile your gsl with -DHAVE_INLINE 
## 

CFLAGS := -I src -O3 -std=c99 -fopenmp -Wall -Werror -Wextra -ansi -pedantic ${CCFLAGS}
LDFLAGS := -lgsl -lgslcblas -lm

ifdef WITH_GARBAGE_COLLECTOR
LDFLAGS := ${LDFLAGS} -lgc
else
CFLAGS := ${CFLAGS} -DWITHOUT_GARBAGE_COLLECTOR
endif

CC := gcc
COMMON_SOURCES := src/gsl_helper.c src/histogram.c src/debug.c src/utils.c
COMMON := $(COMMON_SOURCES:.c=.o)
MCMC_SOURCES := $(wildcard src/mcmc*.c)
MCMC := $(MCMC_SOURCES:.c=.o)
MARKOV_CHAIN_SOURCES := $(wildcard src/markov_chain*.c)
MARKOV_CHAIN := $(MARKOV_CHAIN_SOURCES:.c=.o)
PARALLEL_TEMPERING_SOURCES := $(wildcard src/parallel_tempering*.c) src/analyse.c
PARALLEL_TEMPERING := $(PARALLEL_TEMPERING_SOURCES:.c=.o)
TEST_SOURCES := $(wildcard tests/*.c)
TEST := $(TEST:.c=.o)
MAIN := apps/generic_main.o
BENCH_MAIN := apps/benchmark_main.o
EVAL_MAIN := apps/eval_main.o

AR := ar r
LINKLIB := ${CC} -shared ${LDFLAGS}
LIBDIR := ./

## help: this clutter
help:
	@grep -E '^## ' Makefile|grep -v ':'|sed 's,^## ,,g'
	@echo This Makefile has the targets:
	@grep -E '^## [.a-z]{2,}:' Makefile|sed 's,^## *,\t,g' |sed 's,: ,\t,g'

## all: 
all: tests.exe tools simplesin.exe benchmark_simplesin.exe eval_simplesin.exe libapemost.so

tools: histogram_tool.exe random_tool.exe ndim_histogram_tool.exe sum_tool.exe matrix_man.exe peaks.exe

libapemost.a: $(MCMC) $(COMMON) $(MARKOV_CHAIN) $(PARALLEL_TEMPERING)
	$(AR) $@ $^

libapemost.so: libapemost.a
	$(LINKLIB) -o $(LIBS) $@ $^

tests.exe: $(TEST_SOURCES) libapemost.a
	$(CC) $(CFLAGS) $(LDFLAGS) -L. -lapemost $^ -o $@

%.exe: apps/%.c $(MAIN) libapemost.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

benchmark_%.exe: apps/%.c $(BENCH_MAIN) libapemost.a
	$(CC) -pg $(CFLAGS) $(LDFLAGS) $^ -o $@

eval_%.exe: apps/%.c $(EVAL_MAIN) libapemost.a
	$(CC) -pg $(CFLAGS) $(LDFLAGS) $^ -o $@

%.exe: tools/%.c libapemost.a
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

## tests: run the tests
tests: tests.exe
	./tests.exe ${TESTNR}

## clean: 
clean:
	rm -f *.exe *.o src/*.o *.a *.so

## doc: build the documentation
doc: doxygen.config
	doxygen doxygen.config
	make -C doc/

.PHONY: clean tests run all help doc
