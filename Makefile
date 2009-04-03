## Makefile for mcmc
## 
## You can call it like this if you have segfaults
##  $ LD_PRELOAD=libduma.so.0.0.0 CCFLAGS="-DDEBUG -DSEGV" make tests
## 
## You might want to compile your gsl with -DHAVE_INLINE 
## 

CFLAGS=-O2 -fopenmp -Wall -Werror -Wextra -g -ansi -pedantic ${CCFLAGS}
LDFLAGS=-lgsl -lgslcblas -lm
CC=gcc
MCMC_SOURCES=mcmc.c mcmc_gettersetter.c mcmc_parser.c mcmc_calculate.c \
	mcmc_markov_chain.c mcmc_dump.c 
PARALLEL_TEMPERING_SOURCES=parallel_tempering.c tempering_interaction.c
TEST_SOURCES=run-tests.c tests.c
HELPER_SOURCES=gsl_helper.c 
DEBUG_SOURCE=debug.c
HEADERS=*.h

## help: this clutter
help:
	@grep -E '^## ' Makefile|grep -v ':'|sed 's,^## ,,g'
	@echo This Makefile has the targets:
	@grep -E '^## [.a-z]{2,}:' Makefile|sed 's,^## *,\t,g' |sed 's,: ,\t,g'
	$(CC) --version

## all: 
all: tests.exe simplesin.exe

%.o : %.c
	$(CC) $(CFLAGS) -c $^

tests.exe: $(HEADERS) tests.o run-tests.o $(MCMC_SOURCES) $(DEBUG_SOURCE) $(HELPER_SOURCES) 
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

simplesin.exe: $(HEADERS) simplesin.c $(PARALLEL_TEMPERING_SOURCES) $(MCMC_SOURCES) $(DEBUG_SOURCE) $(HELPER_SOURCES) 
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@
simplesin5.exe: $(HEADERS) simplesin5.c $(PARALLEL_TEMPERING_SOURCES) $(MCMC_SOURCES) $(DEBUG_SOURCE) $(HELPER_SOURCES) 
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

## tests: run the tests
tests: tests.exe
	./tests.exe ${TESTNR}

## clean: 
clean:
	rm -f *.exe *.o

## doc: build the documentation
doc: doxygen.config
	doxygen doxygen.config


.PHONY: clean tests run all help doc
