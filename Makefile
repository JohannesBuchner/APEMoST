# Makefile for mcmc
# 

CFLAGS=-O2 -Wall -Werror -Wextra -g -ansi -pedantic ${CCFLAGS}
LDFLAGS=-lgsl -lgslcblas -lm
# optional (ignore the warning)
DUMA=libduma.so.0.0
# -DHAVE_INLINE 
CC=gcc
MCMC_SOURCES=mcmc.c mcmc_gettersetter.c mcmc_parser.c mcmc_dump.c 
TEST_SOURCES=run-tests.c tests.c
HELPER_SOURCES=gsl_helper.c 
DEBUG_SOURCE=debug.c
HEADERS=*.h

# help: this clutter
help:
	@echo This Makefile has the targets:
	@grep -E '^# [.a-z]{2,}:' Makefile|sed 's,^# *,\t,g' |sed 's,: ,\t,g'

# all: 
all: tests.exe

%.o : %.c
	$(CC) $(CFLAGS) -c $^

tests.exe: $(HEADERS) tests.o run-tests.o $(MCMC_SOURCES) $(DEBUG_SOURCE) $(HELPER_SOURCES) 
	$(CC) $(CFLAGS) $(LDFLAGS) $^ -o $@

# tests: run the tests
tests: tests.exe
	LD_PRELOAD=${DUMA} ./tests.exe

# clean: 
clean:
	rm -f *.exe *.o

# doc: build the documenttion
doc: doxygen.config
	doxygen doxygen.config


.PHONY: clean tests run all help doc
