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
MCMC_SOURCES=src/*.c
TEST_SOURCES=tests/*.c
APPS=apps/*.c
MAIN=src/generic_main.c

## help: this clutter
help:
	@grep -E '^## ' Makefile|grep -v ':'|sed 's,^## ,,g'
	@echo This Makefile has the targets:
	@grep -E '^## [.a-z]{2,}:' Makefile|sed 's,^## *,\t,g' |sed 's,: ,\t,g'

## all: 
all: tests.exe

tests.exe: $(MCMC_SOURCES) $(TEST_SOURCES)
	$(CC) -I src $(CFLAGS) $(LDFLAGS) $^ -o $@

%.exe: apps/%.c $(MCMC_SOURCES) $(MAIN) 
	$(CC) -I src $(CFLAGS) $(LDFLAGS) $^ -o $@

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
