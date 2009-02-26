# Makefile for mcmc
# 

CFLAGS=-O3 -Wall -Werror -Wextra -g -ansi -pedantic
LDFLAGS=-lgsl -lgslcblas -lm
# -DHAVE_INLINE 
CC=gcc

# help: this clutter
help:
	@echo targets:
	@grep -E '^# [.a-z]{2,}:' Makefile|sed 's,^# *,\t,g' |sed 's,: ,\t,g'

all: main.exe

run: main.exe
	./main.exe

tests.exe: Makefile mcmc.h run-tests.c tests.c idl-functions.c
	$(CC) $(CFLAGS) $(LDFLAGS) run-tests.c tests.c idl-functions.c -o $@

# tests: run the tests
tests: tests.exe
	./tests.exe
	
# clean: 
clean:
	rm -f *.exe
	
%.exe: %.c Makefile mcmc.h
	$(CC) $(CFLAGS) $(LDFLAGS) $< -o $@

.PHONY: clean tests run all help
