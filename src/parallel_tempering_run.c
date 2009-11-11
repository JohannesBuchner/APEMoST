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

#include "mcmc.h"
#include "parallel_tempering_run.h"

#include <signal.h>
#include <time.h>

int run = 1;
int dumpflag;

void ctrl_c_handler(int signalnr) {
	printf("\nreceived Ctrl-C (%d). Stopping ... (please be patient)\n\n",
			signalnr);
	run = 0;
}
void sigusr_handler(int signalnr) {
	printf("\nreceived SIGUSR (%d). Will dump at next opportunity.\n\n",
			signalnr);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);
	dumpflag = 1;
}

#include <time.h>

int get_duration() {
	static clock_t stored = 0;
	clock_t new = stored;
	stored = clock();
	return stored - new;
}

void register_signal_handlers() {
	signal(SIGINT, ctrl_c_handler);
	signal(SIGUSR2, sigusr_handler);
	signal(SIGUSR1, sigusr_handler);
}

long unsigned int get_ticks_per_second() {
	return CLOCKS_PER_SEC;
}

