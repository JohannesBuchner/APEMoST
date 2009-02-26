
/* we do testing of the inner works */
#include <ctype.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>

extern int (*tests_registration[])(void);

int count_tests()
{
	int n = -1;
	while(tests_registration[++n] != NULL);
	return n;
}

/* 
 * when testnr in 1..n, execute the specific test.
 * output is conforming to TAP 
 */
int test(int testnr){
	int n = count_tests();
	int r;
	if(testnr <= 0 || testnr > n){
		fprintf(stderr, "Test number out of range (1..%d)\n", n);
		return -1;
	}
	r = tests_registration[testnr-1]();
	if(r == 0)
		printf("ok %d\n", testnr);
	else
		printf("not ok %d - returned %d\n", testnr, r);
	return r;
}
int test_all(){
	int i = 0;
	int n = count_tests();
	int count = 0;
	printf("1..%d\n", n);
	while ( i++ < n){
		count += test(i);
	}
	if(count != 0){
		printf("  %d of %d tests unsuccessful\n", count, n);
	}else{
		printf("  all %d tests successful\n", n);
	}
	return count;
}

void usage(char * progname){
	printf( "SYNOPSIS: \n\n"
		"%s               \trun all tests\n"
		"%s [<testnumber>]\trun a specific test by number\n"
		"\n"
		"Tests available: 1..%d\n\n"
		"", progname, progname, count_tests());
	
}

int main(int argc, char ** argv) {
	int t;
	if(argc==1) {
		return test_all();
	} else if(argc == 2 && isdigit(argv[1][0])) {
		errno = 0;
		t = strtol(argv[1], (char **) NULL, 10);
		if(errno != 0) {
			usage(argv[0]);
		} else {
			test(t);
		}
	} else {
		usage(argv[0]);
	}
	return 0;
}

