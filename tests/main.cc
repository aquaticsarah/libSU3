/* libSU3: Test driver */

#include <stdio.h>

#include "SU3.h"
#include "test.h"

int tests_run;
int tests_passed;

TEST(tmp) void do_test()
{
	DO_TEST(test() == 1, "Expected test() to return 1");
}

int main()
{
	/* TODO: Autogenerate code like this for each top-level test */
	tests_run = 0;
	tests_passed = 0;
	fprintf(stderr, "\nTesting tmp...\n");
	do_test();
	fprintf(stderr, "tmp: %d/%d tests passed\n", tests_passed, tests_run);

	return 0;
}
