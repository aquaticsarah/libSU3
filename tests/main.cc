/* libSU3: Test driver */

#include <stdio.h>
#include <limits.h>

#include "SU3.h"
#include "test.h"

int tests_run;
int tests_passed;

/* TODO: Autogenerate a list of tests */
void test_gcd();
void test_sqrat();
void test_sqrat_arith();

int main()
{
    /* TODO: Autogenerate code like this for each top-level test */
    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting gcd...\n");
    test_gcd();
    fprintf(stderr, "gcd: %d/%d tests passed\n", tests_passed, tests_run);

    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting sqrat...\n");
    test_sqrat();
    fprintf(stderr, "sqrat: %d/%d tests passed\n", tests_passed, tests_run);

    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting sqrat-arithmetic...\n");
    test_sqrat_arith();
    fprintf(stderr, "sqrat-arithmetic: %d/%d tests passed\n", tests_passed, tests_run);

    return 0;
}
