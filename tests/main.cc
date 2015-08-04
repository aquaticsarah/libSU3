/* libSU3: Test driver */

#include <stdio.h>
#include <limits.h>

#include "SU3.h"
#include "test.h"

int tests_run;
int tests_passed;

/* TODO: Autogenerate a list of tests */
void test_sqrat();
void test_sqrat_arithmetic();
void test_dimension();
void test_degeneracy();
void test_rep_sizes();

int main()
{
#if 0
    long p1, q1, p2, q2, p, q;
    /* Check that our recurrence relations work (ie, don't crash)
        for a range of different sets of reps */
    for (p1 = 0; p1 < 5; ++p1)
        for (q1 = 0; q1 < 5; ++q1)
            for (p2 = 0; p2 < 5; ++p2)
                for (q2 = 0; q2 < 5; ++q2)
                    for (p = 0; p < 5; ++p)
                        for (q = 0; q < 5; ++q)
                            isoscalars(p1,q1,p2,q2,p,q);
    fprintf(stderr, "Isoscalar calculations okay\n");
#else
    isoscalars(1,1,1,1,1,1);
    fprintf(stderr, "Isoscalar calculations okay\n");
#endif

    /* TODO: Autogenerate code like this for each top-level test */
    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting sqrat...\n");
    test_sqrat();
    fprintf(stderr, "sqrat: %d/%d tests passed\n", tests_passed, tests_run);

    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting sqrat_arithmetic...\n");
    test_sqrat_arithmetic();
    fprintf(stderr, "sqrat_arithmetic: %d/%d tests passed\n", tests_passed, tests_run);

    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting dimension...\n");
    test_dimension();
    fprintf(stderr, "dimension: %d/%d tests passed\n", tests_passed, tests_run);

    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting degeneracy...\n");
    test_degeneracy();
    fprintf(stderr, "degeneracy: %d/%d tests passed\n", tests_passed, tests_run);

    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\nTesting rep_sizes...\n");
    test_rep_sizes();
    fprintf(stderr, "rep_sizes: %d/%d tests passed\n", tests_passed, tests_run);

    return 0;
}
