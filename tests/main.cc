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

void print_isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (d == 0) return;

    isoarray* isf = isoscalars(p, q, p1, q1, p2, q2);

    long n, k, l, k1, l1, k2, l2;
    char buf[64];

    for (n = 0; n < d; ++n)
    {
        printf("Degenerate rep %ld:\n", n);

        /* For now, only look at the state of highest weight */
        k = p+q; l = 0;

        for (k1 = q1; k1 <= p1+q1; ++k1)
            for (l1 = 0; l1 <= q1; ++l1)
                for (k2 = q2; k2 <= p2+q2; ++k2)
                {
                    l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3 - (k1 + l1 + k2);
                    if ((l2 < 0) || (l2 > q2)) continue;

                    sqrat val = (*isf)(n, k, l, k1, l1, k2, l2);
                    val.tostring(buf, 64);
                    printf("    (%ld,%ld) : (%ld,%ld) x (%ld,%ld) = %s\n",
                        k, l, k1, l1, k2, l2, buf);
                }
        printf("\n");
    }

    delete isf;
}

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
    /* Decompose 3 x 3bar -> 1 + 8 and print coefficients */
    print_isoscalars(0,0,1,0,0,1);
    print_isoscalars(1,1,1,0,0,1);

    /* Also try 8 x 8 -> 8 */
    print_isoscalars(1,1,1,1,1,1);

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
