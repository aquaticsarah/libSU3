/* libSU3: Tests for isoscalar factor calculations */

#include <stdio.h>

#include "SU3.h"
#include "test.h"

/* Check the values for 3 x 3bar -> 1 + 8.
    Note that a lot of combinations are ruled out by hypercharge,
    so we don't check those here. */
TEST(3x3bar)
{
    isoarray* isf;

    /* Check the singlet part of the decomposition */
    isf = isoscalars(0, 0, 1, 0, 0, 1);
                    //   n  k  l k1 l1 k2 l2  num denom
    TEST_EQ_SQRAT((*isf)(0, 0, 0, 0, 0, 1, 1), -1, 3, "(0,0): (0,0) x (1,1)");
    TEST_EQ_SQRAT((*isf)(0, 0, 0, 1, 0, 1, 0),  2, 3, "(0,0): (1,0) x (1,0)");

    delete isf;

    /* Check the octet part of the decomposition */
    isf = isoscalars(1, 1, 1, 0, 0, 1);
                    //   n  k  l k1 l1 k2 l2  num denom
    TEST_EQ_SQRAT((*isf)(0, 1, 0, 0, 0, 1, 0), 1, 1, "(1,0): (0,0) x (1,0)");
    TEST_EQ_SQRAT((*isf)(0, 1, 1, 0, 0, 1, 1), 2, 3, "(1,1): (0,0) x (1,1)");
    TEST_EQ_SQRAT((*isf)(0, 1, 1, 1, 0, 1, 0), 1, 3, "(1,1): (1,0) x (1,0)");
    TEST_EQ_SQRAT((*isf)(0, 2, 0, 0, 0, 1, 1), 0, 1, "(2,0): (0,0) x (1,1)");
    TEST_EQ_SQRAT((*isf)(0, 2, 0, 1, 0, 1, 0), 1, 1, "(2,0): (1,0) x (1,0)");
    TEST_EQ_SQRAT((*isf)(0, 2, 1, 1, 0, 1, 1), 1, 1, "(2,1): (1,0) x (1,1)");

    delete isf;
}

void print_isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (d == 0) return;

    isoarray* isf = isoscalars(p, q, p1, q1, p2, q2);

    long n, k, l, k1, l1, k2, l2;
    char buf[64];

    printf("Isoscalar coefficients for (%ld,%ld): (%ld,%ld) x (%ld,%ld)\n",
        p, q, p1, q1, p2, q2);
    for (n = 0; n < d; ++n)
    {
        printf("Degenerate rep %ld:\n", n);

        for (k = q; k <= p+q; ++k)
            for (l = 0; l <= q; ++l)
                for (k1 = q1; k1 <= p1+q1; ++k1)
                    for (l1 = 0; l1 <= q1; ++l1)
                        for (k2 = q2; k2 <= p2+q2; ++k2)
                        {
                            l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3 - (k1 + l1 + k2 - k - l);
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

TEST(isoscalars)
{
#if 0
    long p1, q1, p2, q2, p, q;
    isoarray* isf;
    /* Check that our recurrence relations work (ie, don't crash)
        for a range of different sets of reps */
    for (p1 = 0; p1 < 5; ++p1)
        for (q1 = 0; q1 < 5; ++q1)
            for (p2 = 0; p2 < 5; ++p2)
                for (q2 = 0; q2 < 5; ++q2)
                    for (p = 0; p < 5; ++p)
                        for (q = 0; q < 5; ++q)
                        {
                            isf = isoscalars(p1,q1,p2,q2,p,q);
                            delete isf;
                        }
    fprintf(stderr, "Isoscalar calculations okay\n");
#else
    print_isoscalars(2,2,1,1,1,1);
    print_isoscalars(1,1,2,2,1,1);
    print_isoscalars(1,1,1,1,2,2);

    fprintf(stderr, "Isoscalar calculations okay\n");
#endif
}
