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

/* A family of tests designed to exercise the code which calculates the ISFs
    for one rep by using the symmetry relations, rather than directly.
    TODO: These tests should pass iff the calculation function doesn't
    throw an exception.
*/
TEST(symmetries)
{
    isoarray* isf;

    /* Can be calculated directly */
    isf = isoscalars(2, 3, 1, 2, 1, 1);
    DO_TEST(isf != NULL, "Empty irrep");
    delete isf;

    /* Needs to use the 1<->3bar symmetry */
    isf = isoscalars(1, 1, 3, 2, 1, 2);
    DO_TEST(isf != NULL, "Empty irrep");
    delete isf;

    /* Needs to use the 2<->3bar symmetry */
    isf = isoscalars(1, 1, 1, 2, 3, 2);
    DO_TEST(isf != NULL, "Empty irrep");
    delete isf;
}

/* Test only that there are no crashes, but for a large number
   of representations */
TEST(isoscalars_no_crashes)
{
    isoarray* isf;

    long p, q, p1, q1, p2, q2;
    for (p = 0; p < 5; ++p)
        for (q = 0; q < 5; ++q)
            for (p1 = 0; p1 < 5; ++p1)
                for (q1 = 0; q1 < 5; ++q1)
                    for (p2 = 0; p2 < 5; ++p2)
                        for (q2 = 0; q2 < 5; ++q2)
                        {
                            isf = isoscalars(p, q, p1, q1, p2, q2);
                            delete isf;
                        }

    DO_TEST(1, "Shouldn't happen");
}
