/* libSU3: Tests for irrep information */

#include <stdlib.h>
#include <time.h>

#include "SU3.h"
#include "test.h"

TEST(dimension)
{
    long dim;

#define TEST_DIMENSION(p, q, expected) \
    do \
    { \
        dim = dimension(p, q); \
        DO_TEST(dim == expected, "Expected dim(%d,%d)=%d, got %ld", \
            p, q, expected, dim); \
    } while (0)

    TEST_DIMENSION(0, 0, 1);
    TEST_DIMENSION(1, 0, 3);
    TEST_DIMENSION(0, 1, 3);
    TEST_DIMENSION(1, 1, 8);
    TEST_DIMENSION(2, 0, 6);
    TEST_DIMENSION(0, 2, 6);
    TEST_DIMENSION(3, 0, 10);
    TEST_DIMENSION(0, 3, 10);
    TEST_DIMENSION(2, 1, 15);
    TEST_DIMENSION(1, 2, 15);
    TEST_DIMENSION(2, 2, 27);

#undef TEST_DIMENSION
}

TEST(degeneracy)
{
    long d;

#define TEST_DEGENERACY(p, q, p1, q1, p2, q2, expected) \
    do \
    { \
        d = degeneracy(p, q, p1, q1, p2, q2); \
        DO_TEST(d == expected, \
            "Expected degeneracy(%d,%d,%d,%d,%d,%d)=%d, got %ld", \
            p, q, p1, q1, p2, q2, expected, d); \
    } while (0)

    /* We use the following test because it has most of the important features:
        27 x 27 -> 1 + 2x8 + 10 + 10bar + 3x27 + 28 + 28bar + 2x35 + 2x35bar
                + 2x64 + 81 + 81bar + 125 */
    TEST_DEGENERACY(0, 0, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(1, 1, 2, 2, 2, 2, 2);
    TEST_DEGENERACY(3, 0, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(0, 3, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 2, 2, 3);
    TEST_DEGENERACY(6, 0, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(0, 6, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(4, 1, 2, 2, 2, 2, 2);
    TEST_DEGENERACY(1, 4, 2, 2, 2, 2, 2);
    TEST_DEGENERACY(3, 3, 2, 2, 2, 2, 2);
    TEST_DEGENERACY(5, 2, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(2, 5, 2, 2, 2, 2, 1);
    TEST_DEGENERACY(0, 6, 2, 2, 2, 2, 1);

    /* Also test some reps which shouldn't be included */
    TEST_DEGENERACY(0, 1, 2, 2, 2, 2, 0);
    TEST_DEGENERACY(4, 0, 2, 2, 2, 2, 0);

    /* Misc. other reps */
    TEST_DEGENERACY(3, 0, 1, 1, 0, 0, 0);
    TEST_DEGENERACY(1, 1, 1, 0, 0, 1, 1);
    TEST_DEGENERACY(1, 0, 2, 0, 0, 1, 1);
    TEST_DEGENERACY(2, 1, 2, 0, 0, 1, 1);
    TEST_DEGENERACY(0, 2, 2, 0, 0, 1, 0);

#undef TEST_DEGENERACY
}

/* Check the total size of all summands after decomposing a product rep */

/* Helper function */
void test_product_rep(long p1, long q1, long p2, long q2)
{
    /* Calculate bounds on p, q */
    long upper = p1+q1+p2+q2;
    long p,q;

    long size = 0;
    long expected = dimension(p1, q1) * dimension(p2, q2);

    for (p = 0; p <= upper; ++p)
        for (q = 0; q <= upper-p; ++q)
        {
            long dim = dimension(p, q);
            long d = degeneracy(p, q, p1, q1, p2, q2);
            if (d == 0) continue;
            size += dim * d;
        }

    DO_TEST(size == expected,
        "Decomposing (%ld,%ld) x (%ld,%ld); "
        "expected total size %ld * %ld = %ld, got %ld",
        p1, q1, p2, q2,
        dimension(p1, q1), dimension(p2, q2), expected, size);
}

TEST(rep_sizes)
{
    srand(time(NULL));

    /* Do some randomised tests */
    int i;
    long p1, q1, p2, q2;
    for (i = 0; i < 10; ++i)
    {
        p1 = RANDRANGE(20);
        q1 = RANDRANGE(20);
        p2 = RANDRANGE(20);
        q2 = RANDRANGE(20);
        test_product_rep(p1, q1, p2, q2);
    }
}
