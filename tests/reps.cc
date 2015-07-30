/* libSU3: Tests for irrep information */

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

#define TEST_DEGENERACY(p1, q1, p2, q2, p, q, expected) \
    do \
    { \
        d = degeneracy(p1, q1, p2, q2, p, q); \
        DO_TEST(d == expected, \
            "Expected degeneracy(%d,%d,%d,%d,%d,%d)=%d, got %ld", \
            p1, q1, p2, q2, p, q, expected, d); \
    } while (0)

    /* We use the following test because it has most of the important features:
       27 x 27 -> 1 + 2x8 + 10 + 10bar + 3x27 + 28 + 28bar + 2x35 + 2x35bar
                + 2x64 + 81 + 81bar + 125 */
    TEST_DEGENERACY(2, 2, 2, 2, 0, 0, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 1, 1, 2);
    TEST_DEGENERACY(2, 2, 2, 2, 3, 0, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 0, 3, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 2, 2, 3);
    TEST_DEGENERACY(2, 2, 2, 2, 6, 0, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 0, 6, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 4, 1, 2);
    TEST_DEGENERACY(2, 2, 2, 2, 1, 4, 2);
    TEST_DEGENERACY(2, 2, 2, 2, 3, 3, 2);
    TEST_DEGENERACY(2, 2, 2, 2, 5, 2, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 2, 5, 1);
    TEST_DEGENERACY(2, 2, 2, 2, 0, 6, 1);

    /* Also test some reps which shouldn't be included */
    TEST_DEGENERACY(2, 2, 2, 2, 0, 1, 0);
    TEST_DEGENERACY(2, 2, 2, 2, 4, 0, 0);

    /* Misc. other reps */
    TEST_DEGENERACY(1, 1, 0, 0, 3, 0, 0);
    TEST_DEGENERACY(1, 0, 0, 1, 1, 1, 1);
    TEST_DEGENERACY(2, 0, 0, 1, 1, 0, 1);
    TEST_DEGENERACY(2, 0, 0, 1, 2, 1, 1);
    TEST_DEGENERACY(2, 0, 0, 1, 0, 2, 0);

#undef TEST_DEGENERACY
}
