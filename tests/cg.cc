/* libSU3: Tests for SU(3) Clebsch-Gordans */

#include "SU3.h"
#include "test.h"

#define TEST_CG(k, l, m, k1, l1, m1, k2, l2, m2, p, q) \
    TEST_EQ_SQRAT((*cg)(0, k, l, m, k1, l1, m1, k2, l2, m2), p, q, \
        "<%d,%d,%d| (|%d,%d,%d> x |%d,%d,%d>)", \
        k, l, m, k1, l1, m1, k2, l2, m2)

TEST(cg)
{
    /* Test 3 x 3bar -> 1 + 8.
        Note that the ISFs for this are tested in tests/isoscalars.cc, so
        if there is a bug here we know it is in the conversion ISFs -> CGs
    */
    cgarray* cg = clebsch_gordans(0, 0, 1, 0, 0, 1);

    //      k  l  m k1 l1 m1 k2 l2 m2     p  q
    TEST_CG(0, 0, 0, 0, 0, 0, 1, 1, 1,   -1, 3);
    TEST_CG(0, 0, 0, 1, 0, 0, 1, 0, 1,   -1, 3);
    TEST_CG(0, 0, 0, 1, 0, 1, 1, 0, 0,    1, 3);

    delete cg;
    cg = clebsch_gordans(1, 1, 1, 0, 0, 1);

    //      k  l  m k1 l1 m1 k2 l2 m2     p  q
    TEST_CG(1, 1, 1, 0, 0, 0, 1, 1, 1,    2, 3);
    TEST_CG(1, 1, 1, 1, 0, 0, 1, 0, 1,   -1, 6);
    TEST_CG(1, 1, 1, 1, 0, 1, 1, 0, 0,    1, 6);

    TEST_CG(2, 0, 1, 0, 0, 0, 1, 1, 1,    0, 1);
    TEST_CG(2, 0, 1, 1, 0, 0, 1, 0, 1,    1, 2);
    TEST_CG(2, 0, 1, 1, 0, 1, 1, 0, 0,    1, 2);

    TEST_CG(2, 1, 1, 1, 0, 0, 1, 1, 1,    1, 1);
    TEST_CG(2, 1, 2, 1, 0, 1, 1, 1, 1,    1, 1);
    TEST_CG(2, 0, 0, 1, 0, 0, 1, 0, 0,    1, 1);
    TEST_CG(2, 0, 2, 1, 0, 1, 1, 0, 1,    1, 1);
    TEST_CG(1, 0, 0, 0, 0, 0, 1, 0, 0,    1, 1);
    TEST_CG(1, 0, 1, 0, 0, 0, 1, 0, 1,    1, 1);

    delete cg;
}
