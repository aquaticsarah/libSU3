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

/* Test only that there are no crashes, but for a large number
    of representations */
TEST(isoscalars_no_crashes)
{
    isoarray* isf;

    int i;
    long p, q, p1, q1, p2, q2;
    for (i = 0; i < 10; ++i)
    {
        p = RANDRANGE(10);
        q = RANDRANGE(10);
        p1 = RANDRANGE(10);
        q1 = RANDRANGE(10);
        p2 = RANDRANGE(10);
        q2 = RANDRANGE(10);
        isf = isoscalars(p, q, p1, q1, p2, q2);
        delete isf;
    }

    DO_TEST(1, "Shouldn't happen");
}

/* Helper: Test if two isoarray objects (with the same parameters)
    contain the same values.
*/
static void check_isfs_equal(isoarray* isf1, isoarray* isf2, const char* msg)
{
    /* Check that the factor and target reps are the same */
    long p = isf1->p, q = isf1->q, p1 = isf1->p1, q1 = isf1->q1,
        p2 = isf1->p2, q2 = isf1->q2, d = isf1->d;

    bool condition = (isf2->p == p) && (isf2->q == q) &&
                     (isf2->p1 == p1) && (isf2->q1 == q1) &&
                     (isf2->p2 == p2) && (isf2->q2 == q2);

    DO_TEST(condition,
            "%s: Rep combinations differ: (%ld,%ld)x(%ld,%ld)->(%ld,%ld) "
            "vs. (%ld,%ld)x(%ld,%ld)->(%ld,%ld)\n",
            msg, p1, q1, p2, q2, p, q,
            isf2->p1, isf2->q1, isf2->p2, isf2->q2, isf2->p, isf2->q);
    if (! condition) return;

    /* Check that the values are the same */
    long n, k, l, k1, l1, k2, l2;

    for (n = 0; n < d; ++n)
        FOREACH_ISF(p, q, p1, q1, p2, q2, k, l, k1, l1, k2, l2)
        {
            if ((*isf1)(n, k, l, k1, l1, k2, l2) != (*isf2)(n, k, l, k1, l1, k2, l2))
            {

                DO_TEST(0, "%s: ISFs differ at (%ld,%ld) : (%ld,%ld) x (%ld,%ld) "
                           "in rep %ld.\n",
                        msg, k, l, k1, l1, k2, l2, d);
            }
        }

    DO_TEST(1, "\n");
}

/* Test that the various symmetry relations do the right thing */
TEST(symmetries)
{
    isoarray* isf1, * isf2, * isf3, * isf_tmp, * isf_tmp2;

    long p=2, q=1, p1=1, q1=1, p2=1, q2=0;

    isf1 = isoscalars(p, q, p1, q1, p2, q2);

    /* Test that the relations are self-inverse */
    isf_tmp = isf1->exch_12();
    isf2 = isf_tmp->exch_12();
    check_isfs_equal(isf1, isf2, "Testing 1<->2 self-inverse");
    delete isf_tmp;
    delete isf2;

    isf_tmp = isf1->exch_13bar();
    isf2 = isf_tmp->exch_13bar();
    check_isfs_equal(isf1, isf2, "Testing 1<->3bar self-inverse");
    delete isf_tmp;
    delete isf2;

    isf_tmp = isf1->exch_23bar();
    isf2 = isf_tmp->exch_23bar();
    check_isfs_equal(isf1, isf2, "Testing 2<->3bar self-inverse");
    delete isf_tmp;
    delete isf2;

    /* Test also that exch_23bar() is equivalent to the sequence
        exch_12(), exch_13bar(), exch_12()
    */
    isf_tmp = isf1->exch_12();
    isf_tmp2 = isf_tmp->exch_13bar();
    isf2 = isf_tmp2->exch_12();
    isf3 = isf1->exch_23bar();
    check_isfs_equal(isf2, isf3, "Testing 2<->3bar identity");
    delete isf_tmp;
    delete isf_tmp2;
    delete isf2;
    delete isf3;

    delete isf1;
}
