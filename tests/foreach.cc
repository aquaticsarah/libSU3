/* libSU3: Test that the FOREACH_* macros do not miss anything.
    We do this by counting the number of nonzero entries, once using the
    FOREACH macros and once by explicitly iterating over every possible
    coupling, then checking that the counts are equal.
*/

#include "SU3.h"
#include "test.h"

TEST(foreach_isf)
{
    long d;
    long n, k, l, k1, l1, k2, l2;
    long p = 2, q = 2, p1 = 2, q1 = 2, p2 = 2, q2 = 2;

    /* Counters for how many nonzero values we have seen. */
    long count_foreach = 0, count_all = 0;

    isoarray* isf = isoscalars(p, q, p1, q1, p2, q2);
    d = degeneracy(p, q, p1, q1, p2, q2);

    /* Count nonzero values using the FOREACH macro */
    for (n = 0; n < d; ++n)
        FOREACH_ISF(p, q, p1, q1, p2, q2, k, l, k1, l1, k2, l2)
            if ((*isf)(n, k, l, k1, l1, k2, l2) != 0)
                ++count_foreach;


    /* Count nonzero values explicitly */
    for (n = 0; n < d; ++n)
    for (k = q; k <= p+q; ++k)
        for (l = 0; l <= q; ++l)
            for (k1 = q1; k1 <= p1+q1; ++k1)
                for (l1 = 0; l1 <= q1; ++l1)
                    for (k2 = q2; k2 <= p2+q2; ++k2)
                        for (l2 = 0; l2 <= q2; ++l2)
                            if ((*isf)(n, k, l, k1, l1, k2, l2) != 0)
                                ++count_all;

    DO_TEST(count_foreach == count_all,
            "FOREACH_ISF does not cover all nonzero couplings: Covers only %ld out of %ld.",
            count_foreach, count_all);

    delete isf;
}

TEST(foreach_cgc)
{
    long d;
    long n, k, l, m, k1, l1, m1, k2, l2, m2;
    long p = 2, q = 2, p1 = 2, q1 = 2, p2 = 2, q2 = 2;

    /* Counters for how many nonzero values we have seen. */
    long count_foreach = 0, count_all = 0;

    cgarray* cg = clebsch_gordans(p, q, p1, q1, p2, q2);
    d = degeneracy(p, q, p1, q1, p2, q2);

    /* Count nonzero values using the FOREACH macro */
    for (n = 0; n < d; ++n)
        FOREACH_CGC(p, q, p1, q1, p2, q2, k, l, m, k1, l1, m1, k2, l2, m2)
            if ((*cg)(n, k, l, m, k1, l1, m1, k2, l2, m2) != 0)
                ++count_foreach;


    /* Count nonzero values explicitly */
    for (n = 0; n < d; ++n)
    for (k = q; k <= p+q; ++k)
        for (l = 0; l <= q; ++l)
            for (m = l; m <= k; ++m)
                for (k1 = q1; k1 <= p1+q1; ++k1)
                    for (l1 = 0; l1 <= q1; ++l1)
                        for (m1 = l1; m1 <= k1; ++m1)
                            for (k2 = q2; k2 <= p2+q2; ++k2)
                                for (l2 = 0; l2 <= q2; ++l2)
                                    for (m2 = l2; m2 <= k2; ++m2)
                                        if ((*cg)(n, k, l, m, k1, l1, m1, k2, l2, m2) != 0)
                                            ++count_all;

    DO_TEST(count_foreach == count_all,
            "FOREACH_CGC does not cover all nonzero couplings: Covers only %ld out of %ld.",
            count_foreach, count_all);

    delete cg;
}
