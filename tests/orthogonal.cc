/* libSU3: Test that the Clebsch-Gordan coefficients form an orthogonal matrix.
    Note that the CGCs are effectively a (real) change of basis, so they should
    form an orthogonal matrix when arranged appropriately
*/

#include <assert.h>

#include "SU3.h"
#include "test.h"

/* Helper: Fill out some number of rows of the SU(3) Clebsch-Gordan matrix */
static void fill_rows(long p, long q, long p1, long q1,
                long p2, long q2, sqrat* cgc_matrix, long& index)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (! d) return;

    cgarray* cgcs = clebsch_gordans(p, q, p1, q1, p2, q2);

    /* Copy all the values across to the main matrix.
        Note: Rows are (effectively) labelled by n,k,l,m
        and columns by k1,l1,m1,k2,l2,m2
    */
    long n, k, l, m, k1, l1, m1, k2, l2, m2;
    for (n = 0; n < d; ++n)
        /* We don't use FOREACH_CGC here, as we explicitly want to iterate over
            *every* CGC, even the ones which definitely can't couple, so that we
            fill out entire rows.
        */
        for (k = q; k <= p+q; ++k)
            for (l = 0; l <= q; ++l)
                for (m = l; m <= k; ++m)
                    for (k1 = q1; k1 <= p1+q1; ++k1)
                        for (l1 = 0; l1 <= q1; ++l1)
                            for (m1 = l1; m1 <= k1; ++m1)
                                for (k2 = q2; k2 <= p2+q2; ++k2)
                                    for (l2 = 0; l2 <= q2; ++l2)
                                        for (m2 = l2; m2 <= k2; ++m2)
                                        {
                                            cgc_matrix[index++] =
                                                (*cgcs)(n, k, l, m, k1, l1, m1, k2, l2, m2);
                                        }

    delete cgcs;
}

static int check_orthogonal(long p1, long q1, long p2, long q2)
{
    /* Allocate matrix */
    long N = dimension(p1, q1) * dimension(p2, q2);
    sqrat* cgc_matrix = new sqrat[N*N];

    /* Fill out matrix by calculating CGCs for each summand */
    long index = 0;

    long upper = p1+q1+p2+q2; // Loose bound on which reps can be summands
    long p, q;

    for (p = 0; p <= upper; ++p)
        for (q = 0; q <= upper-p; ++q)
            fill_rows(p, q, p1, q1, p2, q2, cgc_matrix, index);

    /* Make sure the table is fully filled out */
    assert(index == N*N);

    /* Now check that the matrix is orthogonal. We only need to check
        the inner product of each row with previous rows (and itself)
    */
    long m, n;
    for (m = 0; m < N; ++m)
        for (n = 0; n <= m; ++n)
        {
            long expected = (m == n) ? 1 : 0;

            /* Calculate the inner product of rows m and n */
            sqrat res = 0;
            long idx;

            for (idx = 0; idx < N; ++idx)
                res += cgc_matrix[m*N + idx] * cgc_matrix[n*N + idx];

            /* If the rows are not orthonormal, fail the test */
            if (res != expected)
            {
                delete[] cgc_matrix;
                return 0;
            }
        }

    /* If we get here, pass the test */
    delete[] cgc_matrix;
    return 1;
}

TEST(orthogonal)
{
    DO_TEST(check_orthogonal(0, 0, 0, 1),
            "CGCs for (0,0)x(0,1) not orthonormal");
    DO_TEST(check_orthogonal(1, 0, 0, 1),
            "CGCs for (1,0)x(0,1) not orthonormal");
    DO_TEST(check_orthogonal(1, 0, 2, 0),
            "CGCs for (1,0)x(2,0) not orthonormal");
    DO_TEST(check_orthogonal(1, 1, 1, 1),
            "CGCs for (1,1)x(1,1) not orthonormal");
}
