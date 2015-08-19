/* libSU3: State-of-higest-weight (SHW) determination.
   Note: This also does degeneracy resolution */

#include <stdio.h>
#include <assert.h>
#include <stdexcept>

#include "SU3_internal.h"

isoscalar_context::isoscalar_context(isoarray* isf, long p, long q, long p1,
            long q1, long p2, long q2) : p(p), q(q), p1(p1),
            q1(q1), p2(p2), q2(q2), isf(isf)
{
    d = degeneracy(p, q, p1, q1, p2, q2);
    A = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;
}

/* Use the C and D recursion relations to step along the
   k and l axes within a multiplet.
   We fix l=0 in step_k_down, since this can *only* step down to such states.

   The arguments identify the state to be calclated, *not* the values
   of k1,l1,k2,l2 used in the recursion relation itself. */
void isoscalar_context::step_k_down(long n, long k, long k1, long l1,
                                    long k2, long l2)
{
    sqrat beta, d1, d2, d3;
    d_coefficients(k, 0L, k1, l1, k2, l2, beta, d1, d2, d3);

    /* Calculate the value at (k,0,k1,l1,k2,l2) using surrounding values */
    (*isf)(n, k, 0L, k1, l1, k2, l2) = (d1 * (*isf)(n, k+1, 0L, k1+1, l1, k2, l2)
                                       +d2 * (*isf)(n, k+1, 0L, k1, l1, k2+1, l2)
                                       +d3 * (*isf)(n, k+1, 0L, k1, l1, k2, l2+1)
                                       ) * beta;
}

void isoscalar_context::step_l_up(long n, long k, long l, long k1,
                                    long l1, long k2, long l2)
{
    sqrat alpha, c1, c2, c3, c4;

    c_coefficients(k, l, k1, l1, k2, l2, alpha, c1, c2, c3, c4);

    /* Calculate the value at (k,l,k1,l1,k2,l2) using surrounding values */
    (*isf)(n, k, l, k1, l1, k2, l2) = (c1 * (*isf)(n, k+1, l-1, k1, l1, k2, l2)
                                      +c2 * (*isf)(n, k, l-1, k1, l1-1, k2, l2)
                                      +c3 * (*isf)(n, k, l-1, k1, l1, k2-1, l2)
                                      +c4 * (*isf)(n, k, l-1, k1, l1, k2, l2-1)
                                      ) * alpha;
}

/* Internal function: Calculate the isoscalar factors for a particular
    combination of reps. */
int isoscalar_context::calc_isoscalars()
{
    /* Calculate couplings to the state of highest weight (k=p+q, l=0).
        This may fail, in which case we return to the top-level function
        which will try using symmetries to convert into a solvable problem.
    */
    if (! this->calc_shw())
        return 0;

    /* Then fill in the rest of the couplings */
    long n, k, l, k1, l1, k2, l2;
    for (n = 0; n < d; ++n)
    {
        /* Fill in the rest of the k=p+q row */
        k = p+q;
        for (l = 1; l <= q; ++l)
            for (k1 = q1; k1 <= p1+q1; ++k1)
                for (l1 = 0; l1 <= q1; ++l1)
                    for (k2 = q2; k2 <= p2+q2; ++k2)
                    {
                        l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3 - (k1 + l1 + k2 - k - l);
                        if ((l2 < 0) || (l2 > q2)) continue;

                        step_l_up(n, k, l, k1, l1, k2, l2);
                    }

        /* Fill in couplings to one state on this row */
        for (k = p+q-1; k >= q; --k)
        {
            for (k1 = q1; k1 <= p1+q1; ++k1)
                for (l1 = 0; l1 <= q1; ++l1)
                    for (k2 = q2; k2 <= p2+q2; ++k2)
                    {
                        l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3 - (k1 + l1 + k2 - k);
                        if ((l2 < 0) || (l2 > q2)) continue;

                        step_k_down(n, k, k1, l1, k2, l2);
                    }

            /* Fill in the rest of the row */
            for (l = 1; l <= q; ++l)
                for (k1 = q1; k1 <= p1+q1; ++k1)
                    for (l1 = 0; l1 <= q1; ++l1)
                        for (k2 = q2; k2 <= p2+q2; ++k2)
                        {
                            l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3 - (k1 + l1 + k2 - k - l);
                            if ((l2 < 0) || (l2 > q2)) continue;

                            step_l_up(n, k, l, k1, l1, k2, l2);
                        }
        }
    }

    return 1;
}

/* Main calculation function */
isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (! d) return NULL; /* Ignore reps of zero degeneracy */

    isoarray* isf = new isoarray(p,q,p1,q1,p2,q2);
    isoscalar_context* ctx = new isoscalar_context(isf,p,q,p1,q1,p2,q2);

    if (ctx->calc_isoscalars())
    {
        delete ctx;
        return isf;
    }

    /* If the direct algorithm fails, we try using the 1 <-> 3bar symmetry.
        This relates the isoscalar factors for r1 x r2 -> R to those for
        Rbar x r2 -> r1bar.
    */
    delete ctx;

#define SIGN(v) ((((v) % 2) == 0) ? 1 : -1)

    /* TODO: Check that degeneracy is the same under symmetry relations */
    assert(d == degeneracy(q1, p1, q, p, p2, q2));
    isoarray* new_isf = new isoarray(q1, p1, q, p, p2, q2);
    ctx = new isoscalar_context(new_isf, q1, p1, q, p, p2, q2);

    if (ctx->calc_isoscalars())
    {
        /* Use the symmetry relations to fill out the isoscalar factors for
           the reps we wanted originally */
        long n, k, l, k1, l1, k2, l2;
        for (n = 0; n < d; ++n)
            for (k = q; k <= p+q; ++k)
                for (l = 0; l <= q; ++l)
                    for (k1 = q1; k1 <= p1+q1; ++k1)
                        for (l1 = 0; l1 <= q1; ++l1)
                            for (k2 = q2; k2 <= p2+q2; ++k2)
                            {
                                l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3 - (k1 + l1 + k2 - k - l);
                                if ((l2 < 0) || (l2 > q2)) continue;

                                (*isf)(n, k, l, k1, l1, k2, l2) =
                                    SIGN(l2)
                                  * sqrat((p+1)*(q+1)*(p+q+2)*(k1-l1+1), (p1+1)*(q1+1)*(p1+q1+2)*(k-l+1))
                                  * (*new_isf)(n, p1+q1-l1, p1+q1-k1, p+q-l, p+q-k, k2, l2);
                            }

        delete new_isf;
        return isf;
    }

    /* If that fails, combine with the r1 <-> r2 exchange symmetry.
        This results in a relation between the ISFs for r1 x r2 -> R
        and those for Rbar x r1 -> r2bar */
    delete ctx;

    assert(d == degeneracy(q2, p2, q, p, p1, q1));
    new_isf = new isoarray(q2, p2, q, p, p1, q1);
    ctx = new isoscalar_context(new_isf, q2, p2, q, p, p1, q1);

    if (ctx->calc_isoscalars())
    {
        /* For this case, we need to calculate an overall phase factor */
        long x = p1 + p2 - p;
        long y = q1 + q2 - q;
        long xi_1 = SIGN(x + y + max(x, y));

        /* Use the symmetry relations to fill out the isoscalar factors for
           the reps we wanted originally */
        long n, k, l, k1, l1, k2, l2;
        for (n = 0; n < d; ++n)
            for (k = q; k <= p+q; ++k)
                for (l = 0; l <= q; ++l)
                    for (k1 = q1; k1 <= p1+q1; ++k1)
                        for (l1 = 0; l1 <= q1; ++l1)
                            for (k2 = q2; k2 <= p2+q2; ++k2)
                            {
                                l2 = (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3 - (k1 + l1 + k2 - k - l);
                                if ((l2 < 0) || (l2 > q2)) continue;

                                (*isf)(n, k, l, k1, l1, k2, l2) =
                                    xi_1
                                  * SIGN((k1 - l1 + k2 - l2 - k + l)/2)
                                  * SIGN(l1)
                                  * sqrat((p+1)*(q+1)*(p+q+2)*(k2-l2+1), (p2+1)*(q2+1)*(p2+q2+2)*(k-l+1))
                                  * (*new_isf)(n, p2+q2-l2, p2+q2-k2, p+q-l, p+q-k, k1, l1);
                            }

        delete new_isf;
        return isf;
    }

    /* If we get here, nothing has worked, so throw an error */
    delete isf;
    throw std::logic_error("Calculation of ISFs failed");
}
