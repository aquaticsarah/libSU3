/* libSU3: Calculation of isoscalar factors.
    This file does everything except couplings to the state of highest weight.
    For that, see shw.cc */

#include <stdio.h>
#include <assert.h>
#include <stdexcept>

#include "SU3_internal.h"

isoscalar_context::isoscalar_context(long p, long q, long p1,
            long q1, long p2, long q2, long d, sqrat* coefficients)
            : p(p), q(q), p1(p1), q1(q1), p2(p2), q2(q2), d(d),
            coefficients(coefficients)
{
    A = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;
}

/* Functions to get/set particular isoscalar factors.

   Notes on indexing:
   - For each of the three reps involved, we have the following ranges:
      q <= k <= p+q (for a total of p+1 possible values of k)
      0 <= l <= q   (for a total of q+1 possible values of l)
     The value of 'l' can be used as an index directly, but that for
     'k' needs to be shifted down by q.

   - In order to simplify the usage of recursion relations, we allow code to
     request values one space "off the edge" of the valid range (eg, with
     l1=-1). We only allow this when getting values, not when setting them.

   - All reps in a degenerate set need to be processed at once, so we
     allocate storage all at once.

   - We can infer the value of l2 from the values of k,l,k1,l1,k1 using
     hypercharge conservation. As such, we don't need an axis for l2.
     We do check that the value provided is correct; TODO: make this optional.
*/
sqrat isoscalar_context::isf(long n, long k, long l, long k1, long l1,
                            long k2, long l2)
{
    /* Bounds checks; here we allow one space extra around the valid range */
    assert((n >= 0) && (n < d));
    assert((k >= q-1) && (k <= p+q+1));
    assert((l >= -1) && (l <= q+1));
    assert((k1 >= q1-1) && (k1 <= p1+q1+1));
    assert((l1 >= -1) && (l1 <= q1+1));
    assert((k2 >= q2-1) && (k2 <= p2+q2+1));
    assert((l2 >= -1) && (l2 <= q2+1));

    /* Check hypercharge conservation */
    assert(k1+l1+k2+l2-k-l == (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3);

    /* Return zero for values outside of range */
    if (   (k  < q ) || (k  > p +q ) || (l  < 0) || (l  > q )
        || (k1 < q1) || (k1 > p1+q1) || (l1 < 0) || (l1 > q1)
        || (k2 < q2) || (k2 > p2+q2) || (l2 < 0) || (l2 > q2))
        return 0;

    /* Otherwise, get the value from our coefficient array */
    size_t index = ((((n * (p+1) + k-q) * (q+1) + l) * (p1+1) + k1-q1)
                    * (q1+1) + l1) * (p2+1) + k2-q2;
    return coefficients[index];
}

void isoscalar_context::set_isf(long n, long k, long l, long k1, long l1,
                            long k2, long l2, sqrat value)
{
    /* Bounds checks */
    assert((n >= 0) && (n < d));
    assert((k >= q) && (k <= p+q));
    assert((l >= 0) && (l <= q));
    assert((k1 >= q1) && (k1 <= p1+q1));
    assert((l1 >= 0) && (l1 <= q1));
    assert((k2 >= q2) && (k2 <= p2+q2));
    assert((l2 >= 0) && (l2 <= q2));

    /* Check hypercharge conservation */
    assert(k1+l1+k2+l2-k-l == (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3);

    size_t index = ((((n * (p+1) + k-q) * (q+1) + l) * (p1+1) + k1-q1)
                    * (q1+1) + l1) * (p2+1) + k2-q2;
    coefficients[index] = value;
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
    d_coefficients(k, k1, l1, k2, l2, beta, d1, d2, d3);

    /* Calculate the value at (k,0,k1,l1,k2,l2) using surrounding values */
    sqrat res = (d1 * isf(n, k+1, 0L, k1+1, l1, k2, l2)
                +d2 * isf(n, k+1, 0L, k1, l1, k2+1, l2)
                +d3 * isf(n, k+1, 0L, k1, l1, k2, l2+1)) * beta;
    set_isf(n, k, 0L, k1, l1, k2, l2, res);
}

void isoscalar_context::step_l_up(long n, long k, long l, long k1,
                                    long l1, long k2, long l2)
{
    sqrat alpha, c1, c2, c3, c4;

    c_coefficients(k, l, k1, l1, k2, l2, alpha, c1, c2, c3, c4);

    /* Calculate the value at (k,l,k1,l1,k2,l2) using surrounding values */
    sqrat res = (c1 * isf(n, k+1, l-1, k1, l1, k2, l2)
                +c2 * isf(n, k, l-1, k1, l1-1, k2, l2)
                +c3 * isf(n, k, l-1, k1, l1, k2-1, l2)
                +c4 * isf(n, k, l-1, k1, l1, k2, l2-1)) * alpha;
    set_isf(n, k, l, k1, l1, k2, l2, res);
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

    size_t size = d * (p+1) * (q+1) * (p1+1) * (q1+1) * (p2+1);
    sqrat* coefficients = new sqrat[size];
    isoscalar_context* ctx = new isoscalar_context(p, q, p1, q1, p2, q2,
                                                    d, coefficients);

    /* This is what we will return */
    isoarray* res = new isoarray(p, q, p1, q1, p2, q2, d, coefficients);

    if (ctx->calc_isoscalars())
    {
        delete ctx;
        return res;
    }

    /* If the direct algorithm fails, we try using the 1 <-> 3bar symmetry.
        This relates the isoscalar factors for r1 x r2 -> R to those for
        Rbar x r2 -> r1bar.
    */

    /* TODO: Check that degeneracy is the same under symmetry relations */
    assert(d == degeneracy(q1, p1, q, p, p2, q2));
    size_t alt_size = d * (q1+1) * (p1+1) * (q+1) * (p+1) * (p2+1);
    sqrat* alt_coefficients = new sqrat[alt_size];
    isoscalar_context* alt_ctx = new isoscalar_context(q1, p1, q, p, p2, q2,
                                                        d, alt_coefficients);

    if (alt_ctx->calc_isoscalars())
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

                                ctx->set_isf(n, k, l, k1, l1, k2, l2,
                                    SIGN(l2)
                                  * sqrat((p+1)*(q+1)*(p+q+2)*(k1-l1+1), (p1+1)*(q1+1)*(p1+q1+2)*(k-l+1))
                                  * alt_ctx->isf(n, p1+q1-l1, p1+q1-k1, p+q-l, p+q-k, k2, l2));
                            }

        delete alt_ctx;
        delete[] alt_coefficients;
        delete ctx;
        return res;
    }

    /* If that fails, combine with the r1 <-> r2 exchange symmetry.
        This results in a relation between the isfs for r1 x r2 -> R
        and those for Rbar x r1 -> r2bar */
    delete alt_ctx;
    delete[] alt_coefficients;

    assert(d == degeneracy(q2, p2, q, p, p1, q1));
    alt_size = d * (q2+1) * (p2+1) * (q+1) * (p+1) * (p1+1);
    alt_coefficients = new sqrat[alt_size];
    alt_ctx = new isoscalar_context(q2, p2, q, p, p1, q1, d, alt_coefficients);

    if (alt_ctx->calc_isoscalars())
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

                                ctx->set_isf(n, k, l, k1, l1, k2, l2,
                                    xi_1
                                  * SIGN((k1 - l1 + k2 - l2 - k + l)/2)
                                  * SIGN(l1)
                                  * sqrat((p+1)*(q+1)*(p+q+2)*(k2-l2+1), (p2+1)*(q2+1)*(p2+q2+2)*(k-l+1))
                                  * alt_ctx->isf(n, p2+q2-l2, p2+q2-k2, p+q-l, p+q-k, k1, l1));
                            }

        delete alt_ctx;
        delete[] alt_coefficients;
        delete ctx;
        return res;
    }

    /* If we get here, nothing has worked, so throw an error */
    delete alt_ctx;
    delete[] alt_coefficients;
    delete ctx;
    delete res; // Also deletes 'coefficients'
    throw std::logic_error("Calculation of ISFs failed");
}
