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
    Indexing is done just like in src/isoarray.cc, with the exception that we
    allow values which are one space "off the edge" (eg, with l=-1), returning
    0 for those couplings. This is because doing so greatly simplifies the
    main calculation code.
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

    /* Check hypercharge conservation.
        Note that, unlike isoarray/cgarray, we error out if this isn't
        satisfied. This is because library-internal code should *always*
        obey hypercharge conservation, but external code need not do so.
    */
    assert(k1+l1+k2+l2-k-l == (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3);

    /* If compiled with -DNDEBUG, this line is to remove a compiler warning
        about l2 being unused */
    (void)l2;

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

    /* If compiled with -DNDEBUG, this line is to remove a compiler warning
        about l2 being unused */
    (void)l2;

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
void isoscalar_context::calc_isoscalars()
{
    /* Calculate couplings to the state of highest weight (k=p+q, l=0). */
    this->calc_shw();

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
}

/* Internal: Can we calculate the ISFs for a given set of reps directly? */
static int can_calculate(long p, long q, long p1, long q1, long p2, long q2,
                         long d)
{
    long A  = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;
    long Ar = (2*q1 + 2*q2 + 4*p1 + 4*p2 + q - p)/3;
    long smax = min(A, Ar);
    long smin = max(p + q, abs(2*q1 + 2*q2 - A));

    /* The only way a calculation can fail is if the step down to s = smax-2*d fails.
        This happens iff all of the following conditions fail. Each corresponds to
        some reason the calculation at that level succeeds.
    */
    if (    (smax - 2*d < smin)                     // Nothing to calculate
        || ((A + (smax - 2*d))/2 - (p2+q2) < q1)    // step_k1_up succeeds
        || ((A - (smax - 2*d))/2 > q1)              // step_l1_down succeeds
        || ((A + (smax - 2*d))/2 - q2 < p1+q1)      // step_k1_down succeeds
        || ((A - (smax - 2*d))/2 - q2 > 0))         // step_l1_up succeeds
            return 1;
    else
        return 0;
}

/* Internal: Calculate values for one irrep combination, without trying
    the symmetry relations. Returns NULL on failure.
*/
isoarray* isoscalars_single(long p, long q, long p1, long q1,
                                    long p2, long q2, long d)
{
    /* Check that a direct calculation will succeed, before we do it */
    if (! can_calculate(p, q, p1, q1, p2, q2, d))
        return NULL;

    size_t size = d * (p+1) * (q+1) * (p1+1) * (q1+1) * (p2+1);
    sqrat* coefficients = new sqrat[size];
    isoscalar_context* ctx = new isoscalar_context(p, q, p1, q1, p2, q2,
                                                    d, coefficients);
    isoarray* isf = new isoarray(p, q, p1, q1, p2, q2, d, coefficients);

    ctx->calc_isoscalars();
    delete ctx;
    return isf;
}

/* Main calculation function */
isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (! d) return NULL; /* Ignore reps of zero degeneracy */

    /* Try to calculate directly */
    isoarray* isf = isoscalars_single(p, q, p1, q1, p2, q2, d);
    if (isf)
        return isf;

    /* If the direct algorithm fails, we try using the 1 <-> 3bar symmetry.
        This relates the isoscalar factors for r1 x r2 -> R to those for
        Rbar x r2 -> r1bar.
    */
    isf = isoscalars_single(q1, p1, q, p, p2, q2, d);
    if (isf)
    {
        isoarray* new_isf = isf->exch_13bar();
        delete isf;
        return new_isf;
    }

    /* If that fails, combine with the r1 <-> r2 exchange symmetry.
        This results in a relation between the isfs for r1 x r2 -> R
        and those for Rbar x r1 -> r2bar */
    isf = isoscalars_single(q2, p2, p1, q1, q, p, d);
    if (isf)
    {
        isoarray* new_isf = isf->exch_23bar();
        delete isf;
        return new_isf;
    }

    /* If we get here, nothing has worked, so throw an error */
    throw std::logic_error("Calculation of ISFs failed. "
                            "please report this as a bug in libSU3.");
}

/* Wrapper around the above to provide an array of Clebsch-Gordans instead */
cgarray* clebsch_gordans(long p, long q, long p1, long q1, long p2, long q2)
{
    isoarray* isf = isoscalars(p, q, p1, q1, p2, q2);
    cgarray* cg = isf->to_cgarray();

    delete isf;
    return cg;
}
