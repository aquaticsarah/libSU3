/* libSU3: Couplings to the state of highest weight
    (that is, the state of highest I) */

#include <stdio.h>
#include <stdexcept>

#include "SU3_internal.h"

/* Use the A and B recursion relations to step along the
    k1 and l1 axes within a plane.

    The arguments identify the state to be calclated, *not* the values
    of k1,l1,k2,l2 used in the recursion relation itself.
*/
void isoscalar_context::step_k1_up(long n, long s, long k1, long l1)
{
    long k2 = (A+s)/2 - k1, l2 = (A-s)/2 - l1;

    /* Calculate coefficients for the A recursion relation */
    sqrat a1, a2, a3, a4;
    a_coefficients(k1, l1, k2+1, l2, a1, a2, a3, a4);

    /* Calculate the value at (k1,l1,k2,l2) using surrounding values */
    sqrat res = (-a1 * isf(n, p+q, 0, k1-1, l1, k2+1, l2)
                 -a3 * isf(n, p+q, 0, k1, l1-1, k2+1, l2)
                 -a4 * isf(n, p+q, 0, k1, l1, k2+1, l2-1)) / a2;
    set_isf(n, p+q, 0, k1, l1, k2, l2, res);
}

void isoscalar_context::step_k1_down(long n, long s, long k1, long l1)
{
    long k2 = (A+s)/2 - k1, l2 = (A-s)/2 - l1;

    /* Calculate coefficients for the A recursion relation */
    sqrat a1, a2, a3, a4;
    a_coefficients(k1+1, l1, k2, l2, a1, a2, a3, a4);

    /* Calculate the value at (k1,l1,k2,l2) using surrounding values */
    sqrat res = (-a2 * isf(n, p+q, 0, k1+1, l1, k2-1, l2)
                 -a3 * isf(n, p+q, 0, k1+1, l1-1, k2, l2)
                 -a4 * isf(n, p+q, 0, k1+1, l1, k2, l2-1)) / a1;
    set_isf(n, p+q, 0, k1, l1, k2, l2, res);
}

void isoscalar_context::step_l1_up(long n, long s, long k1, long l1)
{
    long k2 = (A+s)/2 - k1, l2 = (A-s)/2 - l1;

    /* Calculate coefficients for the B recursion relation */
    sqrat b1, b2, b3, b4;
    b_coefficients(k1, l1-1, k2, l2, b1, b2, b3, b4);

    /* Calculate the value at (k1,l1,k2,l2) using surrounding values */
    sqrat res = (-b1 * isf(n, p+q, 0, k1+1, l1-1, k2, l2)
                 -b2 * isf(n, p+q, 0, k1, l1-1, k2+1, l2)
                 -b4 * isf(n, p+q, 0, k1, l1-1, k2, l2+1)) / b3;
    set_isf(n, p+q, 0, k1, l1, k2, l2, res);
}

void isoscalar_context::step_l1_down(long n, long s, long k1, long l1)
{
    long k2 = (A+s)/2 - k1, l2 = (A-s)/2 - l1;

    /* Calculate coefficients for the B recursion relation */
    sqrat b1, b2, b3, b4;
    b_coefficients(k1, l1, k2, l2-1, b1, b2, b3, b4);

    /* Calculate the value at (k1,l1,k2,l2) using surrounding values */
    sqrat res = (-b1 * isf(n, p+q, 0, k1+1, l1, k2, l2-1)
                 -b2 * isf(n, p+q, 0, k1, l1, k2+1, l2-1)
                 -b3 * isf(n, p+q, 0, k1, l1+1, k2, l2-1)) / b4;
    set_isf(n, p+q, 0, k1, l1, k2, l2, res);
}

/* Step down from one plane (at s+2) to the next plane (at s).
    In principle this can fail, in which case you need to use the exchange
    symmetries in order to calculate the values. This should never happen with
    this function, however, as isoscalars() has logic for deciding whether this
    can happen *before* it picks which ISFs to calculate.

    Note: When doing this, there are two independent choices we can make:

    * We can try to calculate the value at (k1max,l1min) or at (k1min,l1max)
        (these require different steps afterwards)

    * We can try recursion relation A or recursion relation B
        (these require the same steps afterwards)

    Each of the four combinations works in different cases, so we intelligently
    choose between them.

    Internally, we use the step_k1_up/down and step_l1_up/down functions,
    but we step "from" a non-existent state (with k1,l1 not in the valid range)
    to the state we want. This works because 'isoscalar_array' allows us to
    request (certain) non-existent states and just returns 0 for the coupling
    coefficient. This is exactly what we need for the stepping to work properly.
*/
void isoscalar_context::step_s_down(long n, long s)
{
    /* Calculate unconstrained versions of the min/max values */
    long k1min_u = (A + s)/2 - (p2+q2);
    long k1max_u = (A + s)/2 - q2;
    long l1min_u = (A - s)/2 - q2;
    long l1max_u = (A - s)/2;

    long k1min = max(q1, k1min_u);
    long k1max = min(p1+q1, k1max_u);
    long l1min = max(0, l1min_u);
    long l1max = min(q1, l1max_u);

    /* Each recursion relation has conditions for being valid. We test these here. */
    if (k1min_u < q1)
        step_k1_up(n, s, k1min, l1max);
    else if (l1max_u > q1)
        step_l1_down(n, s, k1min, l1max);
    else
    {
        /* Stepping to (k1min, l1max) failed, so try (k1max, l1min) */
        if (k1max_u < p1+q1)
            step_k1_down(n, s, k1max, l1min);
        else if (l1min_u > 0)
            step_l1_up(n, s, k1max, l1min);
        else
            throw std::logic_error("Couldn't use any recursion relations. "
                                    "Please report this as a bug in libSU3.\n");

        /* If we get here, we succeded at (k1max, l1min).
            So fill out the rest of this plane from this point. */
        long k1, l1;
        for (l1 = l1min+1; l1 <= l1max; ++l1)
            step_l1_up(n, s, k1max, l1);

        for (k1 = k1max-1; k1 >= k1min; --k1)
        {
            step_k1_down(n, s, k1, l1min);
            for (l1 = l1min+1; l1 <= l1max; ++l1)
                step_l1_up(n, s, k1, l1);
        }

        /* Avoid falling through to the code below */
        return;
    }

    /* If we get here, we succeded at (k1min, l1max).
        So fill out the rest of this plane from this point. */
    long k1, l1;
    for (l1 = l1max-1; l1 >= l1min; --l1)
        step_l1_down(n, s, k1min, l1);

    for (k1 = k1min+1; k1 <= k1max; ++k1)
    {
        step_k1_up(n, s, k1, l1max);
        for (l1 = l1max-1; l1 >= l1min; --l1)
            step_l1_down(n, s, k1, l1);
    }
}

/* Calculate the inner product of two sets of isoscalar factors */
sqrat isoscalar_context::inner_product(long m, long n)
{
    long k1, l1, k2, l2;
    sqrat result = sqrat(0);

    for (k1 = q1; k1 <= p1+q1; ++k1)
        for (l1 = 0; l1 <= q1; ++l1)
            for (k2 = q2; k2 <= p2+q2; ++k2)
            {
                /* Some states cannot couple due to hypercharge conservation;
                    ignore those */
                l2 = A - (k1+l1+k2);
                if ((l2 < 0) || (l2 > q2)) continue;

                result += isf(m, p+q, 0, k1, l1, k2, l2)
                        * isf(n, p+q, 0, k1, l1, k2, l2);
            }

    return result;
}

/* Calculate couplings to the state of highest weight.
    This can throw std::logic_error if we can't calculate directly. This should
    never happen, however, as isoscalars() has logic to avoid those cases.
*/
void isoscalar_context::calc_shw()
{
    long smax = min(A, (2*q1 + 2*q2 + 4*p1 + 4*p2 + q - p)/3);
    long smin = max(p + q, abs(2*q1 + 2*q2 - A));

    long k1min, k1max, l1min, l1max;
    long k1, l1, k2, l2;

    /* Fill out the topmost d planes for each of the degenerate reps */
    long m, n, s;
    for (m = 0; m < d; ++m)
    {
        s = smax - 2*m;
        k1min = max(q1, (A + s)/2 - (p2+q2));
        k1max = min(p1+q1, (A + s)/2 - q2);
        l1min = max(0, (A - s)/2 - q2);
        l1max = min(q1, (A - s)/2);

        /* Set one ISF in one particular irrep (leaving the same ISF
            in the other irreps as zero) */
        set_isf(m, p+q, 0, k1min, l1min, (A+s)/2 - k1min, (A-s)/2 - l1min, 1);

        for (n = 0; n < d; ++n)
        {
            /* Use recursion relations (possibly involving the plane
                above the current one, which will already have been filled)
                to fill out the rest of this plane */

            /* First fill across, from (k1min, l1min) to (k1min, l1max) */
            for (l1 = l1min+1; l1 <= l1max; ++l1)
                step_l1_up(n, s, k1min, l1);

            /* Now step upwards through the rows, from k1min to k1max */
            for (k1 = k1min+1; k1 <= k1max; ++k1)
            {
                step_k1_up(n, s, k1, l1min);
                for (l1 = l1min+1; l1 <= l1max; ++l1)
                    step_l1_up(n, s, k1, l1);
            }
        }
    }

    /* Now we have filled out the topmost d planes, step down
        through the rest of them */
    for (s = smax - 2*d; s >= smin; s -= 2)
        for (n = 0; n < d; ++n)
            step_s_down(n, s);

    /* Orthonormalise the ISFs for different representations.
        We orthogonalise each rep against *later* reps in order to get equivalent
        results to the algorithm described in our references. */
    for (n = d-1; n >= 0; --n)
    {
        sqrat v;
        for (m = n+1; m < d; ++m)
        {
            /* Factor to multiply rep 'm' by when subtracting from rep 'n'.
                Note that, when we get here, rep 'm' is already normalised.
            */
            v = inner_product(m, n);

            for (k1 = q1; k1 <= p1+q1; ++k1)
                for (l1 = 0; l1 <= q1; ++l1)
                    for (k2 = q2; k2 <= p2+q2; ++k2)
                    {
                        /* Some states cannot couple due to hypercharge conservation;
                            ignore those */
                        l2 = A - (k1+l1+k2);
                        if ((l2 < 0) || (l2 > q2)) continue;

                        set_isf(n, p+q, 0, k1, l1, k2, l2,
                            isf(n, p+q, 0, k1, l1, k2, l2)
                            - v * isf(m, p+q, 0, k1, l1, k2, l2));
                    }
        }

        /* Normalisation.
            Note that, at the same time, we enforce the sign convention that
            F(p+q, 0; p1+q1, 0, k2max, l2min) > 0.
            Here k2max means "The highest k2 which couples the state
            (p1+q1, 0) in rep 1 to (p+q, 0) in the target rep".
            This works out to be determined by the following:
        */
        long B = (-p1 + 2*p2 + q1 + 4*q2 + p - q)/3;
        long k2max = min(p2+q2, B);
        long l2min = max(0, B - p2 - q2);

        v = sqrt(inner_product(n, n));

        /* Step through until we find a state which couples */
        while (isf(n, p+q, 0, p1+q1, 0, k2max, l2min) == 0)
        {
            k2max -= 1;
            l2min += 1;
        }

        /* If this state has negative coupling, we need to negate the rep,
            in order to match the sign convention */
        if (isf(n, p+q, 0, p1+q1, 0, k2max, l2min) < 0)
            v = -v;

        for (k1 = q1; k1 <= p1+q1; ++k1)
            for (l1 = 0; l1 <= q1; ++l1)
                for (k2 = q2; k2 <= p2+q2; ++k2)
                {
                    /* Some states cannot couple due to hypercharge conservation;
                        ignore those */
                    l2 = A - (k1+l1+k2);
                    if ((l2 < 0) || (l2 > q2)) continue;

                    set_isf(n, p+q, 0, k1, l1, k2, l2,
                        isf(n, p+q, 0, k1, l1, k2, l2) / v);
                }
    }
}
