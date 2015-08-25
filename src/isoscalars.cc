/* libSU3: State-of-higest-weight (SHW) determination.
   Note: This also does degeneracy resolution */

#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>

#include "SU3.h"

/* A class for storing a bunch of useful values during our calculations.
    All functions are run as methods of an object of this class, so we have
    easy access to those values */
class isoscalar_context
{
    /* Target and factor reps */
    long p, q, p1, q1, p2, q2;

    /* Degeneracy */
    long d;

    isoarray* isf;

    void a_coefficients(long k1, long l1, long k2, long l2,
                        sqrat& a1, sqrat& a2, sqrat& a3, sqrat& a4);
    void b_coefficients(long k1, long l1, long k2, long l2,
                        sqrat& b1, sqrat& b2, sqrat& b3, sqrat& b4);

    void step_k1_up(long n, long k1, long l1, long k2, long l2);
    void step_k1_down(long n, long k1, long l1, long k2, long l2);
    void step_l1_up(long n, long k1, long l1, long k2, long l2);
    void step_l1_down(long n, long k1, long l1, long k2, long l2);

    int step_s_down_topleft(long n, long s,
            long k1min, long k1max, long l1min, long l1max);
    int step_s_down_bottomright(long n, long s,
            long k1min, long k1max, long l1min, long l1max);

protected:
    isoscalar_context(isoarray* isf, long p, long q, long p1,
                    long q1, long p2, long q2);
    int calc_shw();

    /* Allow the top-level driver function to interact with
        objects of this class */
    friend void isoscalars(long p, long q, long p1, long q1, long p2, long q2);
};

isoscalar_context::isoscalar_context(isoarray* isf, long p, long q, long p1,
            long q1, long p2, long q2) : p(p), q(q), p1(p1),
            q1(q1), p2(p2), q2(q2), isf(isf)
{
    d = degeneracy(p, q, p1, q1, p2, q2);
}

/* Functions to calculate the coefficients for each of the four recursion
   relations. Each stores the coefficients in its last four arguments.
   The expressions are based on arXiv:nucl-th/9511025.

   Note: The function step_s_down often makes steps "from" locations which
   do not correspond to valid states. In those cases, we know that the
   corresponding coefficient will always be ignored. However, in those cases,
   the coefficient is sometimes indeterminate or infinite, and this causes
   problems. Specifically, this occurs if one of k1-l1 or k2-l2 is zero.

   Because the coefficient is ignored in those cases, we can safely
   special-case it to be zero instead.
*/
void isoscalar_context::a_coefficients(long k1, long l1, long k2, long l2,
                                    sqrat& a1, sqrat& a2, sqrat& a3, sqrat& a4)
{
    long numerator, denominator;
    long s = k1 - l1 + k2 - l2; /* = 2(I_1 + I_2) */
    long t = k1 - l1 - k2 + l2; /* = 2(I_1 - I_2) */

    if (k1 == l1)
        a1 = sqrat(0);
    else
    {
        numerator = 2 * (k1+1) * (k1-q1) * (p1+q1-k1+1) * (p+q+s+3) * (p+q+t+1);
        denominator = (k1-l1) * (k1-l1+1);
        a1 = sqrat(numerator, denominator);
    }

    if (k2 == l2)
    {
        a2 = sqrat(0);
    }
    else
    {
        numerator = 2 * (k2+1) * (k2-q2) * (p2+q2-k2+1) * (p+q+s+3) * (p+q-t+1);
        denominator = (k2-l2) * (k2-l2+1);
        a2 = sqrat(numerator, denominator);
    }

    numerator = -2 * l1 * (q1-l1+1) * (p1+q1-l1+2) * (-p-q+s+1) * (p+q-t+1);
    denominator = (k1-l1+1) * (k1-l1+2);
    a3 = sqrat(numerator, denominator);

    numerator = 2 * l2 * (q2-l2+1) * (p2+q2-l2+2) * (-p-q+s+1) * (p+q+t+1);
    denominator = (k2-l2+1) * (k2-l2+2);
    a4 = sqrat(numerator, denominator);
}

void isoscalar_context::b_coefficients(long k1, long l1, long k2, long l2,
                                    sqrat& b1, sqrat& b2, sqrat& b3, sqrat& b4)
{
    long numerator, denominator;
    long s = k1 - l1 + k2 - l2; /* = 2(I_1 + I_2) */
    long t = k1 - l1 - k2 + l2; /* = 2(I_1 - I_2) */

    numerator = 2 * (k1+2) * (k1-q1+1) * (p1+q1-k1) * (-p-q+s+1) * (p+q-t+1);
    denominator = (k1-l1+1) * (k1-l1+2);
    b1 = sqrat(numerator, denominator);

    numerator = -2 * (k2+2) * (k2-q2+1) * (p2+q2-k2) * (-p-q+s+1) * (p+q+t+1);
    denominator = (k2-l2+1) * (k2-l2+2);
    b2 = sqrat(numerator, denominator);

    if (k1 == l1)
        b3 = sqrat(0);
    else
    {
        numerator = 2 * (l1+1) * (q1-l1) * (p1+q1-l1+1) * (p+q+s+3) * (p+q+t+1);
        denominator = (k1-l1) * (k1-l1+1);
        b3 = sqrat(numerator, denominator);
    }

    if (k2 == l2)
        b4 = sqrat(0);
    else
    {
        numerator = 2 * (l2+1) * (q2-l2) * (p2+q2-l2+1) * (p+q+s+3) * (p+q-t+1);
        denominator = (k2-l2) * (k2-l2+1);
        b4 = sqrat(numerator, denominator);
    }
}

/* Use the A and B recursion relations to step along the
   k1 and l1 axes within a plane, or to step down through planes.

   The values k1,l1,k2,l2 correspond to the value that is actually calculated.
   However, the recurrence relations are given relative to a different
   location. As such, we have to readjust one coordinate in each function.
   But this detail is hidden from the caller.
*/
void isoscalar_context::step_k1_up(long n, long k1, long l1, long k2, long l2)
{
    /* Calculate coefficients for the A recursion relation */
    sqrat a1, a2, a3, a4;
    a_coefficients(k1, l1, k2+1, l2, a1, a2, a3, a4);

    /* Calculate the value at (k1,l1,k2-1,l2) using surrounding values */
    sqrat res = (-a1 * (*isf)(n, p+q, 0, k1-1, l1, k2+1, l2)
                 -a3 * (*isf)(n, p+q, 0, k1, l1-1, k2+1, l2)
                 -a4 * (*isf)(n, p+q, 0, k1, l1, k2+1, l2-1)) / a2;
    (*isf)(n, p+q, 0, k1, l1, k2, l2) = res;
}

void isoscalar_context::step_k1_down(long n, long k1, long l1, long k2, long l2)
{
    /* Calculate coefficients for the A recursion relation */
    sqrat a1, a2, a3, a4;
    a_coefficients(k1+1, l1, k2, l2, a1, a2, a3, a4);

    /* Calculate the value at (k1-1,l1,k2,l2) using surrounding values */
    sqrat res = (-a2 * (*isf)(n, p+q, 0, k1+1, l1, k2-1, l2)
                 -a3 * (*isf)(n, p+q, 0, k1+1, l1-1, k2, l2)
                 -a4 * (*isf)(n, p+q, 0, k1+1, l1, k2, l2-1)) / a1;
    (*isf)(n, p+q, 0, k1, l1, k2, l2) = res;
}

void isoscalar_context::step_l1_up(long n, long k1, long l1, long k2, long l2)
{
    /* Calculate coefficients for the B recursion relation */
    sqrat b1, b2, b3, b4;
    b_coefficients(k1, l1-1, k2, l2, b1, b2, b3, b4);

    /* Calculate the value at (k1,l1,k2,l2) using surrounding values */
    sqrat res = (-b1 * (*isf)(n, p+q, 0, k1+1, l1-1, k2, l2)
                 -b2 * (*isf)(n, p+q, 0, k1, l1-1, k2+1, l2)
                 -b4 * (*isf)(n, p+q, 0, k1, l1-1, k2, l2+1)) / b3;
    (*isf)(n, p+q, 0, k1, l1, k2, l2) = res;
}

void isoscalar_context::step_l1_down(long n, long k1, long l1, long k2, long l2)
{
    /* Calculate coefficients for the B recursion relation */
    sqrat b1, b2, b3, b4;
    b_coefficients(k1, l1, k2, l2-1, b1, b2, b3, b4);

    /* Calculate the value at (k1,l1,k2,l2) using surrounding values */
    sqrat res = (-b1 * (*isf)(n, p+q, 0, k1+1, l1, k2, l2-1)
                 -b2 * (*isf)(n, p+q, 0, k1, l1, k2+1, l2-1)
                 -b3 * (*isf)(n, p+q, 0, k1, l1+1, k2, l2-1)) / b4;
    (*isf)(n, p+q, 0, k1, l1, k2, l2) = res;
}

/* In order to step from one plane (at s+2) down to the next plane (at s),
   there are two independent choices we can make:
   * We can try recursion relation A or recursion relation B
   * We can try to calculate the value at (k1max,l1min) or at (k1min,l1max)
   Each of these works in different cases, so we try all four combinations.
   Further, each option in the second choice requires *different* code,
   so we have one function for each (and just try them in turn).
   Each returns 1 on success, 0 on failure.

   In terms of our functions above, we step "from" a non-existent state
   (with k1,l1 not in the valid range) to the state we want. This works
   because 'isoscalar_array' allows us to request (certain) non-existent states
   and just returns 0 for the coupling coefficient. This is exactly what we
   need for the stepping to work properly.
*/
int isoscalar_context::step_s_down_topleft(long n, long s,
        long k1min, long k1max, long l1min, long l1max)
{
    long A = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;

    /* Try to calculate the value at (k1min, l1max) */
    try
    {
        step_k1_up  (n, k1min, l1max, (A+s)/2 - k1min, (A-s)/2 - l1max);
    }
    catch (std::domain_error& e)
    {
        try
        {
            step_l1_down(n, k1min, l1max, (A+s)/2 - k1min, (A-s)/2 - l1max);
        }
        catch (std::domain_error& e)
        {
            return 0;
        }
    }

    /* If we get here, we succeded. So fill out the rest of this plane. */
    long k1, l1;
    for (l1 = l1max-1; l1 >= l1min; --l1)
        step_l1_down(n, k1min, l1, (A+s)/2 - k1min, (A-s)/2 - l1);

    for (k1 = k1min+1; k1 <= k1max; ++k1)
    {
        step_k1_up(n, k1, l1max, (A+s)/2 - k1, (A-s)/2 - l1max);
        for (l1 = l1max-1; l1 >= l1min; --l1)
            step_l1_down(n, k1, l1, (A+s)/2 - k1, (A-s)/2 - l1);
    }

    return 1;
}

int isoscalar_context::step_s_down_bottomright(long n, long s,
        long k1min, long k1max, long l1min, long l1max)
{
    long A = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;

    /* Try to calculate the value at (k1max, l1min) */
    try
    {
        step_k1_down(n, k1max, l1min, (A+s)/2 - k1max, (A-s)/2 - l1min);
    }
    catch (std::domain_error& e)
    {
        try
        {
            step_l1_up  (n, k1max, l1min, (A+s)/2 - k1max, (A-s)/2 - l1min);
        }
        catch (std::domain_error& e)
        {
            return 0;
        }
    }

    /* If we get here, we succeded. So fill out the rest of this plane. */
    long k1, l1;
    for (l1 = l1min+1; l1 <= l1max; ++l1)
        step_l1_up(n, k1max, l1, (A+s)/2 - k1max, (A-s)/2 - l1);

    for (k1 = k1max-1; k1 >= k1min; --k1)
    {
        step_k1_down(n, k1, l1min, (A+s)/2 - k1, (A-s)/2 - l1min);
        for (l1 = l1min+1; l1 <= l1max; ++l1)
            step_l1_up(n, k1, l1, (A+s)/2 - k1, (A-s)/2 - l1);
    }

    return 1;
}

/* Calculate couplings to the state of highest weight.
    Returns 1 on success, 0 if we need to try again with
    the reps conjugated. */
int isoscalar_context::calc_shw()
{
    long A =           (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;
    long smax = min(A, (2*q1 + 2*q2 + 4*p1 + 4*p2 + q - p)/3);
    long smin = max(p + q, 2*q1 + 2*q2 - A);

    long k1min, k1max, l1min, l1max;

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
        (*isf)(m, p+q, 0, k1min, l1min, (A+s)/2 - k1min, (A-s)/2 - l1min) = 1;

        for (n = 0; n < d; ++n)
        {
            /* Use recursion relations (possibly involving the plane
               above the current one, which will already have been filled)
               to fill out the rest of this plane */

            /* First fill across, from (k1min, l1min) to (k1min, l1max) */
            long k1, l1;
            for (l1 = l1min+1; l1 <= l1max; ++l1)
                step_l1_up(n, k1min, l1, (A+s)/2 - k1min, (A-s)/2 - l1);

            /* Now step upwards through the rows, from k1min to k1max */
            for (k1 = k1min+1; k1 <= k1max; ++k1)
            {
                step_k1_up(n, k1, l1min, (A+s)/2 - k1, (A-s)/2 - l1min);
                for (l1 = l1min+1; l1 <= l1max; ++l1)
                    step_l1_up(n, k1, l1, (A+s)/2 - k1, (A-s)/2 - l1);
            }
        }
    }

    /* Now we have filled out the topmost d planes, step down
       through the rest of them */
    for (s = smax - 2*d; s >= smin; s -= 2)
    {
        k1min = max(q1, (A + s)/2 - (p2+q2));
        k1max = min(p1+q1, (A + s)/2 - q2);
        l1min = max(0, (A - s)/2 - q2);
        l1max = min(q1, (A - s)/2);

        fprintf(stderr,
            "(%ld,%ld)x(%ld,%ld) -> (%ld,%ld):\n"
            "Stepping down to s=%ld; k1min=%ld, k1max=%ld, l1min=%ld, l1max=%ld\n\n",
            p1, q1, p2, q2, p, q, s, k1min, k1max, l1min, l1max);

        for (n = 0; n < d; ++n)
        {
            /* Try both options */
            if (   (! step_s_down_topleft(n, s, k1min, k1max, l1min, l1max))
                && (! step_s_down_bottomright(n, s, k1min, k1max, l1min, l1max)))
            {
                /* If we get here, the stepdown algorithm failed, and we need to
                   try again with the reps conjugated */
                return 0;
            }
        }
    }

    /* Orthonormalise */
    /* TODO */

    return 1;
}

void isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (! d) return; /* Ignore reps of zero degeneracy */

    isoarray* isf = new isoarray(p,q,p1,q1,p2,q2);
    isoscalar_context* ctx = new isoscalar_context(isf,p,q,p1,q1,p2,q2);

    ctx->calc_shw();

    /* TODO: Calculate the rest of the ISFs, given those for the SHW */

    delete ctx;
    delete isf;
}
