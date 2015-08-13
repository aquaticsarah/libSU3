/* libSU3: Coefficients for the various recursion relations we use

   Note: The functions step_s_down_* often make steps "from" locations
   which do not correspond to valid states. In those cases, we know that
   the corresponding coefficient will always be ignored. However, in those
   cases, the coefficient is sometimes indeterminate or infinite, and this
   causes problems. Specifically, this occurs if one of k1-l1 or k2-l2 is
   zero.

   Because the coefficient is ignored in those cases, we can safely
   special-case it to be zero instead. */

#include "SU3_internal.h"

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
