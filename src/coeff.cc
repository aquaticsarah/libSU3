/* libSU3: Coefficients for the various recursion relations we use

    Note: The function step_s_down often makes steps "from" locations
    which do not correspond to valid states. In those cases, we know that
    the corresponding coefficient will always be ignored. However, when doing
    this, the coefficient is sometimes indeterminate or infinite, which
    throws an exception.

    Because the coefficient is ignored in those cases, we can safely
    special-case it to be zero instead.

    Aside from this, the expressions are taken from Kaeding and Williams.
*/

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

void isoscalar_context::c_coefficients(long k, long l, long k1, long l1,
                    long k2, long l2, sqrat& alpha, sqrat& c1, sqrat& c2,
                    sqrat& c3, sqrat& c4)
{
    long numerator, denominator;
    long s = k1 - l1 + k2 - l2; /* = 2(I_1 + I_2) */
    long t = k1 - l1 - k2 + l2; /* = 2(I_1 - I_2) */

    numerator = (k-l+2)*(k-l+2);
    denominator = l*(q-l+1)*(p+q-l+2);
    alpha = sqrat(numerator, denominator);

    if (k-l+t+2 == 0)
    {
        c1 = sqrat(0);
        c2 = sqrat(0);
        c3 = sqrat(0);
    }
    else
    {
        numerator = (k+2)*(k-q+1)*(p+q-k)*(s-k+l)*(k-l-t+2);
        denominator = (k-l+2)*(k-l+2)*(k-l+s+4)*(k-l+t+2);
        c1 = sqrat(numerator, denominator);

        numerator = 4*l1*(q1-l1+1)*(p1+q1-l1+2)*(k1-l1+1);
        denominator = (k1-l1+2)*(k-l+s+4)*(k-l+t+2);
        c2 = sqrat(numerator, denominator);

        /* If k2==l2, c3 is infinite or indeterminate. But in that case,
            it is the coefficient of a state with l2>k2, which is impossible
            (ie, there must be zero coupling). Thus we can just replace it by 0.
        */
        if (k2 == l2)
            c3 = sqrat(0);
        else
        {
            numerator = -(k2+1)*(k2-q2)*(s-k+l)*(p2+q2-k2+1);
            denominator = (k2-l2)*(k2-l2+1)*(k-l+t+2);
            c3 = sqrat(numerator, denominator);
        }
    }

    numerator = l2*(q2-l2+1)*(p2+q2-l2+2)*(k-l-t+2);
    denominator = (k2-l2+1)*(k2-l2+2)*(k-l+s+4);
    c4 = sqrat(numerator, denominator);
}

void isoscalar_context::d_coefficients(long k, long k1, long l1,
                long k2, long l2, sqrat& beta, sqrat& d1, sqrat& d2, sqrat& d3)
{
    long numerator, denominator;
    long s = k1 - l1 + k2 - l2; /* = 2(I_1 + I_2) */
    long t = k1 - l1 - k2 + l2; /* = 2(I_1 - I_2) */

    numerator = k+2;
    denominator = (k-q+1)*(p+q-k);
    beta = sqrat(numerator, denominator);

    /* If k+t+2==0, then the state at which we are evaluating the recurrence
        relation is invalid (as it requires I=(I_2 - I_1) - 1, but in fact we
        have I >= I_2 - I_1).
        Hence we need to replace some coefficients by zero */
    if (k+t+2 == 0)
    {
        d1 = sqrat(0);
        d3 = sqrat(0);
    }
    else
    {
        numerator = 4*(k1+2)*(k1-q1+1)*(p1+q1-k1)*(k1-l1+1);
        denominator = (k1-l1+2)*(k+s+4)*(k+t+2);
        d1 = sqrat(numerator, denominator);

        /* If k2==l2, d3 is infinite or indeterminate. But in that case,
            it is the coefficient of a state with l2>k2, which is impossible
            (ie, there must be zero coupling). Thus we can just replace it by 0.
        */
        if (k2 == l2)
            d3 = sqrat(0);
        else
        {
            numerator = (l2+1)*(q2-l2)*(p2+q2-l2+1)*(s-k);
            denominator = (k2-l2)*(k2-l2+1)*(k+t+2);
            d3 = sqrat(numerator, denominator);
        }
    }

    numerator = (k2+2)*(k2-q2+1)*(p2+q2-k2)*(k-t+2);
    denominator = (k2-l2+1)*(k2-l2+2)*(k+s+4);
    d2 = sqrat(numerator, denominator);
}
