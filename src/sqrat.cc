/* libSU3: Type representing (+ or -) the square root of a rational */

#include <stdio.h>

#include "SU3.h"

/* Calculate a GCD, giving the result the same sign as q.
    This is so that, when we reduce a fraction, the numerator can have any
    sign, but the denominator is always >= 0 */
long gcd(long p, long q)
{
    if (p == 0) return q;
    if (q == 0) return p;

    if (p < 0) p = -p;

    /* Make sure to give the result the sign of q */
    if (q < 0) return -gcd(-q, p%(-q));
    else return gcd(q, p%q);
}

/* Reduce a fraction to its lowest terms */
void sqrat::reduce()
{
    long d = gcd(p, q);
    p /= d;
    q /= d;
}

/* Various constructors */
sqrat::sqrat(long num, long denom) : p(num), q(denom)
{
    reduce();
}

/* Note that when initialising with an integer,
   we need to square it because all our fractions are
   implicitly under a square root sign */
sqrat::sqrat(long i) : q(1)
{
    if (i < 0) p = -(i*i);
    else p = i*i;
}

sqrat::sqrat() : p(0), q(1) {}

/* Temporary functions to get components */
long sqrat::numerator()
{
    return p;
}

long sqrat::denominator()
{
    return q;
}

/* Arithmetic. Note that for multiplication by scalars,
   we need to square the other term because all our fractions are
   implicitly under a square root sign */
sqrat sqrat::operator*(sqrat other)
{
    return sqrat(p * other.p, q * other.q);
}

sqrat sqrat::operator*(long other)
{
    return *this * sqrat(other);
}

sqrat sqrat::operator/(sqrat other)
{
    return sqrat(p * other.q, q * other.p);
}

sqrat sqrat::operator/(long other)
{
    return *this / sqrat(other);
}

/* Convert to a string */
char* sqrat::tostring(char* buffer, size_t len)
{
    if (p < 0)
        snprintf(buffer, len, "-sqrt(%ld/%ld)", -p, q);
    else
        snprintf(buffer, len, "sqrt(%ld/%ld)", p, q);
    return buffer;
}
