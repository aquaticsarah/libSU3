/* libSU3: Type representing (+ or -) the square root of a rational */

#include <stdio.h>
#include <stdexcept>

#include "SU3.h"

/* Helper: Calculate sign(x) * x^2 */
static long sign_square(long x)
{
    if (x < 0) return -x*x;
    else return x*x;
}

/* Helper: Calculate the square root of an integer, which must be a square.
    Raises an exception if the values is not a square */
static long isqrt(long x)
{
    if (x < 0) return -1;
    /* TODO: Possibly optimise this? */
    long v = 0;
    while (1)
    {
        if (v*v > x) throw std::domain_error("Value is not a square");
        else if (v*v == x) return v;
        ++v;
    }
}

/* Calculate a GCD, giving the result the same sign as q.
    This is so that, when we reduce a fraction, the numerator can have any
    sign, but the denominator is always >= 0 */
static long gcd(long p, long q)
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
    /* Prevent division by zero */
    if (denom == 0)
        throw std::domain_error("Division by zero when constructing sqrat");

    reduce();
}

/* Note that when initialising with an integer,
   we need to square it because all our fractions are
   implicitly under a square root sign */
sqrat::sqrat(long i) : p(sign_square(i)), q(1) {}

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

/* Unary plus and minus */
sqrat sqrat::operator+()
{
    return sqrat(p, q);
}

sqrat sqrat::operator-()
{
    return sqrat(-p, q);
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

/* Addition and subtraction are a bit more complicated.
    The general form is:
    +- sqrt(p/q) = (+- sqrt(r/s)) +- (+- sqrt(t/u))
    where each of the four +- signs is independent.

    By playing around with overall signs, we can reduce to two cases:
    (q,r,s,t,u all >=0, but the sign of p is to be determined later)
    sqrt(|p|/q) = sqrt(r/s) +- sqrt(t/u)
 => |p|/q = r/s +- 2 sqrt(rt/su) + t/u
 => |p|/q = 1/su (ru +- 2 sqrt(rtsu) + st)
    so we can set |p| = ru +- 2 sqrt(rtsu) + st, q = su and then reduce.

    Now we have to consider the sign of p. As r,s,t,u >= 0, the only way
    p can be negative is if we are subtracting. In that case, it happens
    iff r*u < s*t.

    We have a helper function to deal with these reduced cases, then the
    operator+ and operator- methods just need to deal with signs.
    The 'sign' argument chooses between + (if 1) and - (if 0) */
static sqrat add_internal(sqrat left, sqrat right, int sign)
{
    long r,s,t,u;
    r = left.numerator();
    s = left.denominator();
    t = right.numerator();
    u = right.denominator();

    if (sign)           return sqrat(r*u + 2*isqrt(r*s*t*u) + s*t, s*u);
    else if (r*u < s*t) return sqrat(-r*u + 2*isqrt(r*s*t*u) - s*t, s*u);
    else                return sqrat(r*u - 2*isqrt(r*s*t*u) + s*t, s*u);
}

/* Addition and subtraction */
sqrat sqrat::operator+(sqrat other)
{
    if (p >= 0)
    {
        if (other.p >= 0)
            return add_internal(*this, other, 1);
        else
            return add_internal(*this, -other, 0);
    }
    else
    {
        if (other.p >= 0)
            return -add_internal(-*this, other, 0);
        else
            return -add_internal(-*this, -other, 1);
    }
}

sqrat sqrat::operator+(long other)
{
    return *this + sqrat(other);
}

sqrat sqrat::operator-(sqrat other)
{
    if (p >= 0)
    {
        if (other.p >= 0)
            return add_internal(*this, other, 0);
        else
            return add_internal(*this, -other, 1);
    }
    else
    {
        if (other.p >= 0)
            return -add_internal(-*this, other, 1);
        else
            return -add_internal(-*this, -other, 0);
    }
}

sqrat sqrat::operator-(long other)
{
    return *this - sqrat(other);
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
