/* libSU3: Type representing (+ or -) the square root of a rational */

#include <stdio.h>
#include <math.h>
#include <stdexcept>

#include "SU3_internal.h"

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
    if (x < 0) throw std::domain_error("Value is not a square");

    /* Use doubles to get an approximation to sqrt(x).
        We can guarantee that this approximation is within
        +-1 of the true value */
    long approx = lrint(sqrt((double)x));
    long approx_sq = approx * approx;

    if (approx_sq < x)
    {
        if ((approx+1)*(approx+1) == x)
            return approx+1;
        else throw std::domain_error("Value is not a square");
    }
    else if (approx_sq > x)
    {
        if ((approx-1)*(approx-1) == x)
            return approx-1;
        else throw std::domain_error("Value is not a square");
    }
    else return approx;
}

/* Reduce a fraction to its lowest terms */
void sqrat::reduce()
{
    long d = gcd(p, q);
    p /= d;
    q /= d;

    /* Make sure that q is nonnegative */
    if (q < 0)
    {
        p = -p;
        q = -q;
    }
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

/* Temporary functions to get components; TODO: Remove */
long sqrat::numerator()
{
    return p;
}

long sqrat::denominator()
{
    return q;
}

/* Unary operators */
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

sqrat* sqrat::operator*=(sqrat other)
{
    sqrat res = *this * other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat* sqrat::operator*=(long other)
{
    sqrat res = *this * other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat operator*(long left, sqrat right)
{
    return sqrat(left) * right;
}



sqrat sqrat::operator/(sqrat other)
{
    return sqrat(p * other.q, q * other.p);
}

sqrat sqrat::operator/(long other)
{
    return *this / sqrat(other);
}

sqrat* sqrat::operator/=(sqrat other)
{
    sqrat res = *this / other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat* sqrat::operator/=(long other)
{
    sqrat res = *this / other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat operator/(long left, sqrat right)
{
    return sqrat(left) / right;
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

sqrat* sqrat::operator+=(sqrat other)
{
    sqrat res = *this + other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat* sqrat::operator+=(long other)
{
    sqrat res = *this + other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat operator+(long left, sqrat right)
{
    return sqrat(left) + right;
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

sqrat* sqrat::operator-=(sqrat other)
{
    sqrat res = *this - other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat* sqrat::operator-=(long other)
{
    sqrat res = *this - other;
    p = res.p;
    q = res.q;
    return this;
}

sqrat operator-(long left, sqrat right)
{
    return sqrat(left) - right;
}



sqrat sqrt(sqrat v)
{
    return sqrat(isqrt(v.p), isqrt(v.q));
}

/* Conversions to various types */
char* sqrat::tostring(char* buffer, size_t len)
{
    if (p < 0)
        snprintf(buffer, len, "-sqrt(%ld/%ld)", -p, q);
    else
        snprintf(buffer, len, "sqrt(%ld/%ld)", p, q);
    return buffer;
}

double sqrat::todouble()
{
    return sqrt(p / (double)q);
}
