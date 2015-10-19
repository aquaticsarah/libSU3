/* libSU3: Type representing (+ or -) the square root of a rational */

#include <stdio.h>
#include <math.h>
#include <stdexcept>

#include "SU3_internal.h"

/* Helper: Calculate the square root of an integer, which must be a square.
    Raises an exception if the values is not a square */
static mpq_class sqrt(mpq_class x)
{
    if (x < 0) throw std::domain_error("Value is not a square");

    /* Try to square-root the numerator and denominator independently.
        Note that mpz_root() returns 0 iff the value passed into it
        is not a perfect square.
    */
    mpz_t num, denom;
    mpz_inits(num, denom, NULL);

    if ((! mpz_root(num, x.get_num_mpz_t(), 2))
     || (! mpz_root(denom, x.get_den_mpz_t(), 2)))
    {
        mpz_clears(num, denom, NULL);
        throw std::domain_error("Value is not a square");
    }

    mpq_class res = mpq_class(mpz_class(num), mpz_class(denom));
    mpz_clears(num, denom, NULL);
    return res;
}

/* Component-wise constructors. sqrat(p, q) returns sign(pq) * sqrt(|p|/|q|) */
sqrat::sqrat(mpz_class num, mpz_class denom) : v(num, denom)
{
    v.canonicalize();
}

sqrat::sqrat(long num, long denom) : v(num, denom)
{
    v.canonicalize();
}

/* Casts from other types. sqrat(x) returns a value which should equal x. */
sqrat::sqrat(mpq_class v) : v(v*abs(v))
{
    v.canonicalize();
}

sqrat::sqrat(long v) : v(v*abs(v)) {}
sqrat::sqrat() : v(0) {}

/* Internal constructor: Build a value which is sign(v) * sqrt(|v|) */
static sqrat sqrat_raw(mpq_class v)
{
    return sqrat(v.get_num(), v.get_den());
}

/* Unary operators */
sqrat operator+(const sqrat& value)
{
    return sqrat_raw(value.v);
}

sqrat operator-(const sqrat& value)
{
    return sqrat_raw(-value.v);
}

/* Arithmetic */
sqrat operator*(sqrat left, const sqrat& right)
{
    left *= right;
    return left;
}

sqrat& sqrat::operator*=(const sqrat& other)
{
    v *= other.v;
    return *this;
}

sqrat operator/(sqrat left, const sqrat& right)
{
    left /= right;
    return left;
}

sqrat& sqrat::operator/=(const sqrat& other)
{
    /* Special case: If dividing by 0, we need to raise an exception,
        instead of actually doing the division (which would cause a
        SIGFPE, which we can't catch easily)
    */
    if (other.v == 0)
        throw std::domain_error("sqrat: division by zero");

    v /= other.v;
    return *this;
}

/* Addition and subtraction are a bit more complicated.
    The general form is:
    +- sqrt(|x|) = (+- sqrt(|v|)) +- (+- sqrt(|w|))
    where each of the four +- signs is independent.

    By playing around with overall signs, we can fix v,w to be positive.
    Then we just have two cases:
    sqrt(|x|) = sqrt(v) +- sqrt(w)
 => |x| = v +- 2 sqrt(vw) + w

    Now we have to consider the sign of x. As v,w >= 0, the only way
    x can be negative is if we are subtracting, and if v < w.

    We have a helper function to deal with these reduced cases, then the
    operator+ and operator- methods just need to deal with signs.
    The 'sign' argument chooses between + (if 1) and - (if 0)
*/
static mpq_class add_internal(mpq_class v, mpq_class w, int sign)
{
    if (sign)       return v + 2*sqrt(mpq_class(v*w)) + w;
    else if (v < w) return -v + 2*sqrt(mpq_class(v*w)) - w;
    else            return v - 2*sqrt(mpq_class(v*w)) + w;
}

sqrat operator+(sqrat left, const sqrat& right)
{
    left += right;
    return left;
}

sqrat& sqrat::operator+=(const sqrat& other)
{
    if (v >= 0)
    {
        if (other.v >= 0)
            v = add_internal(v, other.v, 1);
        else
            v = add_internal(v, -other.v, 0);
    }
    else
    {
        if (other.v >= 0)
            v = -add_internal(-v, other.v, 0);
        else
            v = -add_internal(-v, -other.v, 1);
    }

    return *this;
}

sqrat operator-(sqrat left, const sqrat& right)
{
    left -= right;
    return left;
}

sqrat& sqrat::operator-=(const sqrat& other)
{
    if (v >= 0)
    {
        if (other.v >= 0)
            v = add_internal(v, other.v, 0);
        else
            v = add_internal(v, -other.v, 1);
    }
    else
    {
        if (other.v >= 0)
            v = -add_internal(-v, other.v, 1);
        else
            v = -add_internal(-v, -other.v, 0);
    }

    return *this;
}

sqrat sqrt(const sqrat& value)
{
    return sqrat_raw(sqrt(value.v));
}

/* Comparisons */
bool operator<(const sqrat& left, const sqrat& right)
{
    return left.v < right.v;
}

bool operator<=(const sqrat& left, const sqrat& right)
{
    return left.v <= right.v;
}

bool operator==(const sqrat& left, const sqrat& right)
{
    return left.v == right.v;
}

bool operator!=(const sqrat& left, const sqrat& right)
{
    return left.v != right.v;
}

bool operator>(const sqrat& left, const sqrat& right)
{
    return left.v > right.v;
}

bool operator>=(const sqrat& left, const sqrat& right)
{
    return left.v >= right.v;
}

/* Conversions to various types */
char* sqrat::tostring(char* buffer, size_t len)
{
    if (v < 0)
    {
        mpq_class w = -v;
        gmp_snprintf(buffer, len, "-sqrt(%Qd)", w.get_mpq_t());
    }
    else
        gmp_snprintf(buffer, len, "sqrt(%Qd)", v.get_mpq_t());
    return buffer;
}

sqrat::operator double()
{
    double x = v.get_d();

    /* The value we want is sign(x) * sqrt(|x|),
        so correct for that */
    if (std::signbit(x))
        return -sqrt(-x);
    else
        return sqrt(x);
}
