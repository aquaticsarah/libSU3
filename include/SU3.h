/* libSU3: External interface */

#ifndef __SU3_H__
#define __SU3_H__

#include <stddef.h>

/* Various useful functions */
long min(long a, long b);
long min(long a, long b, long c);
long min(long a, long b, long c, long d);
long max(long a, long b);
long max(long a, long b, long c, long d, long e, long f);
long gcd(long, long);

/* Class to represent (+ or -) the square root of a rational.
    The actual value represented is (sign(p) * sqrt(|p| / q))
*/
class sqrat
{
private:
    long p;
    long q;
    void reduce();

public:
    sqrat(long, long);
    sqrat(long);
    sqrat();

    long numerator();
    long denominator();

    sqrat operator+();
    sqrat operator-();

    /* Arithmetic operations */
    sqrat operator*(sqrat);
    sqrat operator*(long);

    sqrat operator/(sqrat);
    sqrat operator/(long);

    sqrat operator+(sqrat);
    sqrat operator+(long);

    sqrat operator-(sqrat);
    sqrat operator-(long);

    /* Convert to a string */
    char* tostring(char*, size_t);
};

/* A class to hold the isoscalar factors for a particular coupling */
class isoarray
{
public:
    /* Target and factor reps */
    long p, q, p1, q1, p2, q2;

    /* Degeneracy of target rep */
    long d;

    sqrat* coefficients;

    isoarray(long p, long q, long p1, long q1, long p2, long q2);
    ~isoarray();

    /* We use operator() instead of operator[] as an easy way to use
    multiple indices */
    sqrat& operator()(long n, long k, long l, long k1, long l1,
                        long k2, long l2);
};

/* Information about representations */
long dimension(long p, long q);
char* repname(char* buffer, size_t len, long p, long q);
long degeneracy(long p, long q, long p1, long q1, long p2, long q2);

/* Main calculation function */
void isoscalars(long p, long q, long p1, long q1, long p2, long q2);

#endif
