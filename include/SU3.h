/* libSU3: External interface */

#ifndef __SU3_H__
#define __SU3_H__

#include <stddef.h>

/* Class to represent (+ or -) the square root of a rational.
    The actual value represented is:
    if (p < 0)  then -sqrt(-p/q)
    if (p >= 0) then  sqrt(p/q)
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

/* Min and max for various numbers of arguments */
long min(long a, long b);
long min(long a, long b, long c);
long min(long a, long b, long c, long d);
long max(long a, long b);
long max(long a, long b, long c, long d, long e, long f);

/* Information about representations */
long dimension(long p, long q);
long degeneracy(long p1, long q1, long p2, long q2, long p, long q);

#endif
