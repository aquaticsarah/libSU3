/* libSU3: External interface */

#ifndef __SU3_H__
#define __SU3_H__

#include <stddef.h>

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

    /* Arithmetic operations */
    sqrat operator+();
    sqrat operator-();

    sqrat operator*(sqrat);
    sqrat operator*(long);
    sqrat* operator*=(sqrat);
    sqrat* operator*=(long);

    sqrat operator/(sqrat);
    sqrat operator/(long);
    sqrat* operator/=(sqrat);
    sqrat* operator/=(long);

    sqrat operator+(sqrat);
    sqrat operator+(long);
    sqrat* operator+=(sqrat);
    sqrat* operator+=(long);

    sqrat operator-(sqrat);
    sqrat operator-(long);
    sqrat* operator-=(sqrat);
    sqrat* operator-=(long);

    friend sqrat sqrt(sqrat);

    /* Comparisons */
    bool operator<(sqrat);
    bool operator<(long);

    bool operator<=(sqrat);
    bool operator<=(long);

    bool operator==(sqrat);
    bool operator==(long);

    bool operator!=(sqrat);
    bool operator!=(long);

    bool operator>(sqrat);
    bool operator>(long);

    bool operator>=(sqrat);
    bool operator>=(long);

    /* Conversions to various types */
    char* tostring(char*, size_t);
    double todouble();
};

/* Reversed arithmetic operators */
sqrat operator*(long, sqrat);
sqrat operator/(long, sqrat);
sqrat operator+(long, sqrat);
sqrat operator-(long, sqrat);

/* Reversed comparisons */
bool operator<(long, sqrat);
bool operator<=(long, sqrat);
bool operator==(long, sqrat);
bool operator!=(long, sqrat);
bool operator>(long, sqrat);
bool operator>=(long, sqrat);

/* A class to hold the isoscalar factors for a particular coupling */
class isoarray
{
public:
    /* Target and factor reps */
    const long p, q, p1, q1, p2, q2;

    /* Degeneracy of target rep */
    const long d;

    const sqrat* coefficients;

    /* Note: This type takes ownership of the array passed in - that is, it will
        delete the array when the isoarray object is deleted. */
    isoarray(long p, long q, long p1, long q1, long p2, long q2, long d,
                sqrat* coefficients);
    ~isoarray();

    /* We use operator() instead of operator[] as an easy way to use
    multiple indices */
    sqrat operator()(long n, long k, long l, long k1, long l1,
                    long k2, long l2);
};

/* Information about representations */
long dimension(long p, long q);
char* repname(char* buffer, size_t len, long p, long q);
long degeneracy(long p, long q, long p1, long q1, long p2, long q2);

/* Main calculation function */
isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2);

#endif
