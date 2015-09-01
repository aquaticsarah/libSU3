/* libSU3: External interface */

#ifndef __SU3_H__
#define __SU3_H__

#include <stddef.h>
#include <gmpxx.h>

/* Class to represent (+ or -) the square root of a rational.
    The actual value represented is sign(v) * |v|
*/
class sqrat
{
public:
    mpq_class v;

    sqrat(mpq_class v);
    sqrat(long, long);
    sqrat(long);
    sqrat();

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

/* Classes to hold SU(3) isoscalar factors and Clebsch-Gordan coefficients.

   Notes on indexing:
   - For each of the three reps involved, we have the following ranges:
      q <= k <= p+q (for a total of p+1 possible values of k)
      0 <= l <= q   (for a total of q+1 possible values of l)
     The value of 'l' can be used as an index directly, but that for
     'k' needs to be shifted down by q.

   - All reps in a degenerate set need to be processed at once, so we
     allocate storage all at once.

   Implementation details:
   - Sometimes, when calculating ISFs, we apply a recursion relation involving
     values which are one space "off the edge" of the valid range (eg, with
     l1=-1). In order to simplify this code, we allow indexes of that form,
     but *only* for the internal-only isoscalar_context class.

   - Given k, l, k1, l1, k2, there is a unique valid value of l2 determined
     by hypercharge conservation. As such, we can save a factor of (q2+1) on
     memory space by not allocating an l2 axis.
     However, we do accept l2 as an argument and, unless -DNDEBUG is specified
     when compiling the library, we check that it is valid.
*/
class isoarray;
class cgarray;

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
