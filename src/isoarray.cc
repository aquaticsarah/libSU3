/* libSU3: Externally-visible container for isoscalar coefficients
    For some notes on how the indexing is done, see include/SU3.h
*/

#include <assert.h>

#include "SU3_internal.h"

/* Note: This type takes ownership of the array passed in - that is, it will
    delete the array when the isoarray object is deleted. */
isoarray::isoarray(long p, long q, long p1, long q1, long p2, long q2, long d,
    sqrat* coefficients) : p(p), q(q), p1(p1), q1(q1), p2(p2), q2(q2),
    d(d), coefficients(coefficients)
{}

isoarray::~isoarray()
{
    delete[] coefficients;
}

/* We use operator() instead of operator[] as an easy way to use
    multiple indices */
sqrat isoarray::operator()(long n, long k, long l, long k1, long l1,
                            long k2, long l2)
{
    /* Bounds checks */
    assert((n >= 0) && (n < d));
    assert((k >= q) && (k <= p+q));
    assert((l >= 0) && (l <= q));
    assert((k1 >= q1) && (k1 <= p1+q1));
    assert((l1 >= 0) && (l1 <= q1));
    assert((k2 >= q2) && (k2 <= p2+q2));
    assert((l2 >= 0) && (l2 <= q2));

    /* Check hypercharge conservation */
    assert(k1+l1+k2+l2-k-l == (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3);

    /* If compiled with -DNDEBUG, this line is to remove a compiler warning
        about l2 being unused */
    (void)l2;

    size_t index = ((((n * (p+1) + k-q) * (q+1) + l) * (p1+1) + k1-q1)
                    * (q1+1) + l1) * (p2+1) + k2-q2;
    return coefficients[index];
}
