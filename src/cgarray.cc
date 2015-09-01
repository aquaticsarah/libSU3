/* libSU3: Externally-visible container for Clebsch-Gordan coefficients.
    This works by storing the corresponding isoscalar coefficients, then
    multiplying them on-demand by appropriate SU(2) Clebsch-Gordans
*/

#include <assert.h>

#include "SU3_internal.h"

/* Note: This type takes ownership of the array passed in - that is, it will
    delete the array when the isoarray object is deleted. */
cgarray::cgarray(long p, long q, long p1, long q1, long p2, long q2, long d,
    sqrat* isf_array) : isf_array(isf_array), p(p), q(q), p1(p1), q1(q1),
    p2(p2), q2(q2), d(d)
{
    size = d * (p+1) * (q+1) * (p1+1) * (q1+1) * (p2+1);
}

cgarray::~cgarray()
{
    delete[] isf_array;
}

/* We use operator() instead of operator[] as an easy way to use
    multiple indices */
sqrat cgarray::operator()(long n, long k, long l, long m,
                            long k1, long l1, long m1,
                            long k2, long l2, long m2)
{
    /* Bounds checks */
    assert((n >= 0) && (n < d));
    assert((k >= q) && (k <= p+q));
    assert((l >= 0) && (l <= q));
    assert((m >= l) && (m <= k));
    assert((k1 >= q1) && (k1 <= p1+q1));
    assert((l1 >= 0) && (l1 <= q1));
    assert((m1 >= l1) && (m1 <= k1));
    assert((k2 >= q2) && (k2 <= p2+q2));
    assert((l2 >= 0) && (l2 <= q2));
    assert((m2 >= l2) && (m2 <= k2));

    /* Check hypercharge conservation */
    assert(k1+l1+k2+l2-k-l == (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3);

    size_t index = ((((n * (p+1) + k-q) * (q+1) + l) * (p1+1) + k1-q1)
                    * (q1+1) + l1) * (p2+1) + k2-q2;
    sqrat isf = isf_array[index];

    mpq_class I  = mpq_class(k -l , 2), Iz  = mpq_class(m -l , 1) - I;
    mpq_class i1 = mpq_class(k1-l1, 2), i1z = mpq_class(m1-l1, 1) - i1;
    mpq_class i2 = mpq_class(k2-l2, 2), i2z = mpq_class(m2-l2, 1) - i2;
    sqrat su2_cg = su2_cgc(I, Iz, i1, i1z, i2, i2z);

    return isf * su2_cg;
}

/* Convert to ISFs */
isoarray* cgarray::to_isoarray()
{
    /* Create a new coefficient array, to be owned by the isoarray we will
        end up returning */
    sqrat* new_isf_array = new sqrat[size];

    /* Explicitly copy each element of the array */
    size_t i;
    for (i = 0; i < size; ++i)
        new_isf_array[i] = isf_array[i];

    isoarray* isf = new isoarray(p, q, p1, q1, p2, q2, d, new_isf_array);

    return isf;
}
