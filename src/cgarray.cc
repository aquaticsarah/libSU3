/* libSU3: Externally-visible container for Clebsch-Gordan coefficients.
    This works by storing the corresponding isoscalar coefficients, then
    multiplying them on-demand by appropriate SU(2) Clebsch-Gordans
*/

#include <assert.h>

#include "SU3_internal.h"

/* Note: This type takes ownership of the isoarray object passed in */
cgarray::cgarray(isoarray* isf) : isf(isf) {}

cgarray::~cgarray()
{
    delete isf;
}

/* We use operator() instead of operator[] as an easy way to use
    multiple indices */
sqrat cgarray::operator()(long n, long k, long l, long m,
                            long k1, long l1, long m1,
                            long k2, long l2, long m2)
{
    /* Bounds checks */
    assert((m >= l) && (m <= k));
    assert((m1 >= l1) && (m1 <= k1));
    assert((m2 >= l2) && (m2 <= k2));

    mpq_class I  = mpq_class(k -l , 2), Iz  = mpq_class(m -l , 1) - I;
    mpq_class i1 = mpq_class(k1-l1, 2), i1z = mpq_class(m1-l1, 1) - i1;
    mpq_class i2 = mpq_class(k2-l2, 2), i2z = mpq_class(m2-l2, 1) - i2;
    sqrat su2_cg = su2_cgc(I, Iz, i1, i1z, i2, i2z);

    return (*isf)(n, k, l, k1, l1, k2, l2) * su2_cg;
}

/* Convert to ISFs. This returns a newly-allocated isoarray object. */
isoarray* cgarray::to_isoarray()
{
    /* Create a new coefficient array, to be owned by the isoarray we will
        end up returning */
    sqrat* new_isf_array = new sqrat[isf->size];

    /* Explicitly copy each element of the array */
    size_t i;
    for (i = 0; i < isf->size; ++i)
        new_isf_array[i] = isf->isf_array[i];

    return new isoarray(isf->p, isf->q, isf->p1, isf->q1,
                            isf->p2, isf->q2, isf->d, new_isf_array);
}
