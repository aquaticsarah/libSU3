/* libSU3: Externally-visible container for isoscalar coefficients

    Notes on indexing:
    - All reps in a degenerate set need to be processed at once, so we
        allocate storage all at once.

    - For each of the three reps involved, we have the following ranges:
        q <= k <= p+q (for a total of p+1 possible values of k)
        0 <= l <= q   (for a total of q+1 possible values of l)
        The value of 'l' can be used as an index directly, but that for
        'k' needs to be shifted down by q.

    - Given k, l, k1, l1, k2, there is a unique valid value of l2 determined
        by hypercharge conservation. As such, we can save a factor of (q2+1) on
        memory space by not allocating an l2 axis.
        However, we do accept l2 as an argument and, unless -DNDEBUG is
        specified when compiling the library, we check that it is valid.
*/

#include <assert.h>

#include "SU3_internal.h"

/* Note: This type takes ownership of the array passed in - that is, it will
    delete the array when the isoarray object is deleted. */
isoarray::isoarray(long p, long q, long p1, long q1, long p2, long q2, long d,
    sqrat* isf_array) : isf_array(isf_array), p(p), q(q), p1(p1), q1(q1),
    p2(p2), q2(q2), d(d)
{
    size = d * (p+1) * (q+1) * (p1+1) * (q1+1) * (p2+1);
}

isoarray::~isoarray()
{
    delete[] isf_array;
}

void isoarray::set_isf(long n, long k, long l, long k1, long l1,
                        long k2, long l2, sqrat v)
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
    assert(index < size);
    isf_array[index] = v;
}

/* We use operator() instead of operator[] as an easy way to use
    multiple indices */
sqrat isoarray::operator()(long n, long k, long l, long k1, long l1,
                            long k2, long l2)
{
    /* Bounds checks - as this is a user-visible function, we don't want
        to crash if an invalid value is passed, just to return zero.
    */
    if (    (n < 0) || (n >= d)
         || (k  < q ) || (k  > p +q ) || (l  < 0) || (l  > q )
         || (k1 < q1) || (k1 > p1+q1) || (l1 < 0) || (l1 > q1)
         || (k2 < q2) || (k2 > p2+q2) || (l2 < 0) || (l2 > q2))
        return 0;

    /* Check hypercharge conservation */
    if(k1+l1+k2+l2-k-l != (2*p1 + 2*p2 + 4*q1 + 4*q2 - 2*p - 4*q)/3)
        return sqrat(0);

    size_t index = ((((n * (p+1) + k-q) * (q+1) + l) * (p1+1) + k1-q1)
                    * (q1+1) + l1) * (p2+1) + k2-q2;
    assert(index < size);
    return isf_array[index];
}

/* Convert to Clebsch-Gordans. This returns a newly-allocated cgarray object. */
cgarray* isoarray::to_cgarray()
{
    /* Create a new coefficient array, to be owned by the cgarray we will
        end up returning */
    sqrat* new_isf_array = new sqrat[size];

    /* Explicitly copy each element of the array */
    size_t i;
    for (i = 0; i < size; ++i)
        new_isf_array[i] = isf_array[i];

    isoarray* isf = new isoarray(p, q, p1, q1, p2, q2, d, new_isf_array);
    return new cgarray(isf);
}

/* Internal: Check that the sign convention is obeyed.
    Note: In the situations where this is called, we know that the sign
    is consistent between degenerate irreps, we just might have an overall
    - sign on everything.
*/
void isoarray::check_sign_convention()
{
    long B = (-p1 + 2*p2 + q1 + 4*q2 + p - q)/3;
    long k2max = min(p2+q2, B);
    long l2min = max(0, B - p2 - q2);

    while ((*this)(0, p+q, 0, p1+q1, 0, k2max, l2min) == 0)
    {
        k2max -= 1;
        l2min += 1;
    }

    /* Check the sign */
    if ((*this)(0, p+q, 0, p1+q1, 0, k2max, l2min) < 0)
    {
        /* If the sign is wrong, enforce the convention */
        long n, k, l, k1, l1, k2, l2;
        for (n = 0; n < d; ++n)
            FOREACH_ISF(p, q, p1, q1, p2, q2, k, l, k1, l1, k2, l2)
            {
                this->set_isf(n, k, l, k1, l1, k2, l2,
                    -(*this)(n, k, l, k1, l1, k2, l2));
            }
    }
}

/* Apply the various symmetry relations.
    The formulas for these relations are adapted from Williams.
*/
isoarray* isoarray::exch_12()
{
    size_t new_size = d * (p+1) * (q+1) * (p2+1) * (q2+1) * (p1+1);
    sqrat* new_isf_array = new sqrat[new_size];
    isoarray* array = new isoarray(p, q, p2, q2, p1, q1, d, new_isf_array);

    /* Fill the new array */
    long xi_1 = phase_exch_12(p, q, p1, q1, p2, q2);
    long n, k, l, k1, l1, k2, l2;

    for (n = 0; n < d; ++n)
        FOREACH_ISF(p, q, p2, q2, p1, q1, k, l, k2, l2, k1, l1)
            array->set_isf(n, k, l, k2, l2, k1, l1,
                SIGN((k-l-k1+l1-k2+l2)/2) * SIGN(n) * xi_1
                * (*this)(n, k, l, k1, l1, k2, l2));

    array->check_sign_convention();

    return array;
}

isoarray* isoarray::exch_13bar()
{
    /* For this symmetry, the array size is unchanged */
    sqrat* new_isf_array = new sqrat[size];
    isoarray* array = new isoarray(q1, p1, q, p, p2, q2, d, new_isf_array);

    /* Fill the new array */
    long n, k, l, k1, l1, k2, l2;

    for (n = 0; n < d; ++n)
        FOREACH_ISF(p, q, p1, q1, p2, q2, k, l, k1, l1, k2, l2)
            array->set_isf(n, p1+q1-l1, p1+q1-k1, p+q-l, p+q-k, k2, l2,
                    SIGN(l2+n)
                  * sqrat((p1+1)*(q1+1)*(p1+q1+2)*(k-l+1), (p+1)*(q+1)*(p+q+2)*(k1-l1+1))
                  * (*this)(n, k, l, k1, l1, k2, l2));

    array->check_sign_convention();

    return array;
}

/* Combination of the above two, for simplicity */
isoarray* isoarray::exch_23bar()
{
    isoarray* tmp1 = this->exch_12();
    isoarray* tmp2 = tmp1->exch_13bar();
    isoarray* result = tmp2->exch_12();

    delete tmp1;
    delete tmp2;
    return result;
}
