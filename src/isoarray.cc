/* libSU3: Custom type to store isoscalar coefficients */

#include <assert.h>

#include "SU3.h"

/* Notes on indexing:
   - For each of the three reps involved, we have the following ranges:
      q <= k <= p+q (for a total of p+1 possible values of k)
      0 <= l <= q   (for a total of q+1 possible values of l)
     The value of 'l' can be used as an index directly, but that for
     'k' needs to be shifted down by q.

   - In order to simplify other code (namely that in isoscalars.cc),
     we allow the code to run one space off the edge of each axis.
     This means we really allocate p+3 (or q+3) values along each axis
     (except for the degeneracy axis).

   - All reps in a degenerate set need to be processed at once, so we
     allocate storage all at once.

   - We can infer the value of l2 from the values of k,l,k1,l1,k1 using
     hypercharge conservation. As such, we don't need an axis for l2.
     We do check that the value provided is correct; TODO: make this optional.
*/

isoarray::isoarray(long p, long q, long p1, long q1, long p2, long q2)
    : p(p), q(q), p1(p1), q1(q1), p2(p2), q2(q2)
{
    d = degeneracy(p, q, p1, q1, p2, q2);

    size_t size = d * (p+3) * (q+3) * (p1+3) * (q1+3) * (p2+3);
    coefficients = new sqrat[size];

    /* Check that the array is properly initialised.
        TODO: Remove */
    assert(coefficients[111].numerator() == 0);
    assert(coefficients[111].denominator() == 1);
}

isoarray::~isoarray()
{
    delete coefficients;
}

/* We use operator() instead of operator[] as an easy way to use
    multiple indices */
sqrat& isoarray::operator()(long n, long k, long l, long k1, long l1,
                            long k2, long l2)
{
    /* Bounds checks */
    assert((n >= 0) && (n < d));
    assert((k >= q-1) && (k < p+q+1));
    assert((l >= -1) && (l < q+1));
    assert((k1 >= q1-1) && (k1 < p1+q1+1));
    assert((l1 >= -1) && (l1 < q1+1));
    assert((k2 >= q2-1) && (k2 < p2+q2+1));
    assert((l2 >= -1) && (l2 < q2+1));

    /* Check hypercharge conservation */
    assert(k1+l1+k2+l2 == (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3);

    /* Adjust indices */
    k = k+1-q;
    l = l+1;
    k1 = k1+1-q1;
    l1 = l1+1;
    k2 = k2+1-q2;
    l2 = l2+1;

    size_t index = ((((n * (p+3) + k) * (q+3) + l) * (p1+3) + k1)
                    * (q1+3) + l1) * (p2+1) + k2;
    return coefficients[index];
}
