/* libSU3: Tests for information about irreps of SU(3) */

#include <stdio.h>

#include "SU3_internal.h"

/* Calculate the dimension of one irrep */
long dimension(long p, long q)
{
    return (p+1)*(q+1)*(p+q+2)/2;
}

/* Display the name of a representation.
    TODO: Deal with primed reps */
char* repname(char* buffer, size_t len, long p, long q)
{
    long dim = dimension(p, q);
    if (q > p)
        snprintf(buffer, len, "%ld*", dim);
    else
        snprintf(buffer, len, "%ld", dim);
    return buffer;
}

/* Calculate the degeneracy of the (p,q) irrep in the decomposition of
   (p1,q1) x (p2,q2). Returns 0 if (p,q) is not a summand in this
   decomposition.
   Notation and conditions are based on arXiv:nucl-th/9511025 */
long degeneracy(long p, long q, long p1, long q1, long p2, long q2)
{
    long x = p1 + p2 - p;
    long y = q1 + q2 - q;

    /* Reps can only appear if x-y is a multiple of three
       (and therefore so is x+2y) */
    if ((x-y) % 3) return 0;

    long a = (x - y) / 3;
    long b = (x + 2*y) / 3;

    /* More conditions which exclude reps from appearing */
    if ((b < 0) || (b > min(q1+q2, p1+q1, p2+q2))) return 0;
    if ((a < -min(q1, q2)) || (a > max(p1, p2))) return 0;
    if ((a+b < 0) || (a+b > min(p1+p2, p1+q1, p2+q2))) return 0;

    /* If we get here, we have a formula for the degeneracy of the
        rep; the rep appears iff this value is positive */
    long d = 1 + min(q2, p1+q1, b, p1-a)
               - max(0, b-q1, b-p2, -a, b-a-q1, a+b-p2);

    return max(0, d);
}
