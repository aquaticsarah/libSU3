/* libSU3: State-of-higest-weight (SHW) determination.
   Note: This also does degeneracy resolution */

#include <stdio.h>

#include "SU3_internal.h"

isoscalar_context::isoscalar_context(isoarray* isf, long p, long q, long p1,
            long q1, long p2, long q2) : p(p), q(q), p1(p1),
            q1(q1), p2(p2), q2(q2), isf(isf)
{
    d = degeneracy(p, q, p1, q1, p2, q2);
    A = (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;
}

/* Main calculation function */
isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (! d) return NULL; /* Ignore reps of zero degeneracy */

    isoarray* isf = new isoarray(p,q,p1,q1,p2,q2);
    isoscalar_context* ctx = new isoscalar_context(isf,p,q,p1,q1,p2,q2);

    ctx->calc_shw();
    /* TODO: Deal with the case where we need to use the conjugated reps */

    /* TODO: Calculate the rest of the ISFs, given those for the SHW */

    delete ctx;
    return isf;
}
