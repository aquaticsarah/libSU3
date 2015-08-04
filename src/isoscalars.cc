/* libSU3: State-of-higest-weight (SHW) determination.
   Note: This also does degeneracy resolution */

#include <stdio.h>
#include <stdlib.h>

#include "SU3.h"

static void calc_shw(isoarray* isf, long d, long p, long q,
                    long p1, long q1, long p2, long q2)
{
    long A =           (2*p1 + 2*p2 + 4*q1 + 4*q2 + p - q)/3;
    long smax = min(A, (2*q1 + 2*q2 + 4*p1 + 4*p2 + q - p)/3);
    long smin = max(p + q, 2*q1 + 2*q2 - A);

    long k1min, k1max, l1min, l1max;

/* Macro to simplify indexing */
#define SHW(n, s, k1, l1) (*isf)(n, p+q, 0, k1, l1, (A+s)/2 - k1, (A-s)/2 - l1)

    /* Fill out the topmost d planes for each of the degenerate reps */
    long m, n, s;
    for (m = 0; m < d; ++m)
    {
        s = smax - 2*m;
        k1min = max(q1, (A + s)/2 - (p2+q2));
        k1max = min(p1+q1, (A + s)/2 - q2);
        l1min = max(0, (A - s)/2 - q2);
        l1max = min(q1, (A - s)/2);

        /* Set one ISF in one particular irrep (leaving the same ISF
           in the other irreps as zero) */
        SHW(m, s, k1min, l1min) = 1;

        for (n = 0; n < d; ++n)
        {
            /* Use recursion relations (possibly involving the plane
               above the current one, which will already have been filled)
               to fill out the rest of this plane */
            (void)k1max; (void)l1max;
            /* TODO */
        }
    }

    /* Now we have filled out the topmost d planes, step down
       through the rest of them */
    for (s = smax - 2*d; s >= smin; s -= 2)
    {
        /* TODO */
    }

    /* Orthonormalise */
    /* TODO */
}

void isoscalars(long p, long q, long p1, long q1, long p2, long q2)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (! d) return; /* Ignore reps of zero degeneracy */

    char name1[16], name2[16], name3[16];
    repname(name1, 16, p1, q1);
    repname(name2, 16, p2, q2);
    repname(name3, 16, p, q);

    printf("Calculating ISFs for %s x %s -> %s, degeneracy %ld\n",
            name1, name2, name3, d);

    isoarray* isf = new isoarray(p,q,p1,q1,p2,q2);

    calc_shw(isf, d, p, q, p1, q1, p2, q2);

    /* TODO: Calculate the rest of the ISFs, given those for the SHW */

    delete isf;
}
