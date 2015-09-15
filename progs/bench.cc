/* libSU3: Benchmarking program */

#include <time.h>

#include "SU3.h"

#define ITERS 50L
#define DELTA(start, end) ((end - start) / (double)CLOCKS_PER_SEC)

int main()
{
    clock_t start, end;
    double overhead, raw, actual;

    long i;
    long d;
    long n,k,l,m,k1,l1,m1,k2,l2,m2;

    /* Time an empty loop, to estimate overhead */
    printf("Timing overhead (%ld iterations)...\n", ITERS);

    start = clock();
    for (i = 0; i < ITERS; ++i)
        asm volatile ("");
    end = clock();
    overhead = DELTA(start, end);

    /* Benchmark the calculation of 27 x 27 -> 27 */
    /* d = degeneracy * (p+1) * (q+1) * (p1+1) * (q1+1) * (p2+1)
        (the q2 axis is ignored - due to hypercharge conservation, only one
        value of l2 is valid for any given k,l,k1,l1,k2)
    */
    d = 3*3*3*3*3*3;
    printf("Timing ISFs for 27x27->27 (%ld iterations, %ld ISFs)...\n", ITERS, d);

    isoarray* isf;
    start = clock();
    for (i = 0; i < ITERS; ++i)
    {
        isf = isoscalars(2, 2, 2, 2, 2, 2);
        delete isf;
    }
    end = clock();
    raw = DELTA(start, end);
    actual = raw - overhead;

    printf("Total time: %7.3fs = %7.3fms/iter\n", actual, actual*1000./ITERS);
    printf("                     = %7.3fus/ISF\n\n", actual*1000000./(ITERS*d));

    /* Benchmark conversion from ISFs to CGCs */
    cgarray* cgc = clebsch_gordans(2, 2, 2, 2, 2, 2);
    d = 3*(3*3*6/2)*(3*3*6/2)*(3*3*6/2);
    printf("Timing ISF->CGC for 27x27->27 (%ld iterations, %ld CGCs)...\n", ITERS, d);

    start = clock();
    for (i = 0; i < ITERS; ++i)
    {
        for (n = 0; n < 2; ++n)
            for (k = 2; k <= 4; ++k)
                for (l = 0; l <= 2; ++l)
                    for (m = l; m <= k; ++m)
                        for (k1 = 2; k1 <= 4; ++k1)
                            for (l1 = 0; l1 <= 2; ++l1)
                                for (m1 = l1; m1 <= k1; ++m1)
                                    for (k2 = 2; k2 <= 4; ++k2)
                                        for (l2 = 0; l2 <= 2; ++l2)
                                            for (m2 = l2; m2 <= k2; ++m2)
                                            {
                                                (*cgc)(n,k,l,m,k1,l1,m1,k2,l2,m2);
                                            }
    }
    end = clock();
    raw = DELTA(start, end);
    actual = raw - overhead;
    delete cgc;

    printf("Total time: %7.3fs = %7.3fms/iter\n", actual, actual*1000./ITERS);
    printf("                     = %7.3fus/CGC\n\n", actual*1000000./(ITERS*d));
}
