/* libSU3: Benchmarking program */

#include <time.h>

#include "SU3.h"

#define ITERS 25L
#define DELTA(start, end) ((end - start) / (double)CLOCKS_PER_SEC)

int main()
{
    clock_t start, end;
    double elapsed;
    long i;

    printf("Running benchmarks; all results are averages over %ld iterations.\n\n", ITERS);
    printf("Timing calculations for small reps...\n");

    isoarray* isf;
    long p, q, p1, q1, p2, q2;
    start = clock();
    for (i = 0; i < ITERS; ++i)
    {
        for (p = 0; p < 3; ++p)
            for (q = 0; q < 3; ++q)
                for (p1 = 0; p1 < 3; ++p1)
                    for (q1 = 0; q1 < 3; ++q1)
                        for (p2 = 0; p2 < 3; ++p2)
                            for (q2 = 0; q2 < 3; ++q2)
                            {
                                isf = isoscalars(p, q, p1, q1, p2, q2);
                                delete isf;
                            }
    }
    end = clock();
    elapsed = DELTA(start, end);

    printf("Total time: %7.3fs = %7.3fms/iter\n\n", elapsed, elapsed*1000./ITERS);

    printf("Timing ISF->CGC conversion for 27x27->27...\n");

    cgarray* cgc = clebsch_gordans(2, 2, 2, 2, 2, 2);
    long n,k,l,m,k1,l1,m1,k2,l2,m2;
    start = clock();
    for (i = 0; i < ITERS; ++i)
    {
        for (n = 0; n < 2; ++n)
            FOREACH_CGC(2, 2, 2, 2, 2, 2, k, l, m, k1, l1, m1, k2, l2, m2)
                (*cgc)(n,k,l,m,k1,l1,m1,k2,l2,m2);
    }
    end = clock();
    elapsed = DELTA(start, end);
    delete cgc;

    printf("Total time: %7.3fs = %7.3fms/iter\n", elapsed, elapsed*1000./ITERS);
}
