/* libSU3: Miscellaneous helper functions */

#include "SU3_internal.h"

/* Min and max for various numbers of arguments */
long min(long a, long b)
{
    return (a < b) ? a : b;
}

long min(long a, long b, long c)
{
    return min(min(a, b), c);
}

long min(long a, long b, long c, long d, long e, long f, long g, long h, long i)
{
    return min(min(a, b, c), min(d, e, f), min(g, h, i));
}

long max(long a, long b)
{
    return (a > b) ? a : b;
}

long max(long a, long b, long c)
{
    return max(max(a, b), c);
}
