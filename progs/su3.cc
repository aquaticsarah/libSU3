/* libSU3: Demo program to print ISFs and CGCs */

#include <stdio.h>
#include <stdlib.h>

#include "SU3.h"

#define BUF_SIZE 128

enum print_mode
{
    MODE_ISF,
    MODE_CGC
};

const char* usage_message = "\
Usage: %s [isf|cgc] (p1,q1)x(p2,q2)[->(p,q)[,n]]\n\
\n\
If the last representation (the part after '->') is not specified, we do the\n\
full decomposition of (p1,q1) x (p2,q2). Otherwise, we only print the values\n\
for the given target rep.\n\
\n\
If n is given, then we only print the nth degenerate representation (n must be\n\
in the range 1,...,d where d is the degeneracy). Otherwise (or if n is given\n\
as -1) we print all degenerate representations.\n\
";

char* progname;

void usage()
{
    printf(usage_message, progname);
    exit(2);
}

/* Display values for a single irrep in a possibly-degenerate set */
void print_isfs(isoarray* isf, long n, long p, long q,
                long p1, long q1, long p2, long q2, long d)
{
    long k, l, k1, l1, k2, l2;
    char val_buf[BUF_SIZE];

    if (d > 1)
        printf("  Degenerate rep %ld/%ld:\n", n, d);

    FOREACH_ISF(p, q, p1, q1, p2, q2, k, l, k1, l1, k2, l2)
    {
        (*isf)(n-1, k, l, k1, l1, k2, l2).tostring(val_buf, BUF_SIZE);
        printf("    (%ld,%ld) : (%ld,%ld) x (%ld,%ld) = %s\n",
            k, l, k1, l1, k2, l2, val_buf);
    }

    printf("\n");
}

void print_cgcs(cgarray* cg, long n, long p, long q,
                long p1, long q1, long p2, long q2, long d)
{
    long k, l, m, k1, l1, m1, k2, l2, m2;
    char val_buf[BUF_SIZE];

    if (d > 1)
        printf("  Degenerate rep %ld/%ld:\n", n, d);

    FOREACH_CGC(p, q, p1, q1, p2, q2, k, l, m, k1, l1, m1, k2, l2, m2)
    {
        (*cg)(n-1, k, l, m, k1, l1, m1, k2, l2, m2).tostring(val_buf, BUF_SIZE);
        printf("    (%ld,%ld,%ld) : (%ld,%ld,%ld) x (%ld,%ld,%ld) = %s\n",
            k, l, m, k1, l1, m1, k2, l2, m2, val_buf);
    }

    printf("\n");
}

/* Print the values for a particular target representation, but possibly with degeneracy >1. */
void do_rep(long p, long q, long p1, long q1, long p2, long q2, long n, enum print_mode mode)
{
    long d = degeneracy(p, q, p1, q1, p2, q2);
    if (d == 0) return;
    if ((n < -1) || (n == 0) || (n > d))
    {
        printf("Error: Degeneracy label %ld is not in valid range 1,...,%ld\n",
                n, d);
        return;
    }

    if (mode == MODE_ISF)
    {
        printf("Isoscalar factors for (%ld,%ld) x (%ld,%ld) -> (%ld,%ld), degeneracy %ld:\n",
                p1, q1, p2, q2, p, q, d);

        isoarray* isf = isoscalars(p, q, p1, q1, p2, q2);
        if (n > 0)
            print_isfs(isf, n, p, q, p1, q1, p2, q2, d);
        else
        {
            /* Print all reps */
            for (n = 1; n <= d; ++n)
                print_isfs(isf, n, p, q, p1, q1, p2, q2, d);
        }

        delete isf;
    }
    else if (mode == MODE_CGC)
    {
        printf("Clebsch-Gordan coefficients for (%ld,%ld) x (%ld,%ld) -> (%ld,%ld), degeneracy %ld:\n",
                p1, q1, p2, q2, p, q, d);

        cgarray* cg = clebsch_gordans(p, q, p1, q1, p2, q2);
        if (n > 0)
            print_cgcs(cg, n, p, q, p1, q1, p2, q2, d);
        else
        {
            /* Print all reps */
            for (n = 1; n <= d; ++n)
                print_cgcs(cg, n, p, q, p1, q1, p2, q2, d);
        }

        delete cg;
    }
}

/* Print the values requested by the given command and mode settings. */
void print_values(char* command, enum print_mode mode)
{
    long p1, q1, p2, q2, p=-1, q=-1, n=-1;

    int res = sscanf(command, "(%ld,%ld)x(%ld,%ld)->(%ld,%ld),%ld",
                        &p1, &q1, &p2, &q2, &p, &q, &n);

    /* We only allow certain numbers of parameters */
    if ((res < 4) || (res == 5) || (res > 7))
        usage();

    if (p == -1)
    {
        /* Crude bounds on which reps can appear */
        long upper = p1+q1+p2+q2;

        /* Print a decomposition at the top */
        int any_printed = 0; // Keeps track of whether we need to print a '+'
        printf("Clebsch-Gordan series:\n");
        printf("(%ld,%ld) x (%ld,%ld) =", p1, q1, p2, q2);
        for (p = 0; p <= upper; ++p)
            for (q = 0; q <= upper; ++q)
            {
                long d = degeneracy(p, q, p1, q1, p2, q2);

                if (d == 0) continue;
                else if (d == 1)
                {
                    printf("%s (%ld,%ld)", any_printed ? " +" : "", p, q);
                    any_printed = 1;
                }
                else
                {
                    printf("%s %ldx(%ld,%ld)",
                            any_printed ? " +" : "", d, p, q);
                    any_printed = 1;
                }
            }
        printf("\n\n");

        /* Print individual values */
        for (p = 0; p <= upper; ++p)
            for (q = 0; q <= upper; ++q)
                do_rep(p, q, p1, q1, p2, q2, -1, mode);
    }
    else
        do_rep(p, q, p1, q1, p2, q2, n, mode);
}

/* Main function: Just check that the arguments are in the right format and
    pass them along
*/
int main(int argc, char** argv)
{
    progname = argv[0];

    if ((argc < 2) || (argc > 3))
        usage();

    if (argc == 2)
    {
        print_values(argv[1], MODE_ISF);
    }
    else
    {
        if (! strcmp(argv[1], "isf"))
            print_values(argv[2], MODE_ISF);
        else if (! strcmp(argv[1], "cgc"))
            print_values(argv[2], MODE_CGC);
        else
            usage();
    }

    return 0;
}
