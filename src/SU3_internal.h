/* libSU3: Internal functions shared across multiple files */

#ifndef __SU3_INTERNAL_H__
#define __SU3_INTERNAL_H__

#include "SU3.h"

/* Various useful functions */
long min(long a, long b);
long min(long a, long b, long c);
long min(long a, long b, long c, long d);
long max(long a, long b);
long max(long a, long b, long c, long d, long e, long f);
long gcd(long, long);

/* A class for storing a bunch of useful values during our calculations.
    All functions are run as methods of an object of this class, so we have
    easy access to those values */
class isoscalar_context
{
private:
    /* Values which many functions need access to */
    long p, q, p1, q1, p2, q2; // Target and factor reps
    long d; // Degeneracy
    long A; // = 1/3 (2(p1+p2) + 4(q1+q2) + (p-q))

    isoarray* isf;

    /* Calculate the coefficients for each of the four recursion relations.
       Each stores the coefficients in its last four arguments.
       The expressions are based on arXiv:nucl-th/9511025.
        
       Note: The functions step_s_down_* often make steps "from" locations
       which do not correspond to valid states. In those cases, we know that
       the corresponding coefficient will always be ignored. However, in those
       cases, the coefficient is sometimes indeterminate or infinite, and this
       causes problems. Specifically, this occurs if one of k1-l1 or k2-l2 is
       zero.

       Because the coefficient is ignored in those cases, we can safely
       special-case it to be zero instead.
    */
    void a_coefficients(long k1, long l1, long k2, long l2,
                        sqrat& a1, sqrat& a2, sqrat& a3, sqrat& a4);
    void b_coefficients(long k1, long l1, long k2, long l2,
                        sqrat& b1, sqrat& b2, sqrat& b3, sqrat& b4);

    /* Use the A and B recursion relations to step along the
       k1 and l1 axes within a plane.

       The arguments identify the state to be calclated, *not* the values
       of k1,l1,k2,l2 used in the recursion relation itself.
    */
    void step_k1_up(long n, long s, long k1, long l1);
    void step_k1_down(long n, long s, long k1, long l1);
    void step_l1_up(long n, long s, long k1, long l1);
    void step_l1_down(long n, long s, long k1, long l1);

    /* Step down from one plane (at s+2) to the next plane (at s).
       Sometimes this can fail, in which case you need to conjugate all reps,
       try again (this will always succeed in this case) and then use symmetries
       to extract the original coupling coefficients.
       The function returns 1 on sucess, 0 on failure.
    */
    int step_s_down(long n, long s,
            long k1min, long k1max, long l1min, long l1max);

    /* Calculate the inner product of two sets of isoscalar factors */
    sqrat inner_product(long m, long n);

protected:
    isoscalar_context(isoarray* isf, long p, long q, long p1,
                    long q1, long p2, long q2);
    int calc_shw();

    /* Allow the top-level driver function to interact with
        objects of this class */
    friend isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2);
};

#endif
