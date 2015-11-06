/* libSU3: Internal functions shared across multiple files */

#ifndef __SU3_INTERNAL_H__
#define __SU3_INTERNAL_H__

#include "SU3.h"

/* Various useful functions */
long min(long, long);
long min(long, long, long);
long min(long, long, long, long, long, long, long, long, long);
long max(long, long);
long max(long, long, long);

/* Macro to calculate (-1)^v */
#define SIGN(v) ((((v) % 2) == 0) ? 1 : -1)

/* A class for storing a bunch of useful values during our calculations.
    All functions are run as methods of an object of this class, so we have
    easy access to those values */
class isoscalar_context
{
public:
    /* Functions to get/set particular isoscalar factors. */
    sqrat isf(long n, long k, long l, long k1, long l1, long k2, long l2);
    void set_isf(long n, long k, long l, long k1, long l1, long k2, long l2,
                    sqrat value);

private:
    /* Values which many functions need access to */
    long p, q, p1, q1, p2, q2; // Target and factor reps
    long d; // Degeneracy
    long A; // = 1/3 (2(p1+p2) + 4(q1+q2) + (p-q))

    sqrat* coefficients;

    isoscalar_context(long p, long q, long p1, long q1, long p2, long q2,
                        long d, sqrat* coefficients);

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
    void c_coefficients(long k, long l, long k1, long l1, long k2, long l2,
                sqrat& alpha, sqrat& c1, sqrat& c2, sqrat& c3, sqrat& c4);
    void d_coefficients(long k, long k1, long l1, long k2, long l2,
                sqrat& beta, sqrat& d1, sqrat& d2, sqrat& d3);

    /* Use the A and B recursion relations to step along the
       k1 and l1 axes within a plane of constant s.

       The arguments identify the state to be calclated, *not* the values
       of k1,l1,k2,l2 used in the recursion relation itself.
    */
    void step_k1_up(long n, long s, long k1, long l1);
    void step_k1_down(long n, long s, long k1, long l1);
    void step_l1_up(long n, long s, long k1, long l1);
    void step_l1_down(long n, long s, long k1, long l1);

    /* Use the C and D recursion relations to step along the
       k and l axes within a multiplet.
       We fix l=0 in step_k_down, since this can *only* step down to such states.

       The arguments identify the state to be calclated, *not* the values
       of k1,l1,k2,l2 used in the recursion relation itself. */
    void step_k_down(long n, long k, long k1, long l1, long k2, long l2);
    void step_l_up(long n, long k, long l, long k1, long l1, long k2, long l2);

    /* Step down from one plane (at s+2) to the next plane (at s).
       In principle this can fail, in which case you need to use the exchange
       symmetries in order to calculate the values. This should never happen with
       this function, however, as isoscalars() has logic for deciding whether this
       can happen *before* it picks which ISFs to calculate.

       Note: When doing this, there are two independent choices we can make:

       * We can try to calculate the value at (k1max,l1min) or at (k1min,l1max)
         (these require different steps afterwards)

       * We can try recursion relation A or recursion relation B
         (these require the same steps afterwards)

       Each of the four combinations works in different cases, so we intelligently
       choose between them.

       Internally, we use the step_k1_up/down and step_l1_up/down functions,
       but we step "from" a non-existent state (with k1,l1 not in the valid range)
       to the state we want. This works because 'isoscalar_array' allows us to
       request (certain) non-existent states and just returns 0 for the coupling
       coefficient. This is exactly what we need for the stepping to work properly.
    */
    void step_s_down(long n, long s);

    /* Calculate the inner product of two sets of isoscalar factors */
    sqrat inner_product(long m, long n);

    /* Calculate couplings to the state of highest weight.
        This can throw std::logic_error if we can't calculate directly. This should
        never happen, however, as isoscalars() has logic to avoid those cases.
    */
    void calc_shw();

    /* Fill out each multiplet, assuming that the SHWs have been calculated */
    void calc_isoscalars();

    /* Allow the top-level driver function to interact with
        objects of this class */
    friend isoarray* isoscalars_single(long p, long q, long p1, long q1,
                                        long p2, long q2, long d);
};

#endif
