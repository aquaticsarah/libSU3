/* libSU3: Tests for information about irreps of SU(3) */

#include <stdio.h>
#include <assert.h>

#include "SU3_internal.h"

/* Calculate the dimension of one irrep */
long dimension(long p, long q)
{
    return (p+1)*(q+1)*(p+q+2)/2;
}

/* Calculate the degeneracy of the (p,q) irrep in the decomposition of
    (p1,q1) x (p2,q2). Returns 0 if (p,q) is not a summand in this
    decomposition.
    Note: The intermediate values used in this function are based on those in
    arXiv:hep-th/9509167, but are all multiplied by 3 relative to that paper.
    We also define delta = gamma + sigma.
*/
long degeneracy(long p, long q, long p1, long q1, long p2, long q2)
{
    long gamma = (p1 + p2 - p);
    long sigma = (q1 + q2 - q);
    long delta = gamma + sigma;

    /* Reps can only appear if gamma-sigma is a multiple of three */
    if ((gamma - sigma) % 3) return 0;

    /* Aside: This is invariant under all of the symmetry transformations.
        See docs/degeneracy.md for proof.
    */
    long eta_prime = min(3*p1 + sigma, 3*p2 + sigma, 3*q + sigma,
                         3*q1 + gamma, 3*q2 + gamma, 3*p + gamma,
                         2*delta, 3*(p1+q1) - delta, 3*(p2+q2) - delta);

    /* If we get here, we have a formula for the degeneracy of the
        rep; the rep appears iff this value is positive */
    long eta = max(eta_prime + 3 - max(gamma, sigma), 0);

    /* Eta should now be 3*(degeneracy), and the degeneracy is an integer */
    assert(! (eta % 3));
    return eta/3;
}

/* Phase changes under the 1<->2 symmetry and the conjugation symmetry */
long phase_exch_12(long p, long q, long p1, long q1, long p2, long q2)
{
    long gamma = (p1 + p2 - p);
    long sigma = (q1 + q2 - q);

    /* Note: The phase as given in Williams, converted to our variables, is
        (-1)^(gamma/3 + sigma/3 + max(gamma/3, sigma/3)),
        where the exponent is guaranteed to be an integer. Hence we can multiply
        it by 3 to get an equivalent phase of:
        (-1)^(gamma + sigma + max(gamma, sigma))
     == (-1)^(2*max(gamma, sigma) + min(gamma, sigma))
     == (-1)^min(gamma, sigma)
    */
    return SIGN(min(gamma, sigma));
}

long phase_conj(long p, long q, long p1, long q1, long p2, long q2)
{
    long gamma = (p1 + p2 - p);
    long sigma = (q1 + q2 - q);

    /* Note: Williams' expression for this phase is the same as above,
        except with max replaced by min. A similar argument simplifies
        it to (-1)^max(gamma, sigma)
    */
    return SIGN(max(gamma, sigma));
}
