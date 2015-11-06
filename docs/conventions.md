Conventions
===========

Here we document the conventions used in this library. We attempt to follow
Kaeding and Williams, but it is worth reiterating everything here.

Labelling eigenstates
---------------------

We label irreducible representations of SU(3) by (p, q) throughout.

Rather than the usual variables of hypercharge y, isospin i and isospin
z-component iz, we instead use the following variables:

    k = (p+2q)/3 + y/2 + i
    l = (p+2q)/3 + y/2 - i
    m = (p+2q)/3 + y/2 + iz

These can be inverted to give:

    y = k + l - 2(p+2q)/3
    i = (k - l)/2
    iz = (2m - k - l)/2

The advantage is that k,l,m are always integers, with ranges:

    0 <= l <= q
    q <= k <= p+q
    l <= m <= k

We pick the state of highest weight in each representation to be the state
|SHW> with the largest iz. Equivalently, it is the state such that
V+ |SHW> = 0, U- |SHW> = 0.

This corresponds to the state k = p+q, l = 0.

Phase conventions
-----------------

We use the de Swart phase conventions. This means that:

* All Clebsch-Gordan coefficients and isoscalar factors are real

* For SU(2), we use the Condon-Shortley phase convention
  <J J | j1 j1 j2 (J-j1)> > 0

* The relative phases of states within an SU(3) representation are fixed
  by applying the Condon-Shortley phase convention to both isospin and
  V-spin (ie, between u and s quarks)

* For SU(3), the phase convention for the isoscalar factors
  (dropping p,q labels and only labelling k,l for each representation) is
  F(p+q, 0, p1+q1, 0, k2max, l2min) > 0,
  where k2max, l2min are the largest k2 and smallest l2 (respectively) such that
  the isoscalar factor is nonzero. There is always at least one such set of values.

Degeneracy resolution
---------------------

Sometimes, a tensor product of representations (p1,q1) with (p2,q2)
can include multiple copies of the same representation (p,q). We call
the latter the "target representation".

In this case, we use Kaeding and Williams's method to resolve the degeneracy
(best explained in Williams's paper). This means that the inner product of any
two states from different copies of the target representations are orthogonal.

The only remaining thing to specify is the ordering of degenerate representations.
We order them the same as Kaeding and Williams, which is exactly the reverse order
to the one used by, eg, de Swart.
