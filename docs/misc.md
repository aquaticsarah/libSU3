Miscellaneous notes
===================

Meaning of 's'
--------------

When calculating couplings to the state of highest weight, we define a
variable s = k1-l1+k2-l2. This is used together with k1,l1 to label tensor
product states. Here we try to explain how to interpret this:

Picture a 3d space with horizontal axes (k1,l1) and vertical axis s.

For any particular value of s, we can use hypercharge conservation and the
bounds on k2, l2 (coming from requiring them to correspond to states which
actually exist), we can derive some bounds on k1, l1:

* k1min <= k1 <= k1max
* l1min <= l1 <= l1max

(for the values of these bounds)

We interpret this as a (finite) horizontal plane in our 3d space. There
is a tower of such planes, with some parts of each plane lying vertically
above parts of the planes below.

The lattice points which lie on one of the planes correspond to values of
(k1,l1) which can couple (to some states (k2,l2)) to produce a contribution
to the state of highest weight.

Note that each of the A and B recurrence relations (as named in coeff.cc)
relate two values in one plane (at a particular value of s) to two values
in the plane immediately below (at s-2). Sometimes the values which would
be related do not lie on a plane. In that case the relevant coefficient
is zero - essentially, that state cannot couple.

The recursion relations can be used as long as they relate at least two
different values. In the file derivations.md, we derive conditions under
which they allow us to successfully step down through all of the planes
and thereby calculate all couplings to the state of highest weight.

Indirect calculation
--------------------

The direct method of calculation sometimes fails. The ArXiV version
of Kaeding and Williams states that "In every such case, the algorithm
succeeds for the conjugated ISFs, [...]".

This line is slightly incorrect - using the conjugated ISFs does not
necessarily help us. Instead, we use the exchange relations, ie. the
ones relating the ISFs for one combination of irreps to those for
a different combination of the same irreps and their conjugates.

Specifically, we use the relation which Williams calls (1 <-> 3bar),
and a similar one which we call (2 <-> 3bar) which we derive
by performing the sequence (1 <-> 2), (1 <-> 3bar), (1 <-> 2).

We have not proven that these two alternatives (as well as the direct
calculation method) are sufficient in all cases, but the library is
programmed to throw an exception if it turns out to be insufficient.
