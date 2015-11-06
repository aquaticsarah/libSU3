In Williams, the quantity s = k1 - l1 + k2 - l2 is defined, and is used along
with hypercharge conservation (k1 + l1 + k2 + l2 = constant, which we will
call A) to specify a state in the target representation using the variables
(s, k1, l1). There is then a comment:

"In this description, smax is the maximum value of s for which a non-vanishing
ISF occurs for a coupling of fixed p,q,p1,q1,p2,q2 values:
k1min, l1min, k1max, l1max are the minimum and maximum values of the k1 and l1
variables for a particular value of s."

However, derivations of these values are not given. Here we detail how to derive
them.

Deriving min/max values
=======================

We can rearrange the defining equations for A, s to get:

    k1 = (A+s)/2 - k2
    l1 = (A-s)/2 - l2

Then we use the valid ranges for k2, l2: We have saturated bounds

    q2 <= k2 <= p2+q2
    0 <= l2 <= q2

which imply:

    (A+s)/2 - (p2+q2) <= k1 <= (A+s)/2 - q2
    (A-s)/2 - q2 <= l1 <= (A-s)/2

These upper/lower bounds are what we will call the "unrestricted" bounds on
k1, l1, and we will label them with "_u" at the end. That is, we write

    k1min_u = (A+s)/2 - (p2+q2)
    k1max_u = (A+s)/2 - q2
    l1min_u = (A-s)/2 - q2
    l1max_u = (A-s)/2

(+): Note in particular that, when we step down between planes, we decrease s
by two. Thus we decrease k1{min/max}* by 1 and increase l1{min/max}* by 1.
This will be important later.

Since we only have valid states in our first representation when
q1 <= k1 <= p1+q1, 0 <= l1 <= q1, we then derive the true values

    k1min = max(q1, k1min_u)
    k1max = min(p1+q1, k1max_u)
    l1min = max(0, l1min_u)
    l1max = min(q2, l1max_u)

Deriving s_min, s_max
=====================

Further to the above, we clearly need to have at least one valid state in
the first representation to couple to. This translates to a set of conditions:

    k1min_u <= p1+q1, which implies (A+s)/2 <= p1+q1+p2+q2
    k1max_u >= q1, which implies (A+s)/2 >= q1+q2
    l1min_u <= q1, which implies (A-s)/2 <= q1+q2
    l1max_u >= 0, which implies (A-s)/2 >= 0

We also have isospin conservation, which implies:

abs(k1-l1-k2+l2) <= k-l <= k1-l1+k2-l2

Note that the RHS is s, and we are working with the state of highest
weight in the target representation - that is, k=p+q and l=0.
Thus we also need s >= p+q.

Rearranging all of the above for s and combining them in an expression
of the form smin <= s <= smax gives:

    smax = min(2*(p1+q1+p2+q2) - A, A) = min(Ar, A)
    smin = max(p+q, 2*(q1+q2) - A, A - 2*(q1+q2)) = max(p+q, abs(A - 2*(q1+q2)))

where Ar is the same expression as A but with 'p's and 'q's swapped.

When does the direct method work?
=================================

Here we wish to answer: When can we calculate the ISFs for a particular
combination of irreps directly, rather than having to use symmetry relations?

To answer this question, note the following:
* The calculation can only fail when we try to "step down" through the planes
  of constant s. Within a plane, once we have one value we can always fill out
  the rest of that plane.
* If the degeneracy is d, our degeneracy resolution scheme fills in one value
  on each of the highest d planes automatically. So these planes can never fail.
* We have four different ways we can try to step down (see step_s_down in shw.cc),
  each of which has different success criteria.

So we need to find these success criteria. By looking carefully at the values
used by each recursion relation, using the relation (+) between the limits
of k1,l1 on adjacent planes, and considering the fixed bounds on k1,l1,
we find that:

* Using step_k1_up to find the value at (k1min, l1max) works iff the plane
  above has a value at (k1min, l1max) or at (k1min, l1max-1).
  This happens iff k1min_u < q1

* Using step_l1_down to find the value at (k1min, l1max) works iff the plane
  above has a value at (k1min, l1max) or at (k1min+1, l1max).
  This happens iff l1max_u > q1

* Using step_k1_down to find the value at (k1max, l1min) works iff the plane
  above has a value at (k1max+1, l1min) or at (k1max+1, l1min-1).
  This happens iff k1max_u < p1+q1

* Using step_l1_up to find the value at (k1max, l1min) works iff the plane
  above has a value at (k1max, l1min-1) or at (k1max+1, l1min-1).
  This happens iff l1min_u > 0

Further, we note that since k1{min/max} decrease by 1 and l1{min/max} increase
by 1 each time we step down by a plane, if we successfully fill one plane
(which is not filled by our degeneracy resolution) then we are guaranteed to
succeed on all lower planes.

Hence the stepping down process will succeed iff any of the following five
conditions hold:

    smax - 2*d == smin, ie. the degeneracy resolution fills *every* plane
    k1min_u < q1 at s = smax - 2*d
    l1max_u > q1 at s = smax - 2*d
    k1max_u < p1+q1 at s = smax - 2*d
    l1min_u > 0 at s = smax - 2*d

We use these criteria in isoscalars.cc to determine what calculational
method to use before actually performing any calculations. This simplifies
the SHW calculation code, since it is then guaranteed to succeed.
