Here we prove that the degeneracy is invariant under symmetry transformations.

From Williams, we have the expressions:

    gamma = (p1 + p2 - p)/3, sigma = (q1 + q2 - q)/3

    eta = max(eta' + 1 - max(gamma, sigma), 0)

    eta' = min(p1 + sigma, p2 + sigma, q + sigma,
               q1 + gamma, q2 + gamma, p + gamma,
               2(sigma+gamma), p1+q1-gamma-sigma, p2+q2-gamma-sigma)

We now investigate what happens under each of the symmetry transforms

1 <-> 2
-------

Under this relation, everything (gamma, sigma, eta' and thus eta) is invariant.

1 <-> 3bar
----------

    gamma_new = (q + p2 - q1)/3 = gamma + (p+q-p1-q1)/3

    sigma_new = (p + q2 - p1)/3 = sigma + (p+q-p1-q1)/3

    eta'_new = min(q + sigma_new, p2 + sigma_new, p1 + sigma_new,
                   p + gamma_new, q2 + gamma_new, q1 + gamma_new,
                   2(sigma_new+gamma_new), q+p-gamma_new-sigma_new, p2+q2-gamma_new-sigma_new)

             = min(q + sigma, p2 + sigma, p1 + sigma,
                   p + gamma, q2 + gamma, q1 + gamma,
                   2(sigma+gamma)+(p+q-p1-q2)/3, q+p-gamma-sigma-(p+q-p1-q1), p2+q2-gamma-sigma-(p+q-p1-q1)) + (p+q-p1-q1)/3

             = min(q + sigma, p2 + sigma, p1 + sigma,
                   p + gamma, q2 + gamma, q1 + gamma,
                   p2+q2-gamma-sigma, p1+q1-gamma-sigma, 2(gamma+sigma)) + (p+q-p1-q1)/3

             = eta' + (p+q-p1-q1)/3

    eta_new = max(eta'_new + 1 - max(gamma_new, sigma_new), 0)
            = max(eta' + 1 - max(gamma, sigma), 0) by cancelling the extra terms

2 <-> 3bar
----------

This is a combination of the previous two, so the degeneracy is invariant under this

Conjugation
-----------

This exchanges gamma and sigma, and leaves eta', eta invariant.
Thus the degeneracy is invariant under this.
