/* libSU3: Calculation of SU(2) Clebsch-Gordan coefficients.
    This is based on a formula from arXiv:nucl-th/9511025.
*/

#include <stdexcept>

#include "SU3_internal.h"

/* Extract an integer from an MPQ. Throws an exception if the value
    isn't actually an integer.
*/
static long mpq_to_long(mpq_class x)
{
    mpz_class num = x.get_num(), denom = x.get_den();
    if (denom != 1)
        throw std::domain_error("Isospin value is not a half-integer.");

    return num.get_si();
}

static mpz_class factorial(long x)
{
    if (x < 0) throw std::domain_error("Factorial of a negative number.");

    mpz_class res;
    mpz_fac_ui(res.get_mpz_t(), x);
    return res;
}

/* Calculate a single SU(2) Clebsch-Gordan coefficient.
    All arguments are implicitly doubled - eg, I represents
    2*(the actual isospin).

    The expression used is from Bohm - see README for citation.
*/
sqrat su2_cgc_2i(long I, long Iz, long i1, long i1z,
                    long i2, long i2z)
{
    /* Only states obeying the following conditions can couple:
        Iz = i1z + i2z
        i1 + i2 >= I >= |i1 - i2|
        Note |i1 - i2| = max(i1-i2, i2-i1), and so we can rewrite the last
        inequality as (I >= i1 - i2) && (I >= i2 - i1)
    */
    if (Iz != i1z + i2z) return sqrat(0);
    if ((I > i1 + i2) || (I < i1 - i2) || (I < i2 - i1)) return sqrat(0);

    mpz_class numerator = mpz_class(I + 1) * factorial((I + i1 - i2)/2)
                        * factorial((I - i1 + i2)/2) * factorial((i1 + i2 - I)/2)
                        * factorial((i1 + i1z)/2) * factorial((i1 - i1z)/2)
                        * factorial((i2 + i2z)/2) * factorial((i2 - i2z)/2)
                        * factorial((I + Iz)/2) * factorial((I - Iz)/2);
    mpz_class denominator = factorial((I + i1 + i2)/2 + 1);
    sqrat prefactor = sqrat(numerator, denominator);

    mpq_class sum = 0;
    long zmin = max(0, i2 - i1z - I, i1 + i2z - I)/2;
    long zmax = min(i1 + i2 - I, i1 - i1z, i2 + i2z)/2;
    long z;

    for (z = zmin; z <= zmax; ++z)
    {
        denominator = factorial(z) * factorial((i1 + i2 - I)/2 - z)
                    * factorial((i1 - i1z)/2 - z) * factorial((i2 + i2z)/2 - z)
                    * factorial((I - i2 + i1z)/2 + z) * factorial((I - i1 - i2z)/2 + z);

        sum += mpq_class(SIGN(z), denominator);
    }

    return prefactor * sum;
}

/* Calculate a single SU(2) Clebsch-Gordan coefficient.
    This does *not* take doubled isospins, but instead takes GMP fractions,
    so that half-integer values can be represented.
*/
sqrat su2_cgc(mpq_class I, mpq_class Iz, mpq_class i1, mpq_class i1z,
                mpq_class i2, mpq_class i2z)
{
    /* Eliminate common factors between the numerator and denominator of each
        argument (this is the most convenient place to do this)
    */
    I.canonicalize();
    Iz.canonicalize();
    i1.canonicalize();
    i1z.canonicalize();
    i2.canonicalize();
    i2z.canonicalize();

    return su2_cgc_2i(mpq_to_long(I*2), mpq_to_long(Iz*2),
                        mpq_to_long(i1*2), mpq_to_long(i1z*2),
                        mpq_to_long(i2*2), mpq_to_long(i2z*2));
}
