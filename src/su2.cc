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
        throw std::domain_error("Trying to convert non-integer fraction "
                                "to long");

    return num.get_si();
}

/* Wrappers to allow us to do min/max of MPQs.
    Note that, in this file, we *only* care about MPQs which are
    actually integers
*/
static long min(mpq_class a, mpq_class b, mpq_class c)
{
    return min(mpq_to_long(a), mpq_to_long(b), mpq_to_long(c));
}

static long max(mpq_class a, mpq_class b, mpq_class c)
{
    return max(mpq_to_long(a), mpq_to_long(b), mpq_to_long(c));
}

static mpz_class factorial(mpq_class x)
{
    long v = mpq_to_long(x);

    if (v < 0) throw std::domain_error("Factorial of a negative number");

    mpz_class res;
    mpz_fac_ui(res.get_mpz_t(), v);
    return res;
}

/* Calculate a single SU(2) Clebsch-Gordan coefficient.
    The expression used is from Bohm - see README for citation */
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

    /* Only states obeying the following conditions can couple:
        Iz = i1z + i2z
        i1 + i2 >= I >= |i1 - i2|
       Note |i1 - i2| = max(i1-i2, i2-i1), and so we can rewrite the last
       inequality as (I >= i1 - i2) && (I >= i2 - i1)
    */
    if (Iz != i1z + i2z) return sqrat(0);
    if ((I > i1 + i2) || (I < i1 - i2) || (I < i2 - i1)) return sqrat(0);

    mpz_class numerator = mpz_class(2*I + 1) * factorial(I + i1 - i2)
                        * factorial(I - i1 + i2) * factorial(i1 + i2 - I)
                        * factorial(i1 + i1z) * factorial(i1 - i1z)
                        * factorial(i2 + i2z) * factorial(i2 - i2z)
                        * factorial(I + Iz) * factorial(I - Iz);
    mpz_class denominator = factorial(I + i1 + i2 + 1);
    sqrat prefactor = sqrat(numerator, denominator);

    sqrat sum = 0;
    long zmin = max(0, i2 - i1z - I, i1 + i2z - I);
    long zmax = min(i1 + i2 - I, i1 - i1z, i2 + i2z);
    long z;

    for (z = zmin; z <= zmax; ++z)
    {
        denominator = factorial(z) * factorial(i1 + i2 - I - z)
                    * factorial(i1 - i1z - z) * factorial(i2 + i2z - z)
                    * factorial(I - i2 + i1z + z) * factorial(I - i1 - i2z + z);

        /* zhe summand is now ((-1)^z / denominator), *without* a square root.
            So we have to correct for that fact.
        */
        sum += sqrat(SIGN(z), denominator * denominator);
    }

    return prefactor * sum;
}
