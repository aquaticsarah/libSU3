/* libSU3: External interface */

#ifndef __SU3_H__
#define __SU3_H__

#include <stddef.h>
#include <gmpxx.h>

/* Class to represent (+ or -) the square root of a rational.
    The actual value represented is sign(v) * |v|
*/
class sqrat
{
private:
    mpq_class v;

public:
    /* Component-wise constructors. sqrat(p, q) returns sign(pq) * sqrt(|p|/|q|) */
    sqrat(mpz_class, mpz_class);
    sqrat(long, long);

    /* Casts from other types. sqrat(x) returns a value which should equal x. */
    sqrat(mpq_class);
    sqrat(long);
    sqrat();

    /* In-place arithmetic */
    sqrat& operator+=(const sqrat&);
    sqrat& operator*=(const sqrat&);
    sqrat& operator/=(const sqrat&);
    sqrat& operator-=(const sqrat&);

    /* Arithmetic operations.
        Note that +, -, sqrt may throw std::domain_error if the output cannot be
        represented in a valid format (eg, for 1 + sqrt(2), which cannot
        be represented as a sqrat)
    */
    friend sqrat operator+(const sqrat&);
    friend sqrat operator-(const sqrat&);

    friend sqrat operator*(sqrat, const sqrat&);
    friend sqrat operator/(sqrat, const sqrat&);
    friend sqrat operator+(sqrat, const sqrat&);
    friend sqrat operator-(sqrat, const sqrat&);

    friend sqrat sqrt(const sqrat&);

    /* Comparisons */
    friend bool operator<(const sqrat&, const sqrat&);
    friend bool operator<=(const sqrat&, const sqrat&);
    friend bool operator==(const sqrat&, const sqrat&);
    friend bool operator!=(const sqrat&, const sqrat&);
    friend bool operator>(const sqrat&, const sqrat&);
    friend bool operator>=(const sqrat&, const sqrat&);

    /* Conversions to various types */
    char* tostring(char*, size_t);
    explicit operator double();
};

/* Classes to hold SU(3) isoscalar factors and Clebsch-Gordan coefficients.

    Note that it is easiest to calculate a whole degenerate set of reps at
    the same time, so we always return an array containing them all.
*/
class isoarray;
class cgarray;

/* A class to hold the isoscalar factors for a particular coupling */
class isoarray
{
    friend class cgarray;
    
private:
    isoarray* isf;
    size_t size; // Size of the following array
    sqrat* isf_array;

    void set_isf(long n, long k, long l, long k1, long l1,
                    long k2, long l2, sqrat v);

    /* Target and factor reps */
    const long p, q, p1, q1, p2, q2;

    /* Degeneracy of target rep */
    const long d;

public:
    /* Note: This type takes ownership of the array passed in - that is, it will
        delete the array when the isoarray object is deleted. */
    isoarray(long p, long q, long p1, long q1, long p2, long q2, long d,
                sqrat* isf_array);
    ~isoarray();

    /* We use operator() instead of operator[] as an easy way to use
        multiple indices */
    sqrat operator()(long n, long k, long l, long k1, long l1,
                        long k2, long l2);

    /* Convert to Clebsch-Gordans. This returns a newly-allocated cgarray object. */
    cgarray* to_cgarray();

    /* Apply the various symmetry relations */
    isoarray* exch_12();
    isoarray* exch_13bar();
    isoarray* exch_23bar();
};

/* A class to hold the Clebsch-Gordan coefficients for a particular coupling */
class cgarray
{
private:
    isoarray* isf;

public:
    /* Note: This type takes ownership of the isoarray object passed in */
    cgarray(isoarray* isf);
    ~cgarray();

    /* We use operator() instead of operator[] as an easy way to use
        multiple indices */
    sqrat operator()(long n, long k, long l, long m,
                        long k1, long l1, long m1,
                        long k2, long l2, long m2);

    /* Convert to ISFs. This returns a newly-allocated isoarray object. */
    isoarray* to_isoarray();

    /* Apply the various symmetry relations */
    cgarray* exch_12();
    cgarray* exch_13bar();
    cgarray* exch_23bar();
};

/* Calculate the dimension of one irrep */
long dimension(long p, long q);

/* Calculate the degeneracy of the (p,q) irrep in the decomposition of
    (p1,q1) x (p2,q2). Returns 0 if (p,q) is not a summand in this
    decomposition.
*/
long degeneracy(long p, long q, long p1, long q1, long p2, long q2);

/* Phase changes under the 1<->2 symmetry and the conjugation symmetry.
    Note: Kaeding and Williams call these phases xi_1 and xi_3 respectively.
*/
long phase_exch_12(long p, long q, long p1, long q1, long p2, long q2);
long phase_conj(long p, long q, long p1, long q1, long p2, long q2);

/* Calculate a single SU(2) Clebsch-Gordan coefficient.
    All arguments are implicitly doubled - eg, I represents
    2*(the actual isospin).
    This may throw std::domain_error if any of {I,i1,i2} are negative.
*/
sqrat su2_cgc_2i(long I, long Iz, long i1, long i1z,
                    long i2, long i2z);

/* Calculate a single SU(2) Clebsch-Gordan coefficient.
    This does *not* take doubled isospins, but instead takes GMP fractions,
    so that half-integer values can be represented.
    This throws std::domain_error if the inputs are not half-integers,
    and may do so if any {I,i1,i2} are negative.
*/
sqrat su2_cgc(mpq_class I, mpq_class Iz, mpq_class i1, mpq_class i1z,
                mpq_class i2, mpq_class i2z);

/* Main calculation functions.
    Note that if the calculation fails, these can throw std::logic_error.
    However, this should never happen unless there is a bug in the library.
*/
isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2);
cgarray* clebsch_gordans(long p, long q, long p1, long q1, long p2, long q2);

#endif
