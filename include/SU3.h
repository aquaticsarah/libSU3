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
    sqrat(mpq_class v);
    sqrat(long, long);
    sqrat(mpz_class, mpz_class);
    sqrat(long);
    sqrat();

    /* In-place arithmetic */
    sqrat& operator+=(const sqrat&);
    sqrat& operator*=(const sqrat&);
    sqrat& operator/=(const sqrat&);
    sqrat& operator-=(const sqrat&);

    /* Arithmetic operations */
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

   Notes on indexing:
   - For each of the three reps involved, we have the following ranges:
      q <= k <= p+q (for a total of p+1 possible values of k)
      0 <= l <= q   (for a total of q+1 possible values of l)
     The value of 'l' can be used as an index directly, but that for
     'k' needs to be shifted down by q.

   - All reps in a degenerate set need to be processed at once, so we
     allocate storage all at once.

   Implementation details:
   - Sometimes, when calculating ISFs, we apply a recursion relation involving
     values which are one space "off the edge" of the valid range (eg, with
     l1=-1). In order to simplify this code, we allow indexes of that form,
     but *only* for the internal-only isoscalar_context class.

   - Given k, l, k1, l1, k2, there is a unique valid value of l2 determined
     by hypercharge conservation. As such, we can save a factor of (q2+1) on
     memory space by not allocating an l2 axis.
     However, we do accept l2 as an argument and, unless -DNDEBUG is specified
     when compiling the library, we check that it is valid.
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

/* Get the name of a representation. */
char* repname(char* buffer, size_t len, long p, long q);

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
*/
sqrat su2_cgc_2i(long I, long Iz, long i1, long i1z,
                    long i2, long i2z);

/* Calculate a single SU(2) Clebsch-Gordan coefficient.
    This does *not* take doubled isospins, but instead takes GMP fractions,
    so that half-integer values can be represented.
*/
sqrat su2_cgc(mpq_class I, mpq_class Iz, mpq_class i1, mpq_class i1z,
                mpq_class i2, mpq_class i2z);

/* Main calculation functions. */
isoarray* isoscalars(long p, long q, long p1, long q1, long p2, long q2);
cgarray* clebsch_gordans(long p, long q, long p1, long q1, long p2, long q2);

#endif
