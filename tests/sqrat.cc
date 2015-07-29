/* libSU3: Tests for the 'sqrat' type */

#include "SU3.h"
#include "test.h"

/* Helper for tests involving the sqrat type
   (tests if a == sqrat(p,q)). "thing" is what the calculation
   is supposed to represent */
#define TEST_EQ_SQRAT(a, p, q, thing) \
    do \
    { \
        char abuf[64], pqbuf[64]; \
        a.tostring(abuf, 64); \
        sqrat(p,q).tostring(pqbuf, 64); \
        DO_TEST((a.numerator() == p) && (a.denominator() == q), \
            "Expected %s == %s, got %s", thing, pqbuf, abuf); \
    } while (0)

TEST(gcd) void test_gcd()
{
    /* Test GCDs - in particular, check that the result always has the same
        sign as its second argument */
    long d1 = gcd(2, 5);
    long d2 = gcd(10, 5);
    long d3 = gcd(-22, 14);
    long d4 = gcd(3, -10101);
    long d5 = gcd(-55, -89);

    DO_TEST(d1 == 1, "Expected gcd(2,5)=1, got %ld", d1);
    DO_TEST(d2 == 5, "Expected gcd(5,10)=5, got %ld", d2);
    DO_TEST(d3 == 2, "Expected gcd(-22,14)=2, got %ld", d3);
    DO_TEST(d4 == -3, "Expected gcd(3,-10101)=3, got %ld", d4);
    DO_TEST(d5 == -1, "Expected gcd(-55,-89)=1, got %ld", d5);
}

TEST(sqrat) void test_sqrat()
{
    /* Test reduction */
    sqrat a(2, 4);
    TEST_EQ_SQRAT(a, 1, 2, "reduce(sqrt(2/4))");

    a = sqrat(3, -3);
    TEST_EQ_SQRAT(a, -1, 1, "reduce(sqrt(3/-3))");

    a = sqrat(-11);
    TEST_EQ_SQRAT(a, -121, 1, "-11");
}

TEST(sqrat-arithmetic) void test_sqrat_arith()
{
    sqrat a(1, 4);
    sqrat b(2, 3);
    sqrat c;
    sqrat d(-2, 3);

    /* Test multiplication */
    c = a*b;
    TEST_EQ_SQRAT(c, 1, 6, "sqrt(1/4) * sqrt(2/3)");

    c = a*2;
    TEST_EQ_SQRAT(c, 1, 1, "sqrt(1/4) * 2");

    c = a*(-3);
    TEST_EQ_SQRAT(c, -9, 4, "sqrt(1/4) * -3");

    /* Test division */
    c = a/b;
    TEST_EQ_SQRAT(c, 3, 8, "sqrt(1/4) / sqrt(2/3)");

    c = d/4;
    TEST_EQ_SQRAT(c, -1, 24, "-sqrt(2/3) / 4)");

    /* Test integer addition/subtraction */
    c = a + 1;
    TEST_EQ_SQRAT(c, 9, 4, "sqrt(1/4) + 1");

    c = a + (-1);
    TEST_EQ_SQRAT(c, -1, 4, "sqrt(1/4) + (-1)");

    c = a - 2;
    TEST_EQ_SQRAT(c, -9, 4, "sqrt(1/4) - 2");

    c = a - (-3);
    TEST_EQ_SQRAT(c,  49, 4, "sqrt(1/4) - (-3)");

    /* Test other addition/subtraction */
    c = b + d;
    TEST_EQ_SQRAT(c, 0, 1, "sqrt(2/3) + -sqrt(2/3)");

    c = d - b;
    TEST_EQ_SQRAT(c, -8, 3, "-sqrt(2/3) - sqrt(2/3)");
}
