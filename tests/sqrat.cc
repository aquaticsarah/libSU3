/* libSU3: Tests for the 'sqrat' type */

#include "SU3.h"
#include "test.h"

TEST(sqrat)
{
    /* Test reduction */
    sqrat a(2, 4);
    TEST_EQ_SQRAT(a, 1, 2, "reduce(sqrt(2/4))");

    a = sqrat(3, -3);
    TEST_EQ_SQRAT(a, -1, 1, "reduce(sqrt(3/-3))");

    a = sqrat(-11);
    TEST_EQ_SQRAT(a, -121, 1, "-11");
}

TEST(sqrat_arithmetic)
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

    /* Square-roots */
    c = sqrt(sqrat(1, 9));
    TEST_EQ_SQRAT(c, 1, 3, "sqrt(sqrt(1/9))");
}
