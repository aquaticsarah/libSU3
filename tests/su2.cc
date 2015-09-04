/* libSU3: Tests for SU(2) Clebsch-Gordans */

#include "SU3.h"
#include "test.h"

#define TEST_SU2_CG(I, Iz, i1, i1z, i2, i2z, p, q) \
    TEST_EQ_SQRAT(su2_cgc(I, Iz, i1, i1z, i2, i2z), p, q, \
        "<%Qd,%Qd| (|%Qd,%Qd> x |%Qd,%Qd>)", \
        mpq_class(I).get_mpq_t(), mpq_class(Iz).get_mpq_t(), \
        mpq_class(i1).get_mpq_t(), mpq_class(i1z).get_mpq_t(), \
        mpq_class(i2).get_mpq_t(), mpq_class(i2z).get_mpq_t())

TEST(su2)
{
    /* Test values for 1/2 x 1/2 -> 0 + 1 */
    mpq_class half = mpq_class(1, 2);
    mpq_class three_half = mpq_class(3, 2);

    //           I  Iz    i1    i1z    i2    i2z      p  q
    TEST_SU2_CG(0,  0, half, -half, half,  half,    -1, 2);
    TEST_SU2_CG(0,  0, half,  half, half, -half,     1, 2);

    TEST_SU2_CG(1, -1, half, -half, half, -half,     1, 1);
    TEST_SU2_CG(1,  0, half, -half, half,  half,     1, 2);
    TEST_SU2_CG(1,  0, half,  half, half, -half,     1, 2);
    TEST_SU2_CG(1,  1, half,  half, half,  half,     1, 1);

    /* Test values for 1 x 1/2 -> 1/2 + 3/2 */
    //              I     Iz i1 i1z    i2    i2z      p  q
    TEST_SU2_CG(half, -half, 1, -1, half,  half,    -2, 3);
    TEST_SU2_CG(half, -half, 1,  0, half, -half,     1, 3);
    TEST_SU2_CG(half,  half, 1,  0, half,  half,    -1, 3);
    TEST_SU2_CG(half,  half, 1,  1, half, -half,     2, 3);

    //                    I           Iz i1 i1z    i2    i2z      p  q
    TEST_SU2_CG(three_half, -three_half, 1, -1, half, -half,     1, 1);
    TEST_SU2_CG(three_half,       -half, 1, -1, half,  half,     1, 3);
    TEST_SU2_CG(three_half,       -half, 1,  0, half, -half,     2, 3);
    TEST_SU2_CG(three_half,        half, 1,  0, half,  half,     2, 3);
    TEST_SU2_CG(three_half,        half, 1,  1, half, -half,     1, 3);
    TEST_SU2_CG(three_half,  three_half, 1,  1, half,  half,     1, 1);

    /* Test  values for 1 x 1 -> 0 + 1 + 2 */
    //           I  Iz i1 i1z i2 i2z      p  q
    TEST_SU2_CG(0,  0, 1, -1, 1,  1,     1, 3);
    TEST_SU2_CG(0,  0, 1,  0, 1,  0,    -1, 3);
    TEST_SU2_CG(0,  0, 1,  1, 1, -1,     1, 3);

    TEST_SU2_CG(1, -1, 1, -1, 1,  0,    -1, 2);
    TEST_SU2_CG(1, -1, 1,  0, 1, -1,     1, 2);
    TEST_SU2_CG(1,  0, 1, -1, 1,  1,    -1, 2);
    TEST_SU2_CG(1,  0, 1,  0, 1,  0,     0, 1);
    TEST_SU2_CG(1,  0, 1,  1, 1, -1,     1, 2);
    TEST_SU2_CG(1,  1, 1,  1, 1,  0,     1, 2);
    TEST_SU2_CG(1,  1, 1,  0, 1,  1,    -1, 2);

    TEST_SU2_CG(2, -2, 1, -1, 1, -1,     1, 1);
    TEST_SU2_CG(2, -1, 1, -1, 1,  0,     1, 2);
    TEST_SU2_CG(2, -1, 1,  0, 1, -1,     1, 2);
    TEST_SU2_CG(2,  0, 1, -1, 1,  1,     1, 6);
    TEST_SU2_CG(2,  0, 1,  0, 1,  0,     2, 3);
    TEST_SU2_CG(2,  0, 1,  1, 1, -1,     1, 6);
    TEST_SU2_CG(2,  1, 1,  1, 1,  0,     1, 2);
    TEST_SU2_CG(2,  1, 1,  0, 1,  1,     1, 2);
    TEST_SU2_CG(2,  2, 1,  1, 1,  1,     1, 1);
}
