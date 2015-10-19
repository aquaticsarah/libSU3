/* libSU3: Test infrastructure */

#ifndef __SU3_TEST_H__
#define __SU3_TEST_H__

#include <stdio.h>
#include <stdlib.h>

extern int tests_run;
extern int tests_passed;

/* Macro for generating random values in [0,n) */
#define RANDRANGE(n) ((long)((random() / (double)RAND_MAX) * n))

/* A label to indicate which functions are top-level tests.
    'module' should be the name of the module being tested (not in quotes) */
#define TEST(module) void test_##module()

/* Test some predicate, printing a message if it fails */
#define DO_TEST(pred, message, args...) \
    do \
    { \
        ++tests_run; \
        if (! (pred)) \
        { \
            fprintf(stderr, "%s:%d:", __FILE__, __LINE__); \
            fprintf(stderr, message, ## args); \
            fprintf(stderr, "\n"); \
        } \
        else \
            ++tests_passed; \
    } while (0)

/* Helper for tests involving the sqrat type
    (tests if a == sqrat(p,q)). "thing" is what the calculation
    is supposed to represent */
#define TEST_EQ_SQRAT(a, p, q, thing, args...) \
    do \
    { \
        char _abuf[64], _pqbuf[64], _thingbuf[64]; \
        (a).tostring(_abuf, 64); \
        sqrat _expected = sqrat(p, q); \
        _expected.tostring(_pqbuf, 64); \
        gmp_snprintf(_thingbuf, 64, thing, ## args); \
        DO_TEST((a == _expected), "Expected %s == %s, got %s", \
                _thingbuf, _pqbuf, _abuf); \
    } while (0)

#endif
