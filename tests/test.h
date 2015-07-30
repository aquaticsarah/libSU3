/* libSU3: Test infrastructure */

#ifndef __SU3_TEST_H__
#define __SU3_TEST_H__

#include <stdio.h>

extern int tests_run;
extern int tests_passed;

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

#endif
