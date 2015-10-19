#!/usr/bin/env python
# libSU3: Script to auto-generate a C file, which runs each test in turn.

import re, glob, sys

# We piece together the output C file using the following templates.
# They are listed here in the same order they come in the output file.
HEADER = """\
/* libSU3: Main test driver.
    Auto-generated, do not modify.
*/

#include <stdio.h>
#include <limits.h>

#include "SU3.h"
#include "test.h"

int tests_run;
int tests_passed;

"""

# We need to forward-declare each test we are going to run.
DECLARE = "void test_%s();\n"

MIDDLE = """\

int main()
{
    srand(time(NULL));
    int total_tests_run = 0, total_tests_passed = 0;

"""

# Run each test and print its results
RUN = """\
    tests_run = 0;
    tests_passed = 0;
    fprintf(stderr, "\\nTesting %s...\\n");
    test_%s();
    fprintf(stderr, "%s: %%d/%%d tests passed\\n", tests_passed, tests_run);
    total_tests_run += tests_run;
    total_tests_passed += tests_passed;

"""

FOOTER = """\
    printf("\\nSummary: %d/%d tests passed\\n", total_tests_passed, total_tests_run);
    if (total_tests_passed < total_tests_run)
        return 1;
    else
        return 0;
}
"""

# Find all of the tests defined in a particular file
def scan_file(path):
    data = open(path).read()
    return re.findall("TEST\((\w+)\)", data)

def main(dest):
    all_tests = []

    # Collect up a list of all test functions
    for path in glob.glob("tests/*.cc"):
        all_tests.extend(scan_file(path))

    # Generate a test-runner file
    f = open(dest, "w")
    f.write(HEADER)

    for module in all_tests:
        f.write(DECLARE % module)

    f.write(MIDDLE)

    for module in all_tests:
        f.write(RUN % (module, module, module))

    f.write(FOOTER)

# Allow this program to be run as a script, which takes one argument
if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: %s DEST" % sys.argv[0])
        sys.exit(2)

    main(sys.argv[1])
