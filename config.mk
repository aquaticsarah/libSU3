# libSU3: Configuration for the build system

# Where to place intermediate files
BUILDDIR := /tmp/libSU3

# Toolchain
CC := g++ -c
LD := g++

# Include lists for library code and for the tests
INCLUDE := -I include/
LIB_INCLUDE := $(INCLUDE) -I src/
TEST_INCLUDE := $(INCLUDE) -I tests/

# External libraries to link against
LIBRARIES := -lgmpxx -lgmp

# Flags to pass to the build tools. The version with no prefix is passed
# when doing a regular build, the version with _DEBUG is passed for
# debug builds (which includes when running the tests).
COMMON_CFLAGS := -std=c++11 -Wall -Wextra -Werror
COMMON_LDFLAGS :=

CFLAGS := $(COMMON_CFLAGS) -DNDEBUG -O2
LDFLAGS := $(COMMON_LDFLAGS)

DEBUG_CFLAGS := $(COMMON_CFLAGS) -ggdb
DEBUG_LDFLAGS := $(COMMON_LDFLAGS) -ggdb

# Profile builds should be as close to normal builds as possible, just with
# an extra argument to the compiler and linker
PROFILE_CFLAGS := $(CFLAGS) -pg
PROFILE_LDFLAGS := $(LDFLAGS) -pg
