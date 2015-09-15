# libSU3: Source file autodetection

SRC := $(wildcard src/*.cc)
PROG_SRC := $(wildcard progs/*.cc)
TEST_SRC := $(wildcard tests/*.cc)
TEST_RUNNER := $(BUILDDIR)/debug/test.cc

# Lists of files to produce
OBJ := $(SRC:%.cc=$(BUILDDIR)/%.o)
PROG_OBJ := $(PROG_SRC:%.cc=$(BUILDDIR)/%.o)
DEBUG_OBJ := $(SRC:%.cc=$(BUILDDIR)/debug/%.o)
PROFILE_OBJ := $(SRC:%.cc=$(BUILDDIR)/prof/%.o)
TEST_OBJ := $(TEST_SRC:%.cc=$(BUILDDIR)/debug/%.o)
TEST_RUNNER_OBJ := $(BUILDDIR)/debug/test.o

DEP := $(SRC:%.cc=$(BUILDDIR)/%.d)
PROG_DEP := $(PROG_SRC:%.cc=$(BUILDDIR)/%.d)
DEBUG_DEP := $(SRC:%.cc=$(BUILDDIR)/debug/%.d)
PROFILE_DEP := $(SRC:%.cc=$(BUILDDIR)/prof/%.d)
TEST_DEP := $(TEST_SRC:%.cc=$(BUILDDIR)/debug/%.d)
TEST_RUNNER_DEP := $(BUILDDIR)/debug/test.d

ALL_OBJ := $(OBJ) $(PROG_OBJ) $(DEBUG_OBJ) $(PROFILE_OBJ) $(TEST_OBJ) $(TEST_RUNNER_OBJ)

ALL_DEP := $(DEP) $(PROG_DEP) $(DEBUG_DEP) $(PROFILE_DEP) $(TEST_DEP) $(TEST_RUNNER_DEP)

DIRS := $(sort $(dir $(ALL_OBJ) $(ALL_DEP)))
