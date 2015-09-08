# libSU3: Instructions for building various files

# Where to place intermediate files
BUILDDIR := /tmp/libSU3

# Autodetect source files
SRC := $(wildcard src/*.cc)
TEST_SRC := $(wildcard tests/*.cc)

# Various files which we will produce
OBJ := $(SRC:src/%.cc=$(BUILDDIR)/%.o)
TEST_OBJ := $(TEST_SRC:tests/%.cc=$(BUILDDIR)/tests/%.o)

DEP := $(SRC:src/%.cc=$(BUILDDIR)/%.d)
TEST_DEP := $(TEST_SRC:tests/%.cc=$(BUILDDIR)/tests/%.d)

TEST_RUNNER := $(BUILDDIR)/tests/run.cc
TEST_RUNNER_OBJ := $(BUILDDIR)/tests/run.o
TEST_RUNNER_DEP := $(BUILDDIR)/tests/run.d

BUILD_FILES := $(OBJ) $(DEP)
TEST_FILES := $(TEST_OBJ) $(TEST_DEP) $(TEST_RUNNER) $(TEST_RUNNER_OBJ) $(TEST_RUNNER_DEP)

# Only build the files we need
ifeq ($(findstring test,$(MAKECMDGOALS)),)
FILES := $(BUILD_FILES)
else
FILES := $(BUILD_FILES) $(TEST_FILES)
endif

DIRS := $(sort $(dir $(FILES)))

# Toolchain
CC := g++ -c
CFLAGS := -Wall -Wextra -Werror -DNDEBUG -ggdb

# Include lists for library code and for the tests
# (each of which has internal headers which shouldn't be used
#  by the other)
INCLUDE := -I include/
LIB_INCLUDE := $(INCLUDE) -I src/
TEST_INCLUDE := $(INCLUDE) -I tests/

# External libraries to link against
LIBRARIES := -lgmpxx -lgmp

LD := g++
LDFLAGS := -ggdb

MKDEP := g++ -MM -MP

# Top-level targets
default: libSU3.a

test: run-tests
	@echo "Running tests"
	@./run-tests

clean:
	@echo "Cleaning up"
	@rm -rf $(BUILDDIR)/

clean-all:
	@echo "Cleaning up everything"
	@rm -rf $(BUILDDIR)/
	@rm -f libSU3.a run-tests

libSU3.a: $(OBJ)
	@echo "AR $@"
	@ar rcsu $@ $(OBJ)

# Intermediate targets
run-tests: $(TEST_OBJ) $(TEST_RUNNER_OBJ) libSU3.a
	@echo "Linking test driver"
	@$(LD) $(LDFLAGS) $(TEST_OBJ) $(TEST_RUNNER_OBJ) libSU3.a $(LIBRARIES) -o $@

# Rules to build object files and dependency information
$(OBJ):$(BUILDDIR)/%.o: src/%.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(LIB_INCLUDE) $< -o $@

$(DEP):$(BUILDDIR)/%.d: src/%.cc | $(DIRS)
	@echo "MKDEP $<"
	@$(MKDEP) $(LIB_INCLUDE) $< -MQ $(BUILDDIR)/$*.o -MQ $@ -MF $@

$(TEST_OBJ):$(BUILDDIR)/tests/%.o: tests/%.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(TEST_INCLUDE) $< -o $@

$(TEST_DEP):$(BUILDDIR)/tests/%.d: tests/%.cc | $(DIRS)
	@echo "MKDEP $<"
	@$(MKDEP) $(TEST_INCLUDE) $< -MQ $(BUILDDIR)/tests/$*.o -MQ $@ -MF $@

# The test runner file needs its own rules
$(TEST_RUNNER): scripts/gen_test_runner.py tests/*.cc | $(DIRS)
	@echo "Generating test runner ($@)..."
	@scripts/gen_test_runner.py $@

$(TEST_RUNNER_OBJ): $(TEST_RUNNER) | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(TEST_INCLUDE) $< -o $@

$(TEST_RUNNER_DEP): $(TEST_RUNNER) | $(DIRS)
	@echo "MKDEP $<"
	@$(MKDEP) $(TEST_INCLUDE) $< -MQ $(BUILDDIR)/tests/$*.o -MQ $@ -MF $@

# Directory tree
$(DIRS):
	@mkdir -p $@

# Include dependency information, but only if relevant
ifeq ($(findstring clean,$(MAKECMDGOALS)),)
-include $(DEP)
endif

ifneq ($(findstring test,$(MAKECMDGOALS)),)
-include $(TEST_DEP) $(TEST_RUNNER_DEP)
endif
