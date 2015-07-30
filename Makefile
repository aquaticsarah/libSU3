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

DIRS := $(sort $(dir $(OBJ) $(DEP) $(TEST_OBJ) $(TEST_DEP)))

# Toolchain
CC := g++ -c
CFLAGS := -Wall -Wextra -Werror
INCLUDE := -I include/
TEST_INCLUDE := $(INCLUDE) -I tests/

LD := g++
LDFLAGS :=

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
run-tests: $(TEST_OBJ) libSU3.a
	@echo "Linking test driver"
	@$(LD) $(LDFLAGS) $(TEST_OBJ) libSU3.a -o $@

# Rules to build object files and dependency information
$(OBJ):$(BUILDDIR)/%.o: src/%.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

$(DEP):$(BUILDDIR)/%.d: src/%.cc | $(DIRS)
	@echo "MKDEP $<"
	@$(MKDEP) $(INCLUDE) $< -MQ $(BUILDDIR)/$*.o -MQ $@ -MF $@

$(TEST_OBJ):$(BUILDDIR)/tests/%.o: tests/%.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(TEST_INCLUDE) $< -o $@

$(TEST_DEP):$(BUILDDIR)/tests/%.d: tests/%.cc | $(DIRS)
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
-include $(TEST_DEP)
endif
