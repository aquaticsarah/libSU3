# Autodetect source files
SRC := $(wildcard src/*.cc)
TEST_SRC := $(wildcard tests/*.cc)

# Various files which we will produce
OBJ := $(SRC:src/%.cc=build/%.o)
TEST_OBJ := $(TEST_SRC:tests/%.cc=build/tests/%.o)

DEP := $(SRC:src/%.cc=build/%.d)
TEST_DEP := $(TEST_SRC:tests/%.cc=build/tests/%.d)

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

test: build/tests/main
	@echo "Running tests"
	@build/tests/main

clean:
	@echo "Cleaning up"
	@rm -rf build/

clean-all:
	@echo "Cleaning up everything"
	@rm -rf build/
	@rm -f libSU3.a

libSU3.a: $(OBJ)
	@echo "AR $@"
	@ar rcsu $@ $(OBJ)

# Intermediate targets
build/tests/main: $(TEST_OBJ) libSU3.a
	@echo "Linking test driver"
	@$(LD) $(LDFLAGS) $(TEST_OBJ) libSU3.a -o $@

# Rules to build object files and dependency information
$(OBJ):build/%.o: src/%.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(INCLUDE) $< -o $@

$(DEP):build/%.d: src/%.cc | $(DIRS)
	@echo "MKDEP $<"
	@$(MKDEP) $(INCLUDE) $< -MQ build/$*.o -MQ $@ -MF $@

$(TEST_OBJ):build/tests/%.o: tests/%.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(TEST_INCLUDE) $< -o $@

$(TEST_DEP):build/tests/%.d: tests/%.cc | $(DIRS)
	@echo "MKDEP $<"
	@$(MKDEP) $(TEST_INCLUDE) $< -MQ build/tests/$*.o -MQ $@ -MF $@

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
