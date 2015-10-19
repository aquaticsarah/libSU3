# libSU3: Instructions for building various files

# Load configuration
include config.mk

# Autodetect source files and fill out relevant variables
include autodetect.mk

# Top-level targets
default: libSU3.a
debug: libSU3-debug.a
profile: libSU3-prof.a

su3: $(BUILDDIR)/progs/su3.o libSU3.a
	@echo "Linking su3"
	@$(LD) $(LDFLAGS) $(BUILDDIR)/progs/su3.o libSU3.a $(LIBRARIES) -o $@

test: run-tests
	@echo "Running tests"
	@./run-tests

bench: run-bench
	@./run-bench

clean:
	@echo "Cleaning up"
	@rm -rf $(BUILDDIR)/

clean-all:
	@echo "Cleaning up everything"
	@rm -rf $(BUILDDIR)/
	@rm -f libSU3*.a run-tests run-bench su3

# Intermediate files
libSU3.a libSU3-debug.a libSU3-prof.a:
	@echo "AR $@"
	@ar rcsu $@ $?

libSU3.a: $(OBJ)
libSU3-debug.a: $(DEBUG_OBJ)

run-tests: $(TEST_OBJ) $(TEST_RUNNER_OBJ) libSU3-debug.a
	@echo "Linking test driver"
	@$(LD) $(DEBUG_LDFLAGS) $^ $(LIBRARIES) -o $@

run-bench: $(BUILDDIR)/progs/bench.o libSU3.a
	@echo "Linking benchmark program"
	@$(LD) $(LDFLAGS) $^ $(LIBRARIES) -o $@

# The object files and dependency information can be done in a uniform way across
# the library and its tests. The only difference is which directories should be
# searched for include files.
$(BUILDDIR)/%.o: %.cc | $(DIRS)
	@echo "CC $<"
	@$(CC) $(CFLAGS) $(THIS_INCLUDE) $< -o $(BUILDDIR)/$*.o \
		-MMD -MQ $(BUILDDIR)/$*.o -MQ $(BUILDDIR)/$*.d -MF $(BUILDDIR)/$*.d

$(BUILDDIR)/debug/%.o: %.cc | $(DIRS)
	@echo "CC $< [DEBUG]"
	@$(CC) $(DEBUG_CFLAGS) $(THIS_INCLUDE) $< -o $(BUILDDIR)/debug/$*.o \
		-MMD -MQ $(BUILDDIR)/debug/$*.o -MQ $(BUILDDIR)/debug/$*.d -MF $(BUILDDIR)/debug/$*.d

$(BUILDDIR)/prof/%.o: %.cc | $(DIRS)
	@echo "CC $< [PROFILE]"
	@$(CC) $(PROFILE_CFLAGS) $(THIS_INCLUDE) $< -o $(BUILDDIR)/prof/$*.o \
		-MMD -MQ $(BUILDDIR)/prof/$*.o -MQ $(BUILDDIR)/prof/$*.d -MF $(BUILDDIR)/prof/$*.d

$(OBJ) $(DEBUG_OBJ) $(PROFILE_OBJ): THIS_INCLUDE=$(LIB_INCLUDE)
$(PROG_OBJ): THIS_INCLUDE=$(INCLUDE)
$(TEST_OBJ): THIS_INCLUDE=$(TEST_INCLUDE)

# The test runner file needs its own rules
$(TEST_RUNNER): scripts/gen_test_runner.py tests/*.cc | $(DIRS)
	@echo "Generating test runner ($@)..."
	@scripts/gen_test_runner.py $@

$(TEST_RUNNER_OBJ): $(TEST_RUNNER) | $(DIRS)
	@echo "CC $< [DEBUG]"
	@$(CC) $(DEBUG_CFLAGS) $(TEST_INCLUDE) $< -o $(TEST_RUNNER_OBJ) \
		-MMD -MQ $(TEST_RUNNER_OBJ) -MQ $(TEST_RUNNER_DEP) -MF $(TEST_RUNNER_DEP)

# Directory tree
$(DIRS):
	@mkdir -p $@

# Include dependency information
-include $(ALL_DEP)
