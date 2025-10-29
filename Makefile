CXXFLAGS  += -I./include -I./include/star -DHAVE_CONFIG_H -DTOOLBOX -Wconversion
LDLIBS	  += -Lsstar -lsstar -Lnode -lnode -Lstd -lstd -lm
CFLAGS    += -O
CXXFLAGS  += -O

EXE	= starev
#DIRS	= node std sstar dstar
DIRS	= node std sstar dstar rdc

all:	lib $(EXE)

test-compiler:
	echo $(CXX)
	echo `which $(CXX)`

Makefile.inc:	Makefile.inc.conf
	@sed s#__BASE__#$(PWD)#g Makefile.inc.conf >Makefile.inc

lib:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done

clean:	Makefile.inc
	@for d in $(DIRS) ; do echo '\nmake' $@ in $$d; make -C $$d $@ ; done
	/bin/rm -f *~ $(EXE) Makefile.inc

# --- Test support ---
CXXFLAGS  += -std=c++23
TESTDIR   = tests
BUILDDIR  = build
TESTBIN   = $(BUILDDIR)/run_tests
TESTSRC   = $(wildcard $(TESTDIR)/*.cpp)
TESTOBJ   = $(patsubst $(TESTDIR)/%.cpp,$(BUILDDIR)/%.o,$(TESTSRC))

GTEST_DIR = third_party/googletest/googletest
GTEST_SRC = $(GTEST_DIR)/src/gtest-all.cc
GTEST_OBJ = $(BUILDDIR)/gtest-all.o
GTEST_INC = -I$(GTEST_DIR)/include -I$(GTEST_DIR)

# Ensure build directory exists
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Compile gtest from source
$(GTEST_OBJ): $(GTEST_SRC) | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -c $< -o $@

# Compile test sources
$(BUILDDIR)/%.o: $(TESTDIR)/%.cpp | $(BUILDDIR)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -I./include -c $< -o $@

# Link final test binary
$(TESTBIN): $(TESTOBJ) $(GTEST_OBJ)
	$(CXX) $(CXXFLAGS) $(GTEST_INC) -I./include $^ $(LDLIBS) -lpthread -o $@

# Run tests
test: $(TESTBIN)
	@echo "Running unit tests..."
	@./$(TESTBIN)

# Clean test artifacts
clean-tests:
	/bin/rm -rf $(BUILDDIR)
