
CPP := g++-9

LIB := -fopenmp

INC := -I$(HOME)/Programming/include -Iinclude

SRCDIR := src
BUILDDIR := build
TARGET := bin/ac4dc2

SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT) )

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CXXFLAGS := -std=c++11 -fopenmp

debug: CXXFLAGS += -DDEBUG -g
release: CXXFLAGS += -O3 -DNDEBUG

# For tests:
TESTS := $(shell find tests -type f -name *.$(SRCEXT) )
TINC := $(INC) -Isrc

$(TARGET): $(OBJECTS)
	@echo " Linking $(TARGET)... ";
	@echo " $(CPP) $^ $(LIB) -o $(TARGET) "; $(CPP) $^ $(LIB) -o $(TARGET)

all: $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin $(BUILDDIR)/Wigner $(BUILDDIR)/$(MAINSUBDIR)
	@echo " $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<

test: $(TESTS)
	@mkdir -p bin/tests tmp testOutput
	@echo " Compiling tests/abm_verif.cpp "
	# $(CPP) -g -std=c++11 tests/abm_verif.cpp $(TINC) -o bin/tests/abm_verif
	# $(CPP) -g tests/binsearch.cpp -o bin/tests/binsearch
	# $(CPP) -g tests/spline_check.cpp src/SplineBasis.cpp $(TINC) -o bin/tests/splinecheck
	# $(CPP) -g -std=c++11 tests/q_ee.cpp src/RateSystem.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp $(TINC) -o bin/tests/q_ee
	$(CPP) -g $(CXXFLAGS) tests/q_eii.cpp src/Constant.cpp src/FreeDistribution.cpp src/SplineIntegral.cpp src/SplineBasis.cpp src/Dipole.cpp $(TINC) -o bin/tests/q_eii

debug: all
release: all

scrub:
	$(RM) build/FreeDistribution.o build/ElectronSolver.o build/SplineBasis.o build/SplineIntegral.o

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR)"; $(RM) -r $(BUILDDIR)

.PHONY: clean, test, all
