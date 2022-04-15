
CPP := g++-11

LIB := -fopenmp

INC := -Iinclude -I/opt/homebrew/include

SRCDIR := src
BUILDDIR := build
TARGET := bin/ac4dc

SRCEXT := cpp

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT) )

OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))
CXXFLAGS := -std=c++17 -fopenmp -MD -g -Wall

# debug: CXXFLAGS += -DDEBUG -Wpedantic
# release: CXXFLAGS += -O3 -DNDEBUG

$(TARGET): $(OBJECTS)
	@echo " Linking $(TARGET)... ";
	@echo " $(CPP) $^ $(LIB) -o $(TARGET) "; $(CPP) $^ $(LIB) -o $(TARGET)

all: $(TARGET)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR) bin $(BUILDDIR)/Wigner $(BUILDDIR)/$(MAINSUBDIR)
	@echo " $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<"; $(CPP) $(CXXFLAGS) $(INC) -c -o $@ $<


debug: all
release: all

scrub:
	$(RM) build/FreeDistribution.o build/ElectronSolver.o build/SplineBasis.o build/SplineIntegral.o

clean:
	@echo " Cleaning...";
	@echo " $(RM) -r $(BUILDDIR)"; $(RM) -r $(BUILDDIR)

.PHONY: clean, test, all
-include $(OBJECTS:.o=.d)
