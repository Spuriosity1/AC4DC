# Makefile to compile and link QSI simulator programs
#
# Created on 01/08/2018
# Copyright (C) 2018 Attila Szabó <as2372@cam.ac.uk>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the Creative Commons Attribution License (CC-BY),
# version 4.0 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# Please find a copy of the Creative Commons CC-BY License on
# <https://creativecommons.org/licenses/by/4.0/>.

# ---------- Compiler and linker directives ----------
CXX = g++-11

CXXFLAGS = -std=c++17 -pedantic -g -O1 -Wall -no-pie 
# LDFLAGS = -Wl,--gc-sections

mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Output directory after linking
BIN_PATH = ../bin
BUILD_PATH = ../build

# ------ Establishing objects and sources ----

# Get all of the .cpp files, however nasty
SOURCES = $(wildcard *.cpp) $(wildcard */*.cpp)

# Change ~/include as necessary for your architecture
INC_PATH = -I /opt/homebrew/include




INC_PATH += -I $(mkfile_path)$(COMMON_DIR) -I $(mkfile_path)$(STRUCT_DIR) -I $(mkfile_path)$(MANAGER_DIR) 

# ---------- Libraries needed for linking ----------

# Change or remove ~/lib as necessary for your architecture
LIB_PATH = -L /opt/homebrew/lib

LIBS_STD = -lstdc++ -lm
LIBS_GSL = -lgsl -lgslcblas
LIBS_FFT = -lfftw3


# ---------- Compilation rules ----------

# Nontrivial header dependences
DEPENDS := $(patsubst %.cc,%.d,$(SOURCES))



# Generic .cc -> .o rule: note we generally don't have .hh files
$(BUILD_PATH)/%.o: %.cc
	mkdir -p $(BUILD_PATH)/$(COMMON_DIR)
	mkdir -p $(BUILD_PATH)/$(MANAGER_DIR)
	mkdir -p $(BUILD_PATH)/$(STRUCT_DIR)
	$(CXX) $(CXXFLAGS) $(INC_PATH) -MMD -MP $< -c -o $@




# ---------- Linking rules ----------

# Generic .o -> exec rule: uses all prerequisites (meant to be .o files)
%: $(BUILD_PATH)/%.o $(COMMON_O) $(MANAGER_O) $(STRUCT_O)
	$(CXX) $(CXXFLAGS) $(INC_PATH) $(LDFLAGS) -MMD $^ \
	$(LIB_PATH) $(LIBS_STD) $(LIBS_GSL) $(LIBS_FFT) \
	-o $(BIN_PATH)/$@

# Don't allow direct compilation 
%: %.cc
# Don't allow removal of object files in chained implicit linking
.SECONDARY:

# Cleanup:
.PHONY: clean build_directories
clean:
	-rm -r $(BUILD_PATH)/**
print:
	@echo $(COMMON_O)
	@echo $(MANAGER_O)
	@echo $(STRUCT_O)

