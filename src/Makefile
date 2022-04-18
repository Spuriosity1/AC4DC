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

CXXFLAGS = -std=c++17 -fopenmp -MD -g -Wall -pedantic
# LDFLAGS = -Wl,--gc-sections

mkfile_path := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Output directory after linking
BIN_PATH = ../bin
BUILD_PATH = ../build

# ------- System Include  --------------
# Change ~/include as necessary for your architecture
INC_PATH = -I/opt/homebrew/include


# -------- Establishing objects and sources ------

RATE_SRC = $(wildcard hartreefock/*.cpp) $(wildcard struct_hartreefock/*.cpp) $(wildcard numerical/*.cpp)
RATE_OBJ = $(addprefix $(BUILD_PATH)/,$(RATE_SRC:.cpp=.o)))

# ------------- Specific Include -----------

RATE_INC = $(wildcard hartreefock/*.hpp) $(wildcard struct_hartreefock/*.hpp) $(wildcard numerical/*.hpp)


# ---------- Libraries needed for linking ----------

# Change or remove ~/lib as necessary for your architecture
LIB_PATH = -L/opt/homebrew/lib 




# Generic .cc -> .o rule: note we generally don't have .hh files
$(BUILD_PATH)/%.o: %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INC_PATH) -MMD -MP $< -c -o $@


# ---------- Linking rules ----------

# Generic .o -> exec rule: uses all prerequisites (meant to be .o files)
ratecalc: $(RATE_OBJ)
	$(CXX) $(CXXFLAGS) $(INC_PATH) $(RATE_INC) $(LDFLAGS) -MMD $^ \
	$(LIB_PATH) -o $(BIN_PATH)/$@

# TODO
# simulate: 

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

