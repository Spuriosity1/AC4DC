/**
 * @file Plotting.h
 * @brief Defines the class Plotting, which handles the creation of _live_plot.png
*/

/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
// (C) Spencer Passmore 2023

// https://stackoverflow.com/questions/8293775/live-graph-for-a-c-application
// https://stackoverflow.com/questions/3286448/calling-a-python-method-from-c-c-and-extracting-its-return-value
// https://stackoverflow.com/questions/28269157/plotting-in-a-non-blocking-way-with-matplotlib


#pragma once

#define PY_SSIZE_T_CLEAN
//#include <Python.h>
#include <pybind11/embed.h>
#include <pybind11/stl.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
//#include <dlfcn.h>


namespace py = pybind11;


struct Plotting{    
    //py::scoped_interpreter guard{};
    Plotting();
    void plot_frame(std::vector<double> energies, std::vector<double> density,std::vector<double> knot_to_plot);
};
