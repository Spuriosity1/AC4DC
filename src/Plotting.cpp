/**
 * @file Plotting.cpp
 * @brief 
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

#include <Plotting.h>
#include <Display.h>
#include <filesystem>
#include <iostream>
#include "config.h"

Plotting::Plotting(){ 
    #ifdef PYBIND
    //py::object current_path = PySys_GetObject("path");//std::getenv("PYTHONPATH");  // https://stackoverflow.com/questions/7137412/python-c-api-modify-search-path
    //string python_path = string(std::filesystem::current_path()) + "/scripts";    
    //python_path += current_path;   
    //std::wstring wide_py_path = std::wstring(python_path.begin(), python_path.end());
    //PySys_SetPath(wide_py_path.c_str());  
    try {
        PyRun_SimpleString(
        "import sys\n"
        "import os\n"
        "script_path = os.getcwd() + '/scripts' \n"
        "print('[ Plotting ] adding scripts directory to python path:', script_path)\n"
        "sys.path.append(script_path)\n"
        );
    }
    catch (const std::exception& e){
        Display::close();
        std::cerr << "Error: " << e.what() << std::endl;
        py::module_ sys = py::module_::import("sys");
        py::print("Python path:",sys.attr("path"));
        throw std::runtime_error("plotting failure1. See above for error message from python.");
    }
    #endif //PYBIND

}
void Plotting::plot_frame(std::vector<double> energies, std::vector<double> density,std::vector<double> knot_to_plot){
    #ifdef PYBIND
    py::list energy_list = py::cast(energies);
    py::list density_list = py::cast(density);    
    py::list knot_to_plot_list = py::cast(knot_to_plot); 
    //auto py_path = PySys_GetObject("path");
    py::module_ sys = py::module_::import("sys");
    try { 
        py::module_::import("scipy"); // test importing
        py::object py_clear_frame = py::module_::import("script_generate_frame").attr("clear_frame_c")();
        py::object py_plot_frame = py::module_::import("script_generate_frame").attr("plot_frame_from_c")(energy_list,density_list,knot_to_plot_list);
    } 
    catch (const std::exception& e) {
        Display::close();
        std::cerr << "Error: " << e.what() << std::endl;
        py::print("Python path:",sys.attr("path"));
        throw std::runtime_error("plotting failure2. See above for error message from python.");
        //TODO skip plotting if fail 
    }   
    #endif //PYBIND
    //py_clear_frame();
    //py_plot_frame(knot,density);        
}