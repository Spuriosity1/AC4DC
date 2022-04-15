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

#include <vector>
#include <iostream>
#ifndef GRIDSPACING_CXX_H
#define GRIDSPACING_CXX_H

struct GridSpacing {
    // GridSpacing(){
    //     mode= linear;
    // }
    // GridSpacing(char c){
    //     mode= c;
    // }
    const static char linear = 0;
    const static char quadratic = 1;
    const static char exponential = 2;
    const static char hybrid = 3;
    const static char powerlaw = 4;
    const static char unknown = 101;
    char mode = unknown;
    size_t num_low = 0; // Used for hybrid spec only
    double transition_e = 0; // Also used to estimate Coulomb logarithm. THis should be split into a separate class
    unsigned zero_degree_0 = 0; // Number of derivatives to set to zero at the origin
    unsigned zero_degree_inf = 0; // Number of derivatives to set to zero at infinity
    double min_coulomb_density=0; // Density below which to ignotre Coulomb interactions 
    // (NOTE: This should have its own class, along with coulomb estimation cutoff)
};

namespace {
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, GridSpacing gs) {
        switch (gs.mode)
        {
        case GridSpacing::linear:
            os << "linear";
            break;
        case GridSpacing::quadratic:
            os << "quadratic";
            break;
        case GridSpacing::exponential:
            os << "exponential";
            break;
        case GridSpacing::hybrid:
            os << "hybrid linear-linear grid, transition at M="<<gs.num_low<<", e="<<gs.transition_e;
            break;
        case GridSpacing::powerlaw:
            os << "power-law grid going through M="<<gs.num_low<<", e="<<gs.transition_e;
            break;
        default:
            os << "Unknown grid type";
            break;
        }
        return os;
    }

    [[maybe_unused]] std::istream& operator>>(std::istream& is, GridSpacing& gs) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No grid type provided, defaulting to linear..."<<std::endl;
            gs.mode = GridSpacing::linear;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'l':
            gs.mode = GridSpacing::linear;
            break;
        case 'q':
            gs.mode = GridSpacing::quadratic;
            break;
        case 'e':
            gs.mode = GridSpacing::exponential;
            break;
        case 'h':
            gs.mode = GridSpacing::hybrid;
            break;
        case 'p':
            gs.mode = GridSpacing::powerlaw;
            break;
        default:
            std::cerr<<"Unrecognised grid type \""<<tmp<<"\""<<std::endl;
            gs.mode = GridSpacing::linear;
            break;
        }
        return is;
    }
}

#endif