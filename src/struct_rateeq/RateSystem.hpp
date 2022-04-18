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

#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include "Constant.h"
#include "FreeDistribution.h"

typedef std::vector<double> bound_t;

// dy/dy = F(y)
// this is the 'y'

// Class responsible for storing the system state.
class state_type
{
public:
    std::vector<bound_t> atomP; // Probabilities of state for all atoms.
    // atomp[atom_idx][orbital_idx];
    Distribution F; // Energy distribution function
    double bound_charge;

    state_type();

    // Critical vector-space devices
    state_type& operator+=(const state_type &s);
    state_type& operator*=(const double x);
    // state_type operator+(const state_type& s2);

    // state_type operator*(double x);
    // convenience members
    state_type& operator=(const double x);
    // state_type& operator=(const state_type& s2);

    double norm() const;

    // Defines number and style of atomP
    // Resizes the container to fit all of the states present in the atom ensemble
    static void set_P_shape(const std::vector<RateData::Atom>& atomsys);
    static void set_P_shape(const std::vector<size_t>& shape) {
        P_sizes = shape;
    }
    static size_t P_size(size_t a) {
        return P_sizes[a];
    }
    static size_t num_atoms() {
        return P_sizes.size();
    }

private:
    static std::vector<size_t> P_sizes;
};

std::ostream& operator<<(std::ostream& os, const state_type& st);
std::ostream& operator<<(std::ostream& os, const bound_t& dist);
std::ostream& operator<<(std::ostream& os, const Distribution& dist);

#endif /* end of include guard: RATESYSTEM_CXX_H */
