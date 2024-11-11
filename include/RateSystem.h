/**
 * @file RateSystem.h
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

#ifndef RATESYSTEM_CXX_H
#define RATESYSTEM_CXX_H


#include <sstream>
#include <assert.h>
#include <iostream>
#include "Constant.h"
#include "FreeDistribution.h"



/**
 * @brief Class responsible for storing the system state at a moment in time.
 * @details y[i] corresponds to t[i] throughout the plasma code, which is somewhat unfortunate.
 * y[i].F and y.atomP are everything you'd want to pop on a plot's axes for time t[i]... perhaps, even, on the y-axis against t? -S.P.
 */
class state_type
{
public:
    /// Probabilities of state for all atoms.
    std::vector<bound_t> atomP; 
    // Tracks sum total of photoionisation for all atoms.
    std::vector<double> cumulative_photo;     
    /// Energy distribution function
    Distribution F;   
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

    double norm(size_t _c) const;

    // Defines number and style of atomP
    // Resizes the container to fit all of the states present in the atom ensemble
    static void set_P_shape(const vector<RateData::Atom>& atomsys);
    static void set_P_shape(const vector<size_t>& shape) {
        P_sizes = shape;
    }
    static size_t P_size(size_t a) {
        return P_sizes[a];
    }
    static size_t num_atoms() {
        return P_sizes.size();
    }

private:
    static vector<size_t> P_sizes;
};

ostream& operator<<(ostream& os, const state_type& st);
ostream& operator<<(ostream& os, const bound_t& dist);
ostream& operator<<(ostream& os, const Distribution& dist);

// All f integrals have the form
// df(e)/dt = Q [f] (e)
// Expand in some basis f = \sum_k a_k(t) f_k
// \sum_k da_k(t)/dt f_k(e) = Q[f] (e)
// This does not make life easier, unless we take the f_k to have compact support.


/*         TODO: Fix these (do not currently habdle multiple P arrays implemented above.)
// Algebra definition for error-controlled steppers
state_type operator/( const state_type &s1 , const state_type &s2 ) {
    std::vector<double> tmpP = s1.P;
    std::vector<double> tmpf = s1.f;

    for (size_t i = 0; i < tmpP.size(); i++) {
        tmpP[i] /= s2.P[i];
    }
    for (size_t i = 0; i < tmpf.size(); i++) {
        tmpf[i] /= s2.f[i];
    }
    return state_type( tmpP, tmpf );
}

state_type abs( const state_type &s) {
    state_type tmp(s.P.size(), s.f.size());

    for (size_t i = 0; i < tmp.P.size(); i++) {
        tmp.P[i] = std::abs(tmp.P[i]);
    }
    for (size_t i = 0; i < tmp.f.size(); i++) {
        tmp.f[i] = std::abs(tmp.f[i]);
    }
    return state_type( tmp );
}

// also only for steppers with error control
namespace boost { namespace numeric { namespace odeint {
template<>
struct vector_space_norm_inf< state_type >
{
    typedef double result_type;
    double operator()( const state_type &p ) const
    {
        using std::abs;
        double max=0.;
        double tmp;
        for (size_t i = 0; i < p.P.size(); i++) {
            tmp = abs(p.P[i]);
            max = (max > tmp) ? max : tmp;
        }
        for (size_t i = 0; i < p.f.size(); i++) {
            tmp = abs(p.f[i]);
            max = (max > tmp) ? max : tmp;
        }
        return max;
    }
};
*/

//
// // Flag no resizing
// namespace boost { namespace numeric { namespace odeint {
// template<>
// struct is_resizeable< state_type >
// {
//     typedef boost::false_type type;
//     const static bool value = type::value;
// };
//
// } } }

// End of definitions

#endif /* end of include guard: RATESYSTEM_CXX_H */
