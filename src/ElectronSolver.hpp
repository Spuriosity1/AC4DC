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
// (C) Alaric Sanders 2020

#ifndef SYS_SOLVER_CXX_H
#define SYS_SOLVER_CXX_H

/*
This file should arguably be called RateEquationSOlver, however, for historical reasons it is not.
*/
// #include <boost/numeric/odeint.hpp>
#include "HybridIntegrator.hpp"
#include "RateSystem.hpp"
#include "Constant.hpp"
#include "MolInp.hpp"
#include "Input.hpp"
#include "Pulse.hpp"
#include <iostream>
#include <stdexcept>

// using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

#define MAX_T_PTS 1e6

class ElectronSolver : private ode::Hybrid<state_type>
{
public:
    /**
     * @brief Construct a new Electron Solver "God Object"
     * 
     * @param filename the path to the input .mol file
     * @param log the path that logfile information will be stored in
     */
    ElectronSolver(const char* filename, std::ofstream& log) :
    Hybrid(3), input_params(filename, log), pf() // (order Adams method)
    {
        pf.set_shape(input_params.pulse_shape);
        pf.set_pulse(input_params.Fluence(), input_params.Width());
        timespan_au = input_params.Width()*4;
    }
    void solve();
    void save(const std::string& folder);
    void compute_cross_sections(std::ofstream& _log, bool recalc=true);
private:
    MolInp input_params;
    Pulse pf;
    double timespan_au; // Atomic units
    // Model parameters
    

    // arrays computed at class initialisation
    std::vector<std::vector<eiiGraph> > RATE_EII;
    std::vector<std::vector<eiiGraph> > RATE_TBR;

    void get_energy_bounds(double& max, double& min);
    void precompute_gamma_coeffs(); // populates above two tensors
    void set_initial_conditions();

    /////// Overrides virtual system state methods
    void sys(const state_type& s, state_type& sdot, const double t); // general dynamics (uses explicit mehtod)
    void sys2(const state_type& s, state_type& sdot, const double t); // electron-electron (uses implicit method)
    /////// 

    bool hasRates = false; // flags whether Store has been populated yet.
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    void saveBound(const std::string& folder);
    state_type get_ground_state();


    bool good_state = true;
    double timestep_reached = 0;
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
