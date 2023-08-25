/**
 * @file ElectronRateSolver.h
 * @authors Alaric Sanders & Spencer Passmore 
 * @brief Defines the ElectronRateSolver class, which executes the high-level operations of the electron ODE solving.
 * @details This is the central hub of Sanders' continuum plasma extension, and is called by main.cpp. This was historically named ElectronSolver.
 * @todo The cross-sectional computations should be computed before the class is called, in main.cpp, rather than within the class, decoupling the parts of the code.
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
// (C) Alaric Sanders 2020

#ifndef SYS_SOLVER_CXX_H
#define SYS_SOLVER_CXX_H

// #include <boost/numeric/odeint.hpp>
#include "HybridIntegrator.hpp"
#include "RateSystem.h"
#include "Constant.h"
#include "MolInp.h"
#include "Input.h"
#include "Pulse.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <chrono>

// using namespace boost::numeric::odeint;

// template<size_t Steps, typename State, typename Value = double,
//          typename Deriv = State, typename Time = Value,
//          typename Algebra = range_algebra,
//          typename Operations = default_operations,
//          typename Resizer = initially_resizer>

#define MAX_T_PTS 1e6

class ElectronRateSolver : protected ode::Hybrid<state_type>
{
public:
    ElectronRateSolver(const char* filename, ofstream& log);
    /// Solve the rate equations
    void solve(ofstream & _log, const string& tmp_data_folder);
    void save(const std::string& folder);
    /// Sets up the rate equations, which requires computing the atomic cross-sections/avg. transition rates to get the coefficients.
    void set_up_grid_and_compute_cross_sections(std::ofstream& _log, bool init,size_t step = 0,bool force_update = false); //bool recalc=true);
    /// creates the tensor of coefficients 
    void initialise_rates();
    void tokenise(std::string str, std::vector<double> &out, const size_t start_idx = 0, const char delim = ' ');

    /// Number of secs taken for simulation to run
    long secs;

    /// What's the time Mr Wolf?
    std::chrono::duration<double, std::milli> 
    display_time, plot_time, dyn_dt_time, backup_time, pre_ode_time, // pre_ode
    dyn_grid_time, user_input_time, post_ode_time,  // post_ode
    pre_tbr_time, transport_time, eii_time, tbr_time,  // sys_bound
    ee_time, apply_delta_time; //sys_ee 

    std::chrono::system_clock::time_point time_of_last_save;   
    #ifndef NO_BACKUP_SAVING
    std::chrono::minutes minutes_per_save{60};
    #else
    std::chrono::minutes minutes_per_save{99999999};
    #endif

protected:
    double IVP_step_tolerance = 5e-3;
    MolInp input_params;  // (Note this is initialised/constructed in the above constructor)  // TODO need to refactor to store variables that we change later rather than alter input_params directly. Currently doing a hybrid of this.
    ManualGridBoundaries elec_grid_regions;
    Cutoffs param_cutoffs; 
    double first_gp_min_E = 4/Constant::eV_per_Ha; // warning: going below ~4 eV leads to a much greater number of steps needed for little benefit. Though there's potential to increase dt once MB grid points go higher I suppose.
    Pulse pf;
    double timespan_au; // Atomic units
    double simulation_start_time;  // [Au]
    double simulation_resume_time; // [Au] same as simulation_start_time unless loading simulation state.
    double simulation_end_time;  // [Au]    
    double fraction_of_pulse_simulated;
    double grid_update_period; // time period between dynamic grid updates.

    void load_filtration_file(){};
    // Model parameters

    // arrays computed at class initialisation
    vector<vector<eiiGraph> > RATE_EII;
    vector<vector<eiiGraph> > RATE_TBR;

    // Dynamic grid
    void dirac_energy_bounds(size_t step, std::vector<double>& maximums, std::vector<double>& minimums, std::vector<double>& peaks, bool allow_shrinkage,size_t num_peaks, double peak_min_density);
    void mb_energy_bounds(size_t step, double& max, double& min, double& peak_density, bool allow_shrinkage);
    void transition_energy(size_t step, double& g_min);
    double approx_nearest_peak(size_t step, double start_energy,double del_energy, size_t min_sequential, double min = -1, double max =-1);  
    double approx_nearest_trough(size_t step, double start_energy,double del_energy, size_t min_sequential, double min = -1, double max =-1,bool accept_negatives = false);  
    double approx_regime_trough(size_t step, double lower_bound, double upper_bound_global,double upper_bound_local, double del_energy);
    double nearest_inflection(size_t step, double start_energy,double del_energy, size_t min_sequential, double min = -1, double max =-1);  
    double approx_regime_bound(size_t step, double start_energy,double del_energy, size_t min_sequential, double min_distance = 40, double min_inflection_fract = 1./4., double min = -1, double max=-1);  
    double approx_regime_peak(size_t step, double lower_bound, double upper_bound, double del_energy);  
    std::vector<double> approx_regime_peaks(size_t step, double lower_bound, double upper_bound, double del_energy, size_t num_peaks = 1, double min_density = 0);
    void precompute_gamma_coeffs(); // populates above two tensors
    void set_initial_conditions();

    // Dynamic time steps
    size_t load_checkpoint_and_decrease_dt(ofstream& _log, size_t current_n, Checkpoint _checkpoint);
    void increase_dt(ofstream& _log, size_t current_n); // bugged


    int steps_per_time_update;    

    /////// Overrides virtual system state methods
    /**
     * @brief         
     * @details 1. Displays general information
     * 2. Handles live plotting
     * 3. Handles dynamic updates to dt (and checkpoints, which also handle case of NaN being encountered by solver).
     */
    void pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update);
    /// general dynamics (uses explicit method)
    void sys_bound(const state_type& s, state_type& sdot, state_type& s_bg, const double t); 
    /// electron-electron (uses implicit method)
    void sys_ee(const state_type& s, state_type& sdot); 
    /**
     * @brief 
     * @details 1. Handles dynamic grid updates 
     * 2. Checks if user wants to exit simulation early, truncating the states y and times t, and returning 1 if so.
     * @return int
     */
    int post_ode_step(ofstream& _log, size_t& n);
    /////// 

    bool hasRates = false; // flags whether Store has been populated yet.
    void copyInput(const std::string& src,const std::string& dir);
    /// Log macros defined in config.h
    void log_config_settings(ofstream& _log);

    //// Saving/Loading
    /// Number of steps outputted, maximum determined by input_params.Out_T_size()
    double min_outputted_points = 25;
    int num_steps_out;
    /// We save last few steps so we can load simulations. This member determines how many of those sequential last steps are outputted
    const int extra_fine_steps_out = 10;
    /// handles deletion of given file
    void file_delete_check(const std::string& fname);
    /// Saves a table of free-electron dynamics to file fname
    void saveFree(const std::string& file);
    void saveFreeRaw(const std::string& fname);
    /// Loads the table of free-electron dynamics at the given time
    void loadFreeRaw_and_times();
    void saveKnots(const std::string& fname);
    void loadKnots();

    void loadBound();
    /// For each atom, saves a table of bound-electron dynamics to folder dir.
    void saveBound(const std::string& folder);
    static double convert_str_time(string str_time);// Helper function to convert string from input file to time.
    static double round_time(double time, const bool ceil = false, const bool floor = false);

    void set_grid_regions(ManualGridBoundaries gb);
    void set_starting_state();
    state_type get_ground_state();
    void update_grid(ofstream& _log, size_t latest_step, bool force_update = false);
    void reload_grid(ofstream& _log, size_t latest_step, std::vector<double> knots, std::vector<state_type> next_ode_states_used);

    //void high_energy_stability_check();
    string its_dinner_time(std::vector<std::chrono::duration<double, std::milli>> times, std::vector<std::string> tags);
    
    /// Folder that saves data periodically (period determined by minutes_per_save)
    string data_backup_folder; 
    bool grid_initialised = false;


    // Rates tracking: (this should be a class esp. due to checkpoint loading but I don't have time to do this properly)
    std::vector<double> bound_transport_rate;
    std::vector<double> photo_rate;
    std::vector<double> fluor_rate;
    std::vector<double> auger_rate;
    std::vector<double> eii_rate;
    std::vector<double> tbr_rate;
};


#endif /* end of include guard: SYS_SOLVER_CXX_H */
