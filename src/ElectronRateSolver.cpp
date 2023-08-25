/**
 * @file ElectronRateSolver.cpp
 * @author Alaric Sanders & Spencer Passmore
 * @brief @copybrief ElectronRateSolver.h
 * @details @copydetail ElectronRateSolver.h
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
// (C) Spencer Passmore 2023


#include "ElectronRateSolver.h"
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include "SplineIntegral.h"
#ifdef NCURSES
#include "Display.h"
#endif
#include <fstream>
#include <iostream>
#include <algorithm>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <math.h>
#include <omp.h>
#include "config.h"


// The Constructor


ElectronRateSolver::ElectronRateSolver(const char* filename, ofstream& log) :
    Hybrid(3), input_params(filename, log), pf() // (order Adams method), argument of Hybrid = order. num steps used for implicit adams-moulton = order - 1.
    {
        log_config_settings(log);

        param_cutoffs = input_params.param_cutoffs;
        
        pf.set_shape(input_params.pulse_shape);
        pf.set_pulse(input_params.Fluence(), input_params.Width());
        
        timespan_au = round_time(input_params.timespan_factor*input_params.Width());
        if (input_params.pulse_shape ==  PulseShape::square){  //-FWHM <= t <= 3FWHM (we've squished the pulse into the first FWHM.)
            if (timespan_au <= 0)  // 0 by default
                timespan_au = round_time(input_params.Width()*4,true); // 4*FWHM, capturing effectively entire pulse.
            simulation_start_time = -input_params.Width(); // e.g. if the timespan_au is 10 fwhm, the first fwhm is the pulse.
        } else { //-1.2FWHM <= t <= 1.2 FWHM
            if (timespan_au <= 0) // 0 by default
                timespan_au = round_time(input_params.Width()*2.4,true);   // > 99% of pulse, (similar to Neutze)
            
            simulation_start_time = -timespan_au/2;
            if (input_params.negative_timespan_factor != 0)
                simulation_start_time = -input_params.negative_timespan_factor*input_params.Width();
        }
        simulation_start_time = round_time(simulation_start_time);
        simulation_end_time =  simulation_start_time + timespan_au; 
        std::cout<<"\033[33m"<<"Imaging from "<<simulation_start_time/input_params.Width() <<" FWHM to " 
        << simulation_end_time/input_params.Width() <<" FWHM"<<"\033[0m"<<std::endl;
        
        if (input_params.Cutoff_Inputted()){
            if (input_params.Simulation_Cutoff() > simulation_end_time){
                std::cout <<"Ignoring cutoff, as cutoff is past the end time."<<endl;
            }
            simulation_end_time = min(input_params.Simulation_Cutoff(),simulation_end_time); 
        }
        simulation_resume_time = simulation_start_time;  // If loading simulation, is overridden by ElectronRateSolver::set_starting_state();
         

        set_grid_regions(input_params.elec_grid_regions);

        grid_update_period = input_params.Grid_Update_Period();  // TODO the grid update period should be made to be at least 3x (probably much more) longer with a gaussian pulse, since early times need it to be updated far less often to avoid instability for such a pulse. 
        steps_per_time_update = max(1 , (int)(input_params.time_update_gap/(timespan_au/input_params.num_time_steps))); 

        if(input_params.Filtration_File() != ""){
            load_filtration_file();
        }

        this->time_of_last_save = std::chrono::high_resolution_clock::now(); 
    }


void ElectronRateSolver::set_starting_state(){
    
    if (input_params.Load_Folder() != ""){ 
        cout << "[ Plasma ] loading sim state from specified files." << endl;
        loadFreeRaw_and_times();
        if (input_params.elec_grid_type.mode == GridSpacing::dynamic)
            loadKnots();
        loadBound();
        simulation_resume_time = t.back();
    }
    else{
        cout << "[ Plasma ] Creating ground state" << endl;
        this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, IVP_step_tolerance);
    }
}

state_type ElectronRateSolver::get_ground_state() {
    state_type initial_condition;
    assert(initial_condition.atomP.size() == input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        initial_condition.atomP[a][0] = input_params.Store[a].nAtoms;
        for(size_t i=1; i<initial_condition.atomP.size(); i++) {
            initial_condition.atomP[a][i] = 0.;
        }
    }
    initial_condition.F=0;
    // std::cout<<"[ Rate Solver ] initial condition:"<<std::endl;
    // std::cout<<"[ Rate Solver ] "<<initial_condition<<std::endl;
    // std::cout << endl;
    return initial_condition;
}
//state_type ElectronRateSolver::set_zero_y(){zero_y = 0*get_ground_state();}

// grid initilisation call order as of 3/04/23 TODO make clearer
// set_up_grid_and_compute_cross_sections(init = true) is called first
//  if static grid, use given static grid  
//  if dynamic grid, use default start
// Then call set_starting_state() (calls in set_up_grid_... now).
//  if loading a simulation:   
//      if static grid, transform loaded data to the given static grid
//      if dynamic grid, transform grid to the loaded data's grid/basis. 
//          Detect the transition energy
// If dynamic grid, call set_up_grid_and_compute_cross_sections(init = false) 
// after each non-zero integer multiple of the update period. 
void ElectronRateSolver::set_up_grid_and_compute_cross_sections(std::ofstream& _log, bool init,size_t step, bool force_update) {//, bool recalc) {
    if (!init && input_params.elec_grid_type.mode != GridSpacing::dynamic){
        return; 
    }
    // recalc being tuned off seems to be broken currently. We only do it once so I'm not bothering to deal with it. -S.P.
    // bool recalc = false;
    // #ifdef RECALC_HARTREE
    // recalc = true;
    // #endif
    bool recalc = true;
    if (init ){
        std::cout << "[ HF ] Computing bound rates for species' allowed orbital configurations" << std::endl;
        input_params.calc_rates(_log, recalc);
        hasRates = true;
    }
    
    

    /// Set up the grid
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        double old_trans_e = param_cutoffs.transition_e;
        if(init){
            // Empirically-based guesses for good starts.
            std::cout << "[ Dynamic Grid ] Starting with initial grid guess"<<std::endl;
            _log << "[ Dynamic Grid ] Starting with initial grid guess"<<endl;
            double e = 1/Constant::eV_per_Ha;
            param_cutoffs.transition_e = 250*e;
            regimes.mb_peak=0; regimes.mb_min=first_gp_min_E; regimes.mb_max=10*e;
            // cast a wide net to ensure we capture all dirac peaks. // TODO? there has to be a much better way than this.
            double photo_peak = -1, photo_min = 1e9, photo_max = -1e9; 
            for(auto& atom : input_params.Store) {
                for(auto& r : atom.Photo) {      
                    //photo_min = min(photo_min,r.energy*5/6);
                    photo_min = max(100*e,min(photo_min,r.energy*5/6));
                    photo_max = max(photo_max,r.energy*7/6);
                }
            }
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                regimes.dirac_minimums[i] = photo_min + i*(photo_max-photo_min)/regimes.num_dirac_peaks;
                regimes.dirac_maximums[i] = photo_min + (i+1)*(photo_max-photo_min)/regimes.num_dirac_peaks;
            }
            Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions, input_params.elec_grid_preset);
        }
        else{ //TODO? reset dirac bounds on first one of these calls (since they are initialised to be very wide) IF shrinkage is to be turned off for dirac regions.
            double old_mb_min = regimes.mb_min, old_mb_max = regimes.mb_max;
            auto old_dirac_mins = regimes.dirac_minimums, old_dirac_maxs = regimes.dirac_maximums;
            
            double last_trans_e = param_cutoffs.transition_e;
            double dirac_peak_cutoff_density = y[step].F(last_trans_e)*last_trans_e*2.5; 
            // TODO this depending on last transition is dangerous for large grid update periods, or low fluences as the peaks tend to spread out.
            if (last_trans_e <= 600/Constant::eV_per_Ha){
                // ad hoc fix - early on density at trans_e point in normalised density dist. is a bit higher than rest of sim since transition energy is stuck at default 250 eV 
                dirac_peak_cutoff_density = y[step].F(last_trans_e)*last_trans_e*1.25; 
            }
            dirac_energy_bounds(step,regimes.dirac_maximums,regimes.dirac_minimums,regimes.dirac_peaks,true,regimes.num_dirac_peaks,dirac_peak_cutoff_density);
            mb_energy_bounds(step,regimes.mb_max,regimes.mb_min,regimes.mb_peak,true); //TODO make allowing shrinkage an input option?
            // Check whether we should update
            if (regimes.mb_min == old_mb_min && regimes.mb_max == old_mb_max){
                if (force_update) _log <<"Forcing grid update"<< endl; // should be unlikely w/o adaptive b spline stuff? 
                else{
                    // We must be still at the start of the simulation, it's unnecessary, even detrimental, to update.
                    _log << "SKIPPED grid update as Maxwell-Boltzmann region's grid point range was the same"<<endl;
                    regimes.dirac_minimums = old_dirac_mins; 
                    regimes.dirac_maximums = old_dirac_maxs;
                    return;
                }
            }             
            transition_energy(step, param_cutoffs.transition_e);
            // [More general basis refinement, but disabled because scope.]           
            //if (y[step].F.seek_basis(_log, input_params.elec_grid_type, step, param_cutoffs)){
                // failed, use regular method.
            Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
                //_log << "B-spline failed to find convergent basis, using backup dynamic grid."; _log.flush();   Disabled
            //}
            //else{
            //     _log << "B-spline basis successfully converged"; _log.flush();
            //}
        } 
        if (old_trans_e != param_cutoffs.transition_e){
            std::cout <<"thermal cutoff energy updated from:\n"
            <<old_trans_e*Constant::eV_per_Ha
            << "\nto:\n"
            << param_cutoffs.transition_e*Constant::eV_per_Ha<< std::endl;
        }
        else std::cout<<"Thermal cutoff energy had NO update." <<std::endl;    

    }
    else{
        _log << "[ Manual Grid ] Setting static grid" <<endl;
        Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
    }


    
    // Set up the container class to have the correct size
    state_type::set_P_shape(input_params.Store);  // TODO should just do on init?


    if (init){
        // Set up the rate equations (setup called from parent Adams_B  M)
        set_starting_state(); // possibly redundant steps in here.
    }

    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ New Knots ] -------------------\n" 
            "Time: "<<t[step]*Constant::fs_per_au <<" fs; "<<"Step: "<<step<<"\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n"; 
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                _log<<"Photo [peak; range]: "<<regimes.dirac_peaks[i]*e<< "; " << regimes.dirac_minimums[i]*e<<" - "<<regimes.dirac_maximums[i]*e<<"\n";
            }
            _log <<"Transition energy: "<<param_cutoffs.transition_e*e<<"\n"  
            << "Grid size: "<<Distribution::size<<"\n"
            << "-----------------------------------------------------" 
            << endl;
        }        
    }

    // create the tensor of coefficients 
    initialise_rates();
}

void ElectronRateSolver::set_grid_regions(ManualGridBoundaries gb){
    elec_grid_regions = gb;
}

void ElectronRateSolver::solve(ofstream & _log, const std::string& tmp_data_folder) {
    assert (hasRates || "Rates weren't calculated!\n");
    auto start = std::chrono::system_clock::now();

    data_backup_folder = tmp_data_folder;

    if (simulation_resume_time > simulation_end_time){cout << "\033[91m[ Error ]\033[0m Simulation is attempting to resume from loaded data with a time after that which it is set to end at." << endl;}


    
    // generate text to display during simulation run
    std::stringstream plasma_header; 
    const string banner = "================================================================================";
    plasma_header<<banner<<"\n\r";
    if (input_params.time_update_gap > 0){
        plasma_header<< "Updating display every " <<input_params.time_update_gap*Constant::fs_per_au<<" fs."<<"\n\r";
    }
    if (simulation_resume_time != simulation_start_time){
        plasma_header/*<<"\033[33m"*/<<"Loaded simulation at:  "/*<<"\033[0m"*/<<(simulation_resume_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    }
    plasma_header/*<<"\033[33m"*/<<"Final time step:  "/*<<"\033[0m"*/<<(simulation_end_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic)
        plasma_header <<"[ sim ] Grid update period: "<<grid_update_period * Constant::fs_per_au<<" fs"<<"\n\r";
    else 
        plasma_header << "[ sim ] Using static grid" << "\n\r";

    plasma_header<<"[ Rate Solver ] Using initial timestep size of "<<this->dt*Constant::fs_per_au<<" fs"<<"\n\r";
    plasma_header<<banner<<"\n\r";

    steps_per_grid_transform =  round(input_params.Num_Time_Steps()*(grid_update_period/timespan_au));


    std::cout << plasma_header.str()<<std::flush; // display in regular terminal, so that it is still visible after end of program

    #ifdef NCURSES
    Display::header += plasma_header.str(); // display this in ncurses screen
    //Display::deactivate();
    #endif

    // Call hybrid integrator to iterate through the time steps (good state)
    good_state = true;
    size_t npoints = (simulation_end_time - simulation_start_time)/this->dt + 1;
    this->iterate(_log,simulation_start_time, npoints, simulation_resume_time, steps_per_time_update); // Inherited from ABM

    
    
    // 12-04 leaving this for now since I'm using it to catch really bad errors, but it now serves no more than that in practice.
    // (original: Uses finer time steps to attempt to resolve NaN encountered in ODE solving.) 
    double time = this->t[0];
    int retries = 0;  // int retries = 1;  <---- Turned off for now
    while (!good_state) {
        throw std::runtime_error("For an unexpected reason, a bad state was not or ceased being attempted to be resolved with finer timesteps.");
    }

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = std::chrono::system_clock::to_time_t(end);
    cout << "[ Solver ] finished computation at " << ctime(&end_time) << endl;
    secs = elapsed_seconds.count();
    
    // I'm sorry for this
    std::vector<std::chrono::duration<double, std::milli>> times{
    display_time, plot_time, dyn_dt_time, backup_time, pre_ode_time,
    dyn_grid_time, user_input_time, post_ode_time,
    pre_tbr_time, transport_time, eii_time, tbr_time,
    ee_time, apply_delta_time
    };
    std::vector<std::string> tags{
    "display", "live plotting", "dt updates", "data backups", "pre_ode()",
    "dynamic grid updates", "user input detection", "post_ode()",
    "misc bound processes", "bound-e transport", "get_Q_eii()", "get_Q_tbr()",
    "get_Q_ee()", "applyDeltaF()"
    };
    
    stringstream solver_times;
    solver_times <<"\n[ Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << "\n";
    solver_times << its_dinner_time(times,tags);
    std::cout<<banner<<"\033[1;32m"<<solver_times.str()<<"\033[0m\n"<<endl;
    _log <<banner<<"\n"<<solver_times.str(); _log.flush();
}

/// What the fuck?
string ElectronRateSolver::its_dinner_time(std::vector<std::chrono::duration<double, std::milli>> times, std::vector<std::string> tags){
    std::stringstream out;
    for(size_t i = 0; i < times.size(); i++){
        auto m = std::chrono::duration_cast<std::chrono::minutes>(times[i]); 
        auto s = std::chrono::duration_cast<std::chrono::seconds>(times[i]);
        out <<"[ Solver ] "<<tags[i]<<" took "<<m.count()<<"m "<<s.count()%60<< "s" << "\n";
    }
    return out.str();
}


// Populates RATE_EII and RATE_TBR based on calls to Distribution::Gamma_eii/ Distribution::Gamma_tbr
// 
void ElectronRateSolver::precompute_gamma_coeffs() {
    std::cout<<"[ Gamma precalc ] Beginning coefficient computation..."<<std::endl;
    size_t N = Distribution::size;
    for (size_t a = 0; a < input_params.Store.size(); a++) {
        std::cout<<"\n[ Gamma precalc ] Atom "<<a+1<<"/"<<input_params.Store.size()<<std::endl;
        auto eiiVec = input_params.Store[a].EIIparams;
        vector<RateData::InverseEIIdata> tbrVec = RateData::inverse(eiiVec);
        size_t counter=1;
        #pragma omp parallel default(none) shared(a, N, counter, tbrVec, eiiVec, RATE_EII, RATE_TBR, std::cout)
		{
			#pragma omp for schedule(dynamic) nowait
            for (size_t n=0; n<N; n++) {
                #pragma omp critical
                {
                    std::cout<<"\r[ Gamma precalc ] "<<counter<<"/"<<N<<" thread "<<omp_get_thread_num()<<std::flush;
                    counter++;
                }
                #ifndef NO_EII
                Distribution::Gamma_eii(RATE_EII[a][n], eiiVec, n);
                #endif
                // Weird indexing exploits symmetry of Gamma_TBR
                // such that only half of the coefficients are stored
                #ifndef NO_TBR
                for (size_t m=n+1; m<N; m++) {
                    size_t k = N + (N*(N-1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                    // // k = N... N(N+1)/2
                    Distribution::Gamma_tbr(RATE_TBR[a][k], tbrVec, n, m);
                }
                Distribution::Gamma_tbr(RATE_TBR[a][n], tbrVec, n, n);
                #endif
                
            }
        }
    }
    std::cout<<"\n[ Gamma precalc ] Done."<<std::endl;
}


// The Big One: The global ODE function, separated into sys_bound and sys_ee.
// Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j] = d/dt(average-atomic-state)
// d/dt f = Q_B[f](t)                      = d/dt(electron-energy-distribution)  // TODO what is Q_B? Q_{Bound}? Q_{ee} contributes though. 

// Non-stiff part of the system. Bound-electron dynamics.
//void ElectronRateSolver::sys_bound(const state_type& s, state_type& sdot, state_type& s_bg ,const double t) {
void ElectronRateSolver::sys_bound(const state_type& s, state_type& sdot, state_type& s_bg ,const double t) {
    const int threads = input_params.Plasma_Threads();
    
    sdot=0;
    Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);

    // Incorrect currently, Need to implement this with the ODE solver somehow :/
    photo_rate = fluor_rate = auger_rate = bound_transport_rate = tbr_rate = eii_rate = {0}; 
    // photo_rate.push_back(0); 
    // fluor_rate.push_back(0); 
    // auger_rate.push_back(0); 
    // bound_transport_rate.push_back(0);
    // tbr_rate.push_back(0); 
    // eii_rate.push_back(0); 
    if (!good_state) return;
    for (size_t a = 0; a < s.atomP.size(); a++) {
        auto t9 = std::chrono::high_resolution_clock::now();
        const bound_t& P = s.atomP[a];
        bound_t& Pdot = sdot.atomP[a];
        
        #ifdef DEBUG_BOUND
        for(size_t i=0;i < Pdot.size();i++){
            assert(Pdot[i] + P[i] >= 0);
        }
        #endif

        // PHOTOIONISATION
        double J = pf(t); // photon flux in atomic units
        for ( auto& r : input_params.Store[a].Photo) {
            double tmp = r.val*J*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(r.energy, r.val*J*P[r.from]);
            sdot.bound_charge +=  tmp;
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*J*P[r.from]);
        }
        photo_rate.back() += sdot.bound_charge;
        double old_bound_charge = sdot.bound_charge;
        #ifdef DEBUG_BOUND
        for(size_t i=0;i < Pdot.size();i++){
            assert(Pdot[i] + P[i] >= 0);
        }
        #endif    
        // FLUORESCENCE
        for ( auto& r : input_params.Store[a].Fluor) {
            double tmp = r.val*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            fluor_rate.back() +=tmp;
            // Energy from optical photon assumed lost
        }
        #ifdef DEBUG_BOUND
        for(size_t i=0;i < Pdot.size();i++){
            assert(Pdot[i] + P[i] >= 0);
        }
        #endif    
        // AUGER
        for ( auto& r : input_params.Store[a].Auger) {
            double tmp = r.val*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(r.energy, r.val*P[r.from]);
            // sdot.F.add_maxwellian(r.energy*2./3., r.val*P[r.from]);
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*P[r.from]);
            sdot.bound_charge +=  tmp;
        }
        auger_rate.back() += sdot.bound_charge - old_bound_charge;
        old_bound_charge = sdot.bound_charge;
        /// CHARGE TRANSPORT
        auto t_transport = std::chrono::high_resolution_clock::now();
        // Every heavy atom state is moved (or not moved) in ratios corresponding to the population of light atom states, and divided by the user-defined charge transport timescale.
        // average charge transferred to heavy atoms is tracked. this charge, divided by the user-defined number of light atom neighbours, is a first order approximation for the local drainage of electrons. 
        // If the local drainage of electrons is 1, then we modify the ratios by assigning the population of a light atom state to the light atom state with 1 less valence electron. If it's 1.5, then we do the same, and then move half of those down to the state corresponding to 2 less valence electrons. (we just go to next orbital if we run out of valence electrons)
        /// Upper bound quick hack version, taking first index, transferring to next 3 indices. Should be relatively accurate for low fluences like in Galli 2015.
        #ifdef BOUND_GD_HACK
        size_t heavy_idx = 0;  // Gd
        std::vector<double> light_fractions = {0,0,0.5,0.5,0}; // Fraction of complex light atoms make up, in index order ({Gd,C,N, O, S} in this case)
        if (a == heavy_idx){
            for (auto& r : input_params.Store[a].EnergyConfig){
                if (r.receiver_index < 0) continue;

                for (size_t a2 = 0; a2 < s.atomP.size(); a2++) {
                    if (light_fractions[a2] == 0 || a2 == a) 
                        continue;
                    double transfer_factor = light_fractions[a2]*P[r.index]/y[0].atomP[a2][0]/dt/10; // VERY lazy OOM approximation corresponding to a very conservatively overestimated lifetime of 0.24 fs.
                    // Iterate through each donor state
                    // valence energy is *negative*
                    for (auto& r2 : input_params.Store[a2].EnergyConfig){
                        if (r.valence_energy >= r2.valence_energy || r2.donator_index < 0) 
                            continue;
                        double tmp = s.atomP[a2][r2.index]*transfer_factor;
                        // Donor state loses an electron
                        sdot.atomP[a2][r2.index] -= tmp; 
                        sdot.atomP[a2][r2.donator_index] += tmp;
                        // Receiver state gains an electron
                        Pdot[r.receiver_index] += tmp;
                        Pdot[r.index] -= tmp;             
                        // track rate
                        bound_transport_rate.back()+=tmp;                   
                    }
                }            
            }
        }
        #endif  //BOUND_GD_HACK
        transport_time += std::chrono::high_resolution_clock::now() - t_transport;        

        #ifdef DEBUG_BOUND
        for(size_t i=0;i < Pdot.size();i++){
            assert(Pdot[i] + P[i] >= 0);
        }
        #endif        

        // EII / TBR bound-state dynamics
        double Pdot_subst [Pdot.size()] = {0};    // subst = substitute.
        double sdot_bound_charge_eii_subst = 0; 
        double sdot_bound_charge_tbr_subst = 0; 
        size_t N = Distribution::size; 
        #pragma omp parallel for num_threads(threads) reduction(+ : Pdot_subst,sdot_bound_charge_eii_subst,sdot_bound_charge_tbr_subst)     
        for (size_t n=0; n<N; n++) {
            double tmp=0; // aggregator
            
            #ifndef NO_EII
            for (size_t init=0;  init<RATE_EII[a][n].size(); init++) {
                for (auto& finPair : RATE_EII[a][n][init]) {
                    tmp = finPair.val*s.F[n]*P[init];
                    Pdot_subst[finPair.idx] += tmp;
                    Pdot_subst[init] -= tmp;
                    sdot_bound_charge_eii_subst += tmp;   //TODO Not a minus sign like others, double check that's intentional. -S.P.
                }
            }
            #endif
            
            
            #ifndef NO_TBR
            // exploit the symmetry: strange indexing engineered to only store the upper triangular part.
            // Note that RATE_TBR has the same geometry as InverseEIIdata.
            for (size_t m=n+1; m<N; m++) {
                size_t k = N + (N*(N-1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                // k = N... N(N+1)/2-1
                // W += RATE_TBR[a][k]*s.F[n]*s.F[m]*2;
                for (size_t init=0;  init<RATE_TBR[a][k].size(); init++) {
                    for (auto& finPair : RATE_TBR[a][k][init]) {
                        tmp = finPair.val*s.F[n]*s.F[m]*P[init]*2;
                        Pdot_subst[finPair.idx] += tmp;
                        Pdot_subst[init] -= tmp;
                        sdot_bound_charge_tbr_subst -= tmp;
                    }
                }
            }
            // the diagonal
            // W += RATE_TBR[a][n]*s.F[n]*s.F[n];
            for (size_t init=0;  init<RATE_TBR[a][n].size(); init++) {
                for (auto& finPair : RATE_TBR[a][n][init]) {
                    tmp = finPair.val*s.F[n]*s.F[n]*P[init];
                    Pdot_subst[finPair.idx] += tmp;
                    Pdot_subst[init] -= tmp;
                    sdot_bound_charge_tbr_subst -= tmp;
                }
            }
            #endif
        }
        // Add parallel containers to their parent containers.
        for(size_t i=0;i < Pdot.size();i++){
            Pdot[i] += Pdot_subst[i];
            #ifdef DEBUG_BOUND
            assert(Pdot[i] + P[i] >= 0);
            // if(Pdot[i] + P[i] < 0){
            //     Pdot[i]= -P[i]*1.00001;
            // }
            #endif
        }
        sdot.bound_charge += sdot_bound_charge_eii_subst + sdot_bound_charge_tbr_subst;
        eii_rate.back() += sdot_bound_charge_eii_subst;
        tbr_rate.back() += sdot_bound_charge_tbr_subst;

        auto t10 = std::chrono::high_resolution_clock::now();
        pre_tbr_time += t10 - t9;

        // Free-electron parts
        #ifdef NO_EII
        #warning No impact ionisation
        #else
        auto t1 = std::chrono::high_resolution_clock::now();
        s.F.get_Q_eii(vec_dqdt, a, P, threads);
        auto t2 = std::chrono::high_resolution_clock::now();
        eii_time += t2 - t1;
        #endif
        #ifdef NO_TBR
        #warning No three-body recombination
        #else
        auto t3 = std::chrono::high_resolution_clock::now();
        s.F.get_Q_tbr(vec_dqdt, a, P, threads);  // Serially, this is the computational bulk of the program - S.P.
        auto t4 = std::chrono::high_resolution_clock::now();
        tbr_time += t4 - t3;
        #endif
    }

    // Add change to distribution
    sdot.F.applyDeltaF(vec_dqdt,threads);
    // if(input_params.Filtration_File() == "")
    //     // No background, use geometric-dependent loss.
    //     sdot.F.addLoss(s.F, input_params.loss_geometry, s.bound_charge);
    // else
    //     // background handles loss, apply photoelectron filtration with background
    //     sdot.F.addFiltration(s.F, s_bg.F,input_params.loss_geometry);
    sdot.F.addLoss(s.F, input_params.loss_geometry, s.bound_charge);
    
    if (isnan(s.norm()) || isnan(sdot.norm())) {
        cerr<<"NaN encountered in ODE iteration."<<endl;
        cerr<< "t = "<<t*Constant::fs_per_au<<"fs"<<endl;
        good_state = false;
        timestep_reached = t*Constant::fs_per_au;
    }    
}


// 'badly-behaved' i.e. stiff part of the system. Electron-electron interactions.
void ElectronRateSolver::sys_ee(const state_type& s, state_type& sdot) {
    sdot=0;
    Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);
    const int threads = input_params.Plasma_Threads(); 
    // compute the dfdt vector
    #ifdef NO_EE
    #warning No electron-electron interactions
    #else
    // Electron-electon repulsions
    auto t5 = std::chrono::high_resolution_clock::now();
    s.F.get_Q_ee(vec_dqdt, threads); 
    auto t6 = std::chrono::high_resolution_clock::now();
    ee_time += t6 - t5;
    #endif
    // 
    // Add change to distribution
    auto t7 = std::chrono::high_resolution_clock::now();
    sdot.F.applyDeltaF(vec_dqdt,threads);
    auto t8 = std::chrono::high_resolution_clock::now();
    apply_delta_time += t8 - t7;
}




//// Mid-simulation Non-ode functions

size_t ElectronRateSolver::load_checkpoint_and_decrease_dt(ofstream &_log, size_t current_n, Checkpoint _checkpoint){    
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    double fact = 1.5; // factor to increase time step density by (past the checkpoint). TODO need to implement max time steps.

    size_t n = _checkpoint.n; // index of step to load at 
    std::vector<double> saved_knots = _checkpoint.knots;
    std::vector<state_type> saved_states = _checkpoint.last_few_states;
    std::vector<double> saved_times = _checkpoint.last_few_times;

    this->regimes = _checkpoint.regimes; // needed for when we next update the regimes (it needs the previous regime for reference)

    std::vector<double> old_knots = Distribution::get_knot_energies();

    const string banner = "================================================================================";
    _log <<banner<<"\nEuler iterations exceeded beyond tolerable error at t="<<t[current_n]*Constant::fs_per_au
    <<". Decreasing remaining time step size from "<<this->dt*Constant::fs_per_au<<" fs to "<< dt/fact*Constant::fs_per_au<<" fs."<<endl;

    // REDUCE time step size by factor, then add on extra steps.
    this->dt/=fact;
    int remaining_steps = t.size() - (n+1);
    assert(remaining_steps > 0 && fact > 1);
    input_params.num_time_steps = ceil(t.size() + (fact - 1)*remaining_steps); // todo separate num steps from input params
    steps_per_time_update = max(1 , (int)(0.5+input_params.time_update_gap/(timespan_au/input_params.num_time_steps))); 
    t.resize(n+1); t.resize(input_params.num_time_steps);

    //TODO implement: size_t npoints = (t_final - t_initial)/this->dt + 1; input_params.num_time_steps is for the full number, not what we want.
    // Set up the t grid       
    for (size_t i=n+1; i<input_params.num_time_steps; i++){   // TODO check potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
        this->t[i] = this->t[i-1] + this->dt;
    }
    steps_per_grid_transform = round(steps_per_grid_transform*fact);

    _log <<"Reloading grid"<<endl;
    // with the states loaded the spline factors are as they were, we just need to load the appropriate knots.
    assert(saved_knots == Distribution::get_knots_from_history(n));
    for(auto elem : saved_states){
        assert(elem.F.container_size() == saved_states.back().F.container_size());
    }
    // Linearly interpolate the previous states as a first order correction. (We should use 4th order Runge-Kutta, but I don't have time to get that working with sys_ee incorporated)
    // last element of saved_states is y[n], so we just iterate up to y[n-1]
    std::vector<state_type> interpolated_states(saved_states.size());
    for(size_t i = 0; i < saved_states.size();i++){
        interpolated_states[i]=saved_states[i];
    }
    for(size_t i = 0; i < saved_states.size() - 1; i++){
        assert(saved_times.back()==t[n]);
        double intra_time = t[n] - (saved_states.size()-1 - i)*dt; 
        size_t loop_step = n;
        while(1){
            int rel_idx = saved_states.size() - (n - loop_step) - 1;
            if(saved_times[rel_idx] == intra_time){
                interpolated_states[i] = saved_states[rel_idx];
                break;
            }
            if(saved_times[rel_idx] < intra_time){
                // if t' = t(n) + u 
                // y(t') = y(n) + [y(n+1) - y(n)] * u/[t(n+1)-t(n)]  = y(n) + dy
                // find dy
                assert(rel_idx < saved_times.size());
                interpolated_states.at(i) = saved_states.at(rel_idx);
                interpolated_states[i] *= -1;
                interpolated_states[i] += saved_states.at(rel_idx +1);
                interpolated_states[i] *= (intra_time-saved_times[rel_idx])/(saved_times[rel_idx+1]-saved_times[rel_idx]);
                // add y(n) to get y(t')
                interpolated_states[i] += saved_states[rel_idx];
                break;
            }
            loop_step--;
            if (saved_states.size() < (n-loop_step)-1)
                throw runtime_error("Unexpected error, could not interpolate times when loading checkpoint.");
        }
    }
    photo_rate.resize(n); 
    fluor_rate.resize(n); 
    auger_rate.resize(n); 
    tbr_rate.resize(n); 
    eii_rate.resize(n); 
    reload_grid(_log, n, saved_knots,interpolated_states);
    

    /* [These ideas need a bit of effort and consideration if it is to be implemented, or even better a more effective adaptive step size algo could be implemented. But questionable if worth it.]
        // It is standard to use a method like fourth-order Runge-Kutta after changing the step size to get us going:
        // (TODO RK4 needs to have free electron dynamics added though).
        // We get indistinguishable results skipping this entirely (for a factor of 1.5 at least).
        
        /#
        n += this->order;
        for (size_t m = n-this->order; m < n; m++) {
                this->step_rk4(m);
        }
        std::cout.clear();
        #/ 
    

        // Remember when time step was decreased, and flag to increase time step size again
        switch (input_params.pulse_shape){

        case PulseShape::gaussian:
            // TODO impact ionisation dominates so we would need a better method to do this. Arguably it's not worth it. 
            // if (t[n] <= 0){
            //     times_to_increase_dt.push_back(-t[n]+0.05/Constant::fs_per_au);
            // }
            break;
        case PulseShape::square:
            break;  
        }   

    */

    _log << "Loaded checkpoint - resuming."<<endl;
    
    return n;
}

/**
 * @brief Increases time step size so as to reduce the number of remaining steps.
 * @param _log 
 * @param current_n 
 * @note current implementation means if a checkpoint was loaded multiple times (as can happen when we have multiple dt decreases within 2*checkpoint period,
 * then this will be called at the same time multiple times. Not critical to make it better though.  
 * An alternative would also be to make dt proportional to the incoming intensity. What be even better would be making it proportional to the
 * net magnitude of delta(density) / delta(t)
 */
void ElectronRateSolver::increase_dt(ofstream &_log, size_t current_n){
    // /* breaks things TODO may no longer break things since I found what was wrong... retry 0.5 fs FWHM mol file after fixing excessive computation time
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    
    size_t n = current_n;

    int remaining_steps = t.size() - (n+1);
    double fact = 1.5; // factor to increase time step size dt by.
    
    _log <<"Step size increase flagged at t="<<t[current_n]*Constant::fs_per_au
    <<". Increasing remaining time step size from "<<this->dt*Constant::fs_per_au<<" fs to "<< dt*fact*Constant::fs_per_au<<" fs."<<endl;
    
    // INCREASE time step size by factor, and reduce number of remaining steps.
    this->dt*=fact;
    assert(remaining_steps > 0 && fact > 1);
    input_params.num_time_steps = input_params.num_time_steps - (1-1/fact)*remaining_steps; // todo separate from input params
    t.resize(n+1); t.resize(input_params.num_time_steps);
    y.resize(n+1); y.resize(input_params.num_time_steps);
    steps_per_grid_transform = round(steps_per_grid_transform/fact);

    //TODO implement: size_t npoints = (t_final - t_initial)/this->dt + 1; input_params.num_time_steps is for the full number, not what we want.
    // Set up the t grid       
    for (int i=n+1; i<input_params.num_time_steps; i++){   // note potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
        this->t[i] = this->t[i-1] + this->dt;
    }
    // */
}


void ElectronRateSolver::update_grid(ofstream& _log, size_t latest_step, bool force_update){
    size_t n = latest_step;
    //// Set up new grid ////        
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    // Latest step is n, so we decide our new grid based on that step, then transform N = "order" of the prior points to the new basis.
    set_up_grid_and_compute_cross_sections(_log,false,n,force_update); // virtual function overridden by ElectronRateSolver
    std::cout.clear();
    //// We need some previous points needed to perform next ode step, so transform them to new basis ////
    std::vector<double> new_energies = Distribution::get_knot_energies();
    // zero_y is used as the starting state for each step, so we need to reset it so it has the right knots.
    this->zero_y = get_ground_state();
    this->zero_y *= 0.; // set it to Z E R O
                
    for (size_t m = n+1 - this->order; m < n+1; m++) {  // TODO turn off if not new knots??
        // We transform this step to the correct basis, but we also need a few steps to get us going, 
        // so we transform a few previous steps 
        // Kinda goofy but it's necessary due to the static variables. 
        assert(this-> order*2 < steps_per_grid_transform); // *2 factor is to guarantee loading works.
        Distribution::load_knots_from_history(n-1); // the n - 1 is correct. transform_basis takes us to basis at step n.
        y[m].F.transform_basis(new_energies);

        #ifdef DEBUG_BOUND
        for(size_t a = 0; a < y[m].atomP.size();a++)
            for(size_t i=0;i < y[m].atomP[a].size();i++){
                assert(y[m].atomP[a][i] >= 0);
            }        
        #endif
    }   
    //TODO rk4 (with sys_ee added) here? 
    initialise_transient_y((int)latest_step);
    // The next containers are made to have the correct size, as the initial state is set to tmp=zero_y and sdot is set to an empty state. 
} 

// Reloads grid using knots and the states needed for the first step of the ode 
// next_ode_states_used - includes current step
// Does not resize containers
void ElectronRateSolver::reload_grid(ofstream& _log, size_t latest_step, std::vector<double> knots, std::vector<state_type> next_ode_states_used){
    size_t n = latest_step;
   
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        // set basis to one given by knots
        while(Distribution::knots_history.back().step >= n && Distribution::knots_history.size() > 0){
            Distribution::knots_history.pop_back();
        }            
        Distribution::set_basis(latest_step, param_cutoffs, regimes,  knots);
    }

    // Set up the container class to have the correct size
    state_type::set_P_shape(input_params.Store);


    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ Reloaded Knots ] -------------------\n" 
            "Time: "<<t[latest_step]*Constant::fs_per_au <<" fs; "<<"Step: "<<latest_step<<"\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n"; 
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                _log<<"Photo [peak; range]: "<<regimes.dirac_peaks[i]*e<< "; " << regimes.dirac_minimums[i]*e<<" - "<<regimes.dirac_maximums[i]*e<<"\n";
            }
            _log <<"Transition energy: "<<param_cutoffs.transition_e*e<<"\n"  
            << "Grid size: "<<Distribution::size<<"\n"
            << "-----------------------------------------------------" 
            << endl;
        }        
    }

    
    this->zero_y = get_ground_state();
    this->zero_y *= 0.; // set it to Z E R O
    // Load distributions and also the previous few states that we need for the first ode step. (Unnecessary now since removed grid updates on checkpoint load, but keeping it here in case we go back to that.)
    assert(next_ode_states_used.size() == order);  // Order previous states + present state
    y.resize(n+1-next_ode_states_used.size());
    for(state_type state : next_ode_states_used){
        y.push_back(state);
    }
    y.resize(input_params.num_time_steps);  
    initialise_transient_y((int)latest_step);
    
    
    // update the tensor of coefficients 
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic)
        initialise_rates();    

    std::cout.clear();

    // The next containers are made to have the correct size, as the initial state is set to tmp=zero_y and sdot is set to an empty state. 
} 


void ElectronRateSolver::initialise_rates(){
    //TODO these aren't constant now - change to lowercase.
    RATE_EII.clear();
    RATE_TBR.clear();    
    RATE_EII.resize(input_params.Store.size());
    RATE_TBR.resize(input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        size_t N = Distribution::num_basis_funcs();
        RATE_EII[a].clear();
        RATE_TBR[a].clear();
        RATE_EII[a].resize(N);
        RATE_TBR[a].resize(N*(N+1)/2);
    }
    precompute_gamma_coeffs();
    Distribution::precompute_Q_coeffs(input_params.Store);    
}

// trivial instantiations
void ElectronRateSolver::pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update){};
int ElectronRateSolver::post_ode_step(ofstream& _log, size_t& n){};




