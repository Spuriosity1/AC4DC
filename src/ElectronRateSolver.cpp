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

        this->time_of_last_save = std::chrono::high_resolution_clock::now(); // This is currently occurring before calculation of rates but this should be fixed once the atomic code is separated from this code. 
    }


void ElectronRateSolver::set_starting_state(){
    // Set up the container class to have the correct size
    state_type::set_P_shape(input_params.Store);    
    if (input_params.Load_Folder() != ""){ 
        load_simulation_state();
    }
    else{         
        cout << "[ Plasma ] Creating ground state" << endl;
        double _dt = this->timespan_au/input_params.Num_Time_Steps();
        this->setup(get_initial_state(), _dt, IVP_step_tolerance);
        this->t[0] = simulation_start_time; // This is a fix so that log shows correct time for initial grid.
    }
    num_steps = round((simulation_end_time - simulation_start_time)/this->dt + 1);
}

state_type ElectronRateSolver::get_initial_state() {
    state_type initial_condition;
    assert(initial_condition.atomP.size() == input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        initial_condition.atomP[a][0] = input_params.Store[a].nAtoms;
        for(size_t i=1; i<initial_condition.atomP.size(); i++) {
            initial_condition.atomP[a][i] = 0.;
        }
    }
    initial_condition.F=0;
    initial_condition.cumulative_photo = std::vector(input_params.Store.size(),0.0);
    // std::cout<<"[ Rate Solver ] initial condition:"<<std::endl;
    // std::cout<<"[ Rate Solver ] "<<initial_condition<<std::endl;
    // std::cout << endl;
    return initial_condition;
}
//state_type ElectronRateSolver::set_zero_y(){zero_y = 0*get_initial_state();}

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
    bool recalc = false;
    if (init ){
        std::cout << "[ HF ] Peforming Hartree-Fock calculations for species' allowed orbital configurations" << std::endl;
        input_params.calc_rates(_log, recalc);
        hasRates = true;
        Distribution::num_continuums = 1;
        #ifndef TRACK_SINGLE_CONTINUUM
        Distribution::num_continuums = 1 + input_params.Store.size(); // Total + 1 for each atom
        #endif
    }
    
    
    
    /// Set up the grid
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        double old_trans_e = param_cutoffs.transition_e;
        if(init){
                if (input_params.Load_Folder()== ""){
                // Empirically-based guesses for good starts.
                std::cout << "[ Dynamic Grid ] Starting with initial grid guess"<<std::endl;
                _log << "[ Dynamic Grid ] Starting with initial grid guess"<<endl;
                }
                else{
                   _log << "[ Dynamic Grid ] Initialising dummy grid"<<endl;
                }
                double e = 1/Constant::eV_per_Ha;
                param_cutoffs.transition_e = 250*e;
                regimes.mb_peak=0; regimes.mb_min=Distribution::get_lowest_allowed_knot(); regimes.mb_max=10*e;
                // cast a wide net to ensure we capture all dirac peaks. // TODO? there has to be a much better way than this.
                double photo_min = 1e9, photo_max = -1e9; 
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
                Distribution::set_basis(step, input_params.elec_grid_type, param_cutoffs, regimes, elec_grid_regions);
        }
        else{ //TODO? reset dirac bounds on first one of these calls (since they are initialised to be very wide) IF shrinkage is to be turned off for dirac regions.
            double old_mb_min = regimes.mb_min, old_mb_max = regimes.mb_max;
            auto old_dirac_mins = regimes.dirac_minimums, old_dirac_maxs = regimes.dirac_maximums;
            
            double last_trans_e = param_cutoffs.transition_e;
            
            #ifndef NO_MIN_DIRAC_DENSITY
            double dirac_peak_cutoff_density = y[step].F(0,last_trans_e)*last_trans_e*2.5;            
            // TODO this depending on last transition is dangerous for large grid update periods, or low fluences as the peaks tend to spread out.
            if (last_trans_e <= 600/Constant::eV_per_Ha){
                // ad hoc fix - early on density at trans_e point in normalised density dist. is a bit higher than rest of sim since transition energy is stuck at default 250 eV 
                dirac_peak_cutoff_density = y[step].F(0,last_trans_e)*last_trans_e*1.25; 
            }
            # else 
            double dirac_peak_cutoff_density = 0;
            #endif // NO_MIN_DIRAC_DENSITY
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


    if (init){         
        // Set up the rate equations (setup called from parent Adams_B  M)
        set_starting_state(); // possibly redundant steps in here.
    }

    if (input_params.elec_grid_type.mode == GridSpacing::dynamic && (init==false || input_params.Load_Folder()== "")){
        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ New Knots ] -------------------\n" 
            "Time: "<<t[step]*Constant::fs_per_au <<" fs; "<<"Step: "<<step<<"\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n"; 
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                _log<<"Photo [peak; range]: "<<regimes.dirac_peaks[i]*e<< "; " << regimes.dirac_minimums[i]*e<<" - "<<regimes.dirac_maximums[i]*e<<"\n";
            }
            _log <<"Transition energy: "<<param_cutoffs.transition_e*e<<" eV\n"  
            << "Grid size: "<<Distribution::size<<"\n"
            << "-----------------------------------------------------" 
            << endl;
        }        
    }

    // create/update the tensor of coefficients corresponding to free distribution grid knots 
    compute_free_grid_rates();
    // zero_y is used as the starting state for each step, so it must match the current knot basis.
    set_zero_y();    
}

void ElectronRateSolver::set_grid_regions(ManualGridBoundaries gb){
    elec_grid_regions = gb;
}

void ElectronRateSolver::execute_solver(ofstream & _log, const std::string& tmp_data_folder) {
    assert (hasRates || "Rates weren't calculated!\n");
    auto start = std::chrono::system_clock::now();

    data_backup_folder = tmp_data_folder;

    if (simulation_resume_time > simulation_end_time){cout << "\033[91m[ Error ]\033[0m Simulation is attempting to resume from loaded data with a time after that which it is set to end at." << endl;}


    
    // generate text to display during simulation run
    std::stringstream plasma_header; 
    const string banner = "================================================================================";
    plasma_header<<banner<<"\n\r";
    if (input_params.time_update_gap > 0){
        plasma_header<< "Updating this display every " <<input_params.time_update_gap*Constant::fs_per_au<<" fs."<<"\n\r";
    }
    if (simulation_resume_time != simulation_start_time){
        plasma_header/*<<"\033[33m"*/<<"Loaded simulation at:  "/*<<"\033[0m"*/<<(simulation_resume_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    }
    plasma_header/*<<"\033[33m"*/<<"Final time step:  "/*<<"\033[0m"*/<<(simulation_end_time)*Constant::fs_per_au<<" fs"<<"\n\r";
    
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        plasma_header <<"[ Grid ] Preset: "<<input_params.elec_grid_preset.name<<"\n\r";
        plasma_header <<"[ Grid ] Update period: "<<grid_update_period * Constant::fs_per_au<<" fs"<<"\n\r";
    }
    else 
        plasma_header << "[ Grid ] Using static grid" << "\n\r";

    plasma_header<<"[ Rate Solver ] Using initial time step size of "<<this->dt*Constant::fs_per_au<<" fs"<<"\n\r";
    plasma_header<<banner<<"\n\r";

    steps_per_grid_transform =  round(num_steps*(grid_update_period/(simulation_end_time-simulation_start_time)));


    std::cout << plasma_header.str()<<std::flush; // display in regular terminal, so that it is still visible after end of program

    #ifdef NCURSES
    Display::header += plasma_header.str(); // display this in ncurses screen
    //Display::deactivate();
    #endif


    // Set up display
    std::stringstream tol;
    tol << "[ sim ] Implicit solver uses relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<"\n\r";
    std::cout << tol.str();  // Display in regular terminal even after ncurses screen is gone.
    Display::header += tol.str(); 
    Display::display_stream = std::stringstream(Display::header, ios_base::app | ios_base::out); // text shown that updates with frequency 1/steps_per_time_update.
    Display::popup_stream = std::stringstream(std::string(), ios_base::app | ios_base::out);  // text displayed during step of special events like grid updates
    Display::create_screen(); 

    // Call hybrid integrator to iterate through the time steps (good state)
    good_state = true;
    assert(num_steps == round((simulation_end_time - simulation_start_time)/this->dt + 1));
    steps_per_time_update = max(1 , (int)(input_params.time_update_gap/this->dt)); 
    this->solve_dynamics(_log,simulation_start_time, simulation_resume_time, steps_per_time_update); // Inherited from ABM

    if (!good_state) {
        _log<<"[ CRITICAL ERROR ] For an unexpected reason, a bad state was not, or ceased being, attempted to be resolved with finer timesteps."<<endl;
        throw std::runtime_error("For an unexpected reason, a bad state was not, or ceased being, attempted to be resolved with finer timesteps.");
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
        if (input_params.Store[a].bound_free_excluded)
            continue;
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

    if (!good_state) return;


    // Bound and bound-free interactions
    for (size_t a = 0; a < s.atomP.size(); a++) {
        auto t9 = std::chrono::high_resolution_clock::now();
        const bound_t& P = s.atomP[a];
        bound_t& Pdot = sdot.atomP[a];
        
        #ifdef DEBUG_BOUND
        for(size_t i=0;i < Pdot.size();i++){
            assert(Pdot[i] + P[i] >= 0);
        }
        #endif
        double old_bound_charge = sdot.bound_charge;
        // PHOTOIONISATION
        double J = pf(t); // photon flux in atomic units
        for ( auto& r : input_params.Store[a].Photo) {
            double tmp = r.val*J*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(a,r.energy, r.val*J*P[r.from]);  // TODO change to tmp?
            sdot.bound_charge +=  tmp;
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*J*P[r.from]);
        }

        sdot.cumulative_photo[a]+=sdot.bound_charge-old_bound_charge;

        #ifndef NO_ELECTRON_SOURCE
        //PHOTOION. SOURCE
        if(t < simulation_start_time + input_params.electron_source_duration*(timespan_au)){
            double injected_density = 0; 
            switch (input_params.electron_source_type){
                case 'c':
                {
                    injected_density = input_params.electron_source_fraction * y[1].cumulative_photo[a] * (J / pf(this->t[1]));
                break;
                }
                case 'p':
                {
                    injected_density = input_params.electron_source_fraction*(sdot.bound_charge-old_bound_charge);  // Additional electrons added at rate proportional to rest of sample. (Note this controls for differences in emission rate that normally would occur in the species producing photoelectrons of a different energy.)
                break;
                }
            }
            if (injected_density > 0)
                sdot.F.addDeltaSpike(a,input_params.electron_source_energy,injected_density);
        }
        #endif //NO_ELECTRON_SOURCE
        
        #ifdef RATES_TRACKING
        photo_rate.back() += sdot.bound_charge;
        #endif
        old_bound_charge = sdot.bound_charge;
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
            #ifdef RATES_TRACKING
            fluor_rate.back() +=tmp;
            #endif
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
            sdot.F.addDeltaSpike(a,r.energy, r.val*P[r.from]);
            // sdot.F.add_maxwellian(r.energy*2./3., r.val*P[r.from]);
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*P[r.from]);
            sdot.bound_charge +=  tmp;
        }
        #ifdef RATES_TRACKING
        auger_rate.back() += sdot.bound_charge - old_bound_charge;
        #endif
        old_bound_charge = sdot.bound_charge;       

        #ifdef DEBUG_BOUND
        for(size_t i=0;i < Pdot.size();i++){
            assert(Pdot[i] + P[i] >= 0);
        }
        #endif        
        // Secondary ionization
        // EII / TBR bound-state dynamics
        #ifdef TRACK_SINGLE_CONTINUUM
        size_t _c = 0; 
        #else
        size_t _c = 1;  // Iterates through the continuum corresponding to each element's initiated cascades and adds separately. // TODO check if having the full continuum contribute to the actual calcs is better.
        #endif 
        for (;_c < Distribution::num_continuums; _c++){
            Eigen::VectorXd vec_dqdt_scndry = Eigen::VectorXd::Zero(Distribution::size);
            if(input_params.Store[a].bound_free_excluded) continue;

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
                        tmp = finPair.val*s.F[_c][n]*P[init];
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
                            #ifdef TRACK_SINGLE_CONTINUUM
                            tmp = finPair.val*s.F[_c][n]*s.F[_c][m]*P[init]*2;   // _c = 0 (total continuum)
                            #else
                            tmp = finPair.val*(s.F[_c][n]*s.F[0][m] + s.F[_c][m]*s.F[0][n])*P[init]; // as we are distinguishing the electron continuum of interest from the full continuum now.
                            #endif
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
                        tmp = finPair.val*s.F[_c][n]*s.F[0][n]*P[init];
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
            #ifdef RATES_TRACKING
            eii_rate.back() += sdot_bound_charge_eii_subst;
            tbr_rate.back() += sdot_bound_charge_tbr_subst;
            #endif
        

            auto t10 = std::chrono::high_resolution_clock::now();
            pre_tbr_time += t10 - t9;

            // Free-electron parts
            #ifdef NO_EII
            #warning No impact ionisation
            #else
            auto t1 = std::chrono::high_resolution_clock::now();
            s.F.get_Q_eii(_c,vec_dqdt_scndry, a, P, threads);
            auto t2 = std::chrono::high_resolution_clock::now();
            eii_time += t2 - t1;
            #endif
            #ifdef NO_TBR
            #warning No three-body recombination
            #else
            auto t3 = std::chrono::high_resolution_clock::now();
            s.F.get_Q_tbr(_c,vec_dqdt_scndry, a, P, threads);  // Serially, this is the computational bulk of the program - S.P.
            auto t4 = std::chrono::high_resolution_clock::now();
            tbr_time += t4 - t3;
            #endif
            // Add secondary ionization to distributions
            auto t7 = std::chrono::high_resolution_clock::now();
            sdot.F.applyDeltaF(_c-1,vec_dqdt_scndry,threads);
            auto t8 = std::chrono::high_resolution_clock::now();
            apply_delta_time += t8 - t7;

        }
        // Add primary ionization to distributions
        auto t7 = std::chrono::high_resolution_clock::now();
        sdot.F.applyDeltaF(a,vec_dqdt,threads);
        sdot.F.addLoss(a,s.F, input_params.loss_geometry, s.bound_charge);
        auto t8 = std::chrono::high_resolution_clock::now();
        apply_delta_time += t8 - t7;
    }

    // if(input_params.Filtration_File() == "")
    //     // No background, use geometric-dependent loss.
    //     sdot.F.addLoss(s.F, input_params.loss_geometry, s.bound_charge);
    // else
    //     // background handles loss, apply photoelectron filtration with background
    //     sdot.F.addFiltration(s.F, s_bg.F,input_params.loss_geometry);
    
    if (isnan(s.norm(0)) || isnan(sdot.norm(0))) {
        cerr<<"NaN encountered in ODE iteration."<<endl;
        cerr<< "t = "<<t*Constant::fs_per_au<<"fs"<<endl;
        good_state = false;
        timestep_reached = t*Constant::fs_per_au;
    }    
}


// 'badly-behaved' i.e. stiff part of the system. Electron-electron interactions.
void ElectronRateSolver::sys_ee(const state_type& s, state_type& sdot) {
    sdot=0;
    #ifdef NO_EE
    #warning No electron-electron interactions
    #else
    #ifdef TRACK_SINGLE_CONTINUUM
    size_t _c = 0; 
    #else
    size_t _c = 1;  // Iterates through the continuum corresponding to each element's initiated cascades and adds separately. // TODO check if having the full continuum contribute to the actual calcs is better.
    #endif 
    for (;_c < Distribution::num_continuums; _c++){
        Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);
        const int threads = input_params.Plasma_Threads(); 
        // compute the dfdt vector
        // Electron-electon repulsions
        auto t5 = std::chrono::high_resolution_clock::now();
        s.F.get_Q_ee(_c, vec_dqdt, threads); 
        auto t6 = std::chrono::high_resolution_clock::now();
        ee_time += t6 - t5;
        // 
        // Add change to distribution
        auto t7 = std::chrono::high_resolution_clock::now();
        /*
        #ifdef TRACK_SINGLE_CONTINUUM
        sdot.F.applyDeltaF(-99,vec_dqdt,threads);  // -99 is a dummy value
        #else
        sdot.F.applyDeltaF_element_scaled(_c, vec_dqdt,threads);
        #endif
        */
        sdot.F.applyDeltaF(_c-1,vec_dqdt,threads);
        auto t8 = std::chrono::high_resolution_clock::now();
        apply_delta_time += t8 - t7;
    }
    #endif // NO_EE
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
    num_steps = ceil(t.size() + (fact - 1)*remaining_steps);
    steps_per_time_update = max(1 , (int)(0.5+input_params.time_update_gap/(this->dt))); 

    _log <<"Reloading grid"<<endl;
    // with the states loaded the spline factors are as they were, we just need to load the appropriate knots.
    assert(saved_knots == Distribution::get_knots_from_history(n));
    for(auto elem : saved_states){
        assert(elem.F.container_size() == saved_states.back().F.container_size());
    }


    // Move current step to first interpolated state step (set below)
    n-= (order-1); 
    // Linearly interpolate the previous states as a first order correction. (We should use 4th order Runge-Kutta, but I don't have time to get that working with sys_ee incorporated)
    // First element of saved_states is y[n-this->order], corresponding to t[n-this->order].
    std::vector<state_type> interpolated_states(saved_states.size());
    for(size_t i = 0; i < saved_states.size();i++){
        interpolated_states[i]=saved_states[i];
    }
    for(size_t i = 0; i < saved_states.size() - 1; i++){
        double intra_time = t[n] +i*dt; 
        size_t loop_step = n;  //
        assert(saved_times.front()==t[loop_step]);
        while(1){
            int rel_idx = loop_step - n; 
            if(saved_times[rel_idx] == intra_time){
                interpolated_states[i] = saved_states[rel_idx];
                break;
            }
            // Check if intra_time is between t[loop_step] and t[loop_step+1]
            if(saved_times[rel_idx] < intra_time){
                // if t' = t(n) + u 
                // y(t') = y(n) + [y(n+1) - y(n)] * u/[t(n+1)-t(n)]  = y(n) + dy
                // find dy
                assert(rel_idx < static_cast<int>(saved_times.size()));
                interpolated_states.at(i) = saved_states.at(rel_idx);
                interpolated_states[i] *= -1;
                interpolated_states[i] += saved_states.at(rel_idx +1);
                interpolated_states[i] *= (intra_time-saved_times[rel_idx])/(saved_times[rel_idx+1]-saved_times[rel_idx]);
                // add y(n) to get y(t')
                interpolated_states[i] += saved_states[rel_idx];
                break;
            }
            loop_step++;
            if (saved_states.size() <= loop_step - n)
                throw runtime_error("Unexpected error, could not interpolate times when loading checkpoint.");
        }
    }

    // Fill in rest of times.
    t.resize(n+1); t.resize(num_steps);
    // Set up the t grid       
    for (int i=n+1; i<num_steps; i++){
        this->t[i] = this->t[i-1] + this->dt;
    }
    steps_per_grid_transform = round(steps_per_grid_transform*fact);

    n = reload_grid(_log, n, saved_knots,interpolated_states);  // Reload grid from the START of the interpolated states
    

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
    num_steps = ceil(t.size() - (1-1/fact)*remaining_steps);
    
    t.resize(n+1); t.resize(num_steps);
    y.resize(n+1); y.resize(num_steps);
    steps_per_grid_transform = round(steps_per_grid_transform/fact);

    // Set up the t grid       
    for (int i=n+1; i<num_steps; i++){   // note potential inconsistency(?) with hybrid's solve_dynamics(): npoints = (t_final - t_initial)/this->dt + 1
        this->t[i] = this->t[i-1] + this->dt;
    }
    // */
}

#ifdef INTERACTIVE
void ElectronRateSolver::pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update){
    auto t_start = std::chrono::high_resolution_clock::now();
    

    ////// Display info ////// (only the regular stuff seen each step, not "popups" from popup_stream)
    auto t_start_disp = std::chrono::high_resolution_clock::now();
    if ((n-this->order)%steps_per_time_update == 0){
        Display::display_stream.str(Display::header); // clear display string
        Display::display_stream<< "\n\r"
        << "--- Press BACKSPACE/DEL to end simulation and save the data ---\n\r"   
        << "[ sim ] Next data backup in "<<(minutes_per_save - std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - time_of_last_save)).count()<<" minute(s).\n\r"  
        << "[ sim ] Current timestep size = "<<this->dt*Constant::fs_per_au<<" fs\n\r"   
        << "[ sim ] t="
        << this->t[n] * Constant::fs_per_au << " fs\n\r" 
        << "[ sim ] " <<Distribution::size << " knots currently active\n\r";
        //<< Distribution::get_knot_energies() << "\n\r"; 
        // << flush; 
        Display::show(Display::display_stream);
    }        
    display_time += std::chrono::high_resolution_clock::now() - t_start_disp;  

    ////// live plotting ////// 
    auto t_start_plot = std::chrono::high_resolution_clock::now();
    if ((n-this->order+1)%input_params.steps_per_live_plot_update == 0){ // TODO implement minimum time. also shld depend on num ministeps
        size_t num_pts = 4000;
        py_plotter.plot_frame(
            Distribution::get_energies_eV(num_pts),
            this->y[n].F.get_densities(0,num_pts,Distribution::get_knot_energies()), 
            Distribution::get_trimmed_knots(Distribution::get_knot_energies())
        );
    }        
    plot_time += std::chrono::high_resolution_clock::now() - t_start_plot;  

    ////// Dynamic time updates //////
    auto t_start_dyn_dt = std::chrono::high_resolution_clock::now();
    // store checkpoint
    if ((n-this->order+1)%checkpoint_gap == 0){  //
        assert(n>this->order);
        old_checkpoint = checkpoint;
        // Find the nearest streak of states in the same knot basis. 
        // This looks for the most recent grid update, BEFORE the current step.
        // if the knot was changed on the current step, stop, as we are only saving the previous steps (which will not have been transformed), and this will correspond to the same point that was stopped by the corresponding previous call of this code. 
        size_t checkpoint_n = n;
        while (static_cast<int>(checkpoint_n) - static_cast<int>(Distribution::most_recent_knot_change_idx(checkpoint_n-1)) < static_cast<int>(this->order)) {
            checkpoint_n--;
            assert(checkpoint_n > this->order);
            assert(static_cast<int>(checkpoint_n) > static_cast<int>(n) - static_cast<int>(this->order)-1); // may trigger if knot changes too close together.
        }           
        std::vector<state_type> check_states;
        std::vector<double> check_times;
        {
            std::vector<state_type>::const_iterator start_vect_idx = y.begin() - order + 1 + checkpoint_n;  
            std::vector<state_type>::const_iterator end_vect_idx = y.begin() + checkpoint_n;  
            check_states = std::vector<state_type>(start_vect_idx, end_vect_idx+1);
        }
        {
            std::vector<double>::const_iterator start_vect_idx = t.begin() - order + 1 + checkpoint_n;  
            std::vector<double>::const_iterator end_vect_idx = t.begin() + checkpoint_n;  
            check_times = std::vector<double>(start_vect_idx, end_vect_idx+1);
        }        
        assert(check_states.size() == order);
        assert(check_times.size() == order);
        assert(check_states.front().F.container_size() == check_states.back().F.container_size());
        checkpoint = {checkpoint_n, Distribution::get_knot_energies(),this->regimes,check_states,check_times};
    }
    if (euler_exceeded || !good_state){     
        bool store_dt_increase = true;   
        if (euler_exceeded){    
            if (!increasing_dt){
                Display::popup_stream <<"\nMax euler iterations exceeded, reloading checkpoint at "<<this->t[old_checkpoint.n]*Constant::fs_per_au<< " fs and decreasing dt.\n\r";
            }
            else{
            Display::popup_stream << "\nReloading checkpoint at "<<t[old_checkpoint.n]*Constant::fs_per_au
            <<" fs, as increasing step size led to euler iterations being exceeded. Will attempt again after "<<2*checkpoint_gap<< "steps.\n\r";
            times_to_increase_dt.back() = t[n] + 2*checkpoint_gap*dt;  
            store_dt_increase = false;
            increasing_dt = 0;                
            }
        }
        else{
            Display::popup_stream <<"\n NaN encountered in ODE iteration. Reloading checkpoint at "<<t[old_checkpoint.n]*Constant::fs_per_au<< " fs and decreasing dt.\n\r";
        }
        Display::show(Display::display_stream,Display::popup_stream);
        // Reload at checkpoint's n, updating the n step, and decreasing time step length
        n = load_checkpoint_and_decrease_dt(_log,n,old_checkpoint);  // virtual function overridden by ElectronRateSolver
        
        if (!store_dt_increase){
            // remove the extra time added.   
            times_to_increase_dt.pop_back();           
        }
        good_state = true;
        euler_exceeded = false;
        checkpoint = old_checkpoint;
    }
    /* Disabled for now.
    else if (!times_to_increase_dt.empty() && times_to_increase_dt.back() < t[n]){
        Display::popup_stream <<"\nAttempting to increase dt\n\r";
        Display::show(Display::display_stream,Display::popup_stream);
        if (increasing_dt > 0){
            increasing_dt--;
        
            if (increasing_dt == 0){
                //successfully increased step size
                times_to_increase_dt.pop_back();
            }
        }
        else{
            // Attempt to increase step size (i.e. decrease num remaining steps)
            this->increase_dt(_log,n);
            // Flag that we are testing to see if still converging.
            size_t num_steps_to_check = 20;
            assert(num_steps_to_check < checkpoint_gap);
            increasing_dt = num_steps_to_check;
        }
    }
    */
    dyn_dt_time += std::chrono::high_resolution_clock::now() - t_start_dyn_dt;
    
    
    /// Save data periodically in case of crash (currently need to manually copy mol file from log folder and move folder to output if want to plot). 
    if(std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - time_of_last_save) > minutes_per_save){
        auto t_start_backup = std::chrono::high_resolution_clock::now();
          _log << "Saving output backup" << endl;
          Display::popup_stream <<"\nSaving output backup, simulation will resume soon...\n\r";
          Display::show(Display::display_stream,Display::popup_stream);        
          time_of_last_save = std::chrono::high_resolution_clock::now();
          size_t old_size = y.size();
          y.resize(n+1); t.resize(n+1);
          std::cout.setstate(std::ios_base::failbit);  // disable character output
          save(data_backup_folder);
          std::cout.clear(); // enable character output
          y.resize(old_size); t.resize(old_size);
          // re-set the future t values       
          for (size_t i=n+1; i<old_size; i++){   // note potential inconsistency(?) with hybrid's solve_dynamics(): npoints = (t_final - t_initial)/this->dt + 1
              this->t[i] = this->t[i-1] + this->dt;
          }          
        backup_time += std::chrono::high_resolution_clock::now() - t_start_backup;
    }
    
    
    pre_ode_time += std::chrono::high_resolution_clock::now() - t_start;  
}

int ElectronRateSolver::post_ode_step(ofstream& _log, size_t& n){
    auto t_start = std::chrono::high_resolution_clock::now();

    //////  Dynamic grid updater ////// 
    #ifndef SWITCH_OFF_ALL_DYNAMIC_UPDATES
    auto t_start_grid = std::chrono::high_resolution_clock::now();
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic && (n-this->order+1)%steps_per_grid_transform == 0){ // TODO if adaptive time step algo is improved would be good to have a variable that this is equal to that is modified to account for changes in time step size. If a dt decreases you push back the grid update. If you increase dt (which currently doesn't happen) you could 'miss' it .
        Display::popup_stream << "\n\rUpdating grid... \n\r"; 
        _log << "[ Dynamic Grid ] Updating grid" << endl;
        Display::show(Display::display_stream,Display::popup_stream);  
        update_grid(_log,n+1,false);
        if (Distribution::reset_on_next_grid_update){
            dyn_grid_time += std::chrono::high_resolution_clock::now() - t_start_grid;  
            Distribution::reset_on_next_grid_update = false;
            reinitialise_solver_with_current_grid(_log);    
            return 1;        
        } 
    }   
    // move from initial grid to dynamic grid shortly after a fresh simulation's start.
    /*
    else if (n-this->order == max(2,(int)(steps_per_grid_transform/10)) && (input_params.Load_Folder() == "") && !grid_initialised){  // TODO make this an input param
        Display::popup_stream << "\n\rPerforming initial grid update... \n\r"; 
        _log << "[ Dynamic Grid ] Performing initial grid update..." << endl;
        Display::show(Display::display_stream,Display::popup_stream); 
        update_grid(_log,n+1,true); 
        //////// 
        // TODO restart simulation with this better grid.
        ////////
    }
    */
    dyn_grid_time += std::chrono::high_resolution_clock::now() - t_start_grid;  
    #endif //SWITCH_OFF_ALL_DYNAMIC_UPDATES
    
    //////  Check if user wants to end simulation early ////// 
    auto t_start_usr = std::chrono::high_resolution_clock::now();
    auto ch = wgetch(Display::win);
    if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){   
        flushinp(); // flush buffered inputs
        Display::popup_stream <<"\n\rExiting early... press backspace/del again to confirm or any other key to cancel and resume the simulation \n\r";
        Display::show(Display::display_stream,Display::popup_stream);
        nodelay(Display::win, false);
        ch = wgetch(Display::win);  // note implicitly refreshes screen
        if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){
            y.resize(n);
            t.resize(n);    
            return 1;
        }
        nodelay(Display::win, true);
        Display::show(Display::display_stream);
    }     
    
    user_input_time += std::chrono::high_resolution_clock::now() - t_start_usr;  
    
    if (Display::popup_stream.rdbuf()->in_avail() != 0){
        // clear stream
        Display::popup_stream.str(std::string());
    }   
    
    
    post_ode_time += std::chrono::high_resolution_clock::now() - t_start;    
    return 0;
}
#endif

void ElectronRateSolver::update_grid(ofstream& _log, size_t latest_step, bool force_update){
    size_t n = latest_step;
    //// Set up new grid ////        
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    // Latest step is n, so we decide our new grid based on that step, then transform N = "order" of the prior points to the new basis.
    set_up_grid_and_compute_cross_sections(_log,false,n,force_update); // virtual function overridden by ElectronRateSolver
    std::cout.clear();
    //// We need some previous points needed to perform next ode step, so transform them to new basis ////
    std::vector<double> new_energies = Distribution::get_knot_energies();
                
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
// first_step = next_ode_states_used[0] 
//
// Does not resize containers
size_t ElectronRateSolver::reload_grid(ofstream& _log, size_t& load_step, std::vector<double> knots, std::vector<state_type> next_ode_states_used){
    assert(next_ode_states_used.size() == order);  



    size_t n = load_step;
    assert(y[n].atomP == next_ode_states_used[0].atomP);

    std::cout.setstate(std::ios_base::failbit);  // disable character output
    
    bool rates_uninitialised = false;
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        
        // // Clear knot history up to this step.
        // If this is turned on (may be useful if change made to update grid on checkpoint load, need to remove lines in HybridIntegrator.hpp and ElectronRateSolver.cpp corresponding to comment: // may trigger if knot changes too close together.
        // while(Distribution::knots_history.back().step >= n && Distribution::knots_history.size() > 0){
        //     Distribution::knots_history.pop_back();
        // }           

        // Clear knot history past this step.
        while(Distribution::knots_history.back().step > n && Distribution::knots_history.size() > 0){
            Distribution::knots_history.pop_back();
            rates_uninitialised = true;
        }              
        // set basis to one given by knots
        if (rates_uninitialised)
            Distribution::set_basis(0, param_cutoffs, regimes, knots, false);
    }


    // Load distributions and also the next few states that we need for the first ode step.
    y.resize(n); // Removes all steps corresponding to and past the steps to load    
    for(state_type state : next_ode_states_used){
        y.push_back(state);
        n++;
    }    
    n-=1; // Set step back to the last updated step.
    y.resize(num_steps);  
    initialise_transient_y((int)n);
   

    if (rates_uninitialised == false){
        if(_log.is_open()){
            _log << "------------------- [ Reloaded Step ] -------------------\n" 
            "Time: "<<t[load_step]*Constant::fs_per_au <<"-"<<t[n]*Constant::fs_per_au<<" fs; "<<"Step: "<<load_step<<"-"<<n<<"\n" 
            << endl;
        }           
        std::cout.clear();
        return n;
    }
    //////// Update rates and container sizes etc. as necessary for when grid is changed ///////////
    {
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        if(_log.is_open()){
            double e = Constant::eV_per_Ha;
            _log << "------------------- [ Reloaded Step + Knots ] -------------------\n" 
            "Time: "<<t[load_step]*Constant::fs_per_au <<"-"<<t[n]*Constant::fs_per_au<<" fs; "<<"Step: "<<load_step<<"-"<<n<<"\n" 
            <<"Therm [peak; range]: "<<regimes.mb_peak*e<< "; "<< regimes.mb_min*e<<" - "<<regimes.mb_max*e<<"\n"; 
            for(size_t i = 0; i < regimes.num_dirac_peaks;i++){
                _log<<"Photo [peak; range]: "<<regimes.dirac_peaks[i]*e<< "; " << regimes.dirac_minimums[i]*e<<" - "<<regimes.dirac_maximums[i]*e<<"\n";
            }
            _log <<"Transition energy: "<<param_cutoffs.transition_e*e<<" eV\n"  
            << "Grid size: "<<Distribution::size<<"\n"
            << "-----------------------------------------------------" 
            << endl;
        }        
    }

  
    // update the tensor of coefficients 
    compute_free_grid_rates();    

    // Set zero_y for new grid
    set_zero_y();    

    std::cout.clear();
    }
    return n;

    // The next containers are made to have the correct size, as the initial state is set to tmp=zero_y and sdot is set to an empty state. 
} 

void ElectronRateSolver::reinitialise_solver_with_current_grid(ofstream& _log){
    _log << "[ Dynamic Grid ] Restarting solver with new grid"<<endl;
    t.resize(1);
    y.resize(1);
    std::cout.setstate(std::ios_base::failbit);  // disable character output
    set_starting_state();
    std::cout.clear(); // enable character output    

    Distribution::knots_history.resize(0);
    Distribution::knots_history.push_back(indexed_knot{0,Distribution::get_knot_energies()});

    solve_dynamics(_log,simulation_start_time, simulation_start_time, steps_per_time_update);
}



void ElectronRateSolver::compute_free_grid_rates(){
    //Clear any old rates. //TODO these aren't constant now - change to lowercase.
    RATE_EII.clear();
    RATE_TBR.clear();
    RATE_EII.resize(input_params.Store.size());
    RATE_TBR.resize(input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        size_t N = Distribution::num_basis_funcs();
        RATE_EII[a].clear();
        RATE_TBR[a].clear();
        if (input_params.Store[a].bound_free_excluded){
            RATE_EII[a].resize(0);
            RATE_TBR[a].resize(0);
        }
        else{
            RATE_EII[a].resize(N);
            RATE_TBR[a].resize(N*(N+1)/2);
        }
    }
    precompute_gamma_coeffs();
    Distribution::precompute_Q_coeffs(input_params.Store);    
}

#ifdef INTERACTIVE
//IOFunctions found in IOFunctions.cpp
void ElectronRateSolver::pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update){
    auto t_start = std::chrono::high_resolution_clock::now();
    
    #ifdef NCURSES
    ////// Display info ////// (only the regular stuff seen each step, not "popups" from popup_stream)
    auto t_start_disp = std::chrono::high_resolution_clock::now();
    if ((n-this->order)%steps_per_time_update == 0){
        Display::display_stream.str(Display::header); // clear display string
        Display::display_stream<< "\n\r"
        << "--- Press BACKSPACE/DEL to end simulation and save the data ---\n\r"   
        << "[ sim ] Next data backup in "
        <<(minutes_per_save - std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - time_of_last_save)).count()<<" minute(s).\n\r"  
        << "[ sim ] Current timestep size = "<<this->dt*Constant::fs_per_au<<" fs\n\r"   
        << "[ sim ] t="
        << this->t[n] * Constant::fs_per_au << " fs\n\r" 
        << "[ sim ] " <<Distribution::size << " knots currently active\n\r";
        //<< Distribution::get_knot_energies() << "\n\r"; 
        // << flush; 
        Display::show(Display::display_stream);
    }        
    display_time += std::chrono::high_resolution_clock::now() - t_start_disp;  
    #endif 

    #ifdef PYBIND
    ////// live plotting ////// 
    auto t_start_plot = std::chrono::high_resolution_clock::now();
    if ((n-this->order+1)%20 == 0){ // TODO implement minimum time. also shld depend on num ministeps
        size_t num_pts = 4000;
        py_plotter.plot_frame(
            Distribution::get_energies_eV(num_pts),
            this->y[n].F.get_densities(num_pts,Distribution::get_knot_energies()), 
            Distribution::get_trimmed_knots(Distribution::get_knot_energies())
        );
    }        
    plot_time += std::chrono::high_resolution_clock::now() - t_start_plot;  
    #endif //PYBIND

    ////// Dynamic time updates //////
    auto t_start_dyn_dt = std::chrono::high_resolution_clock::now();
    // store checkpoint
    if ((n-this->order+1)%checkpoint_gap == 0){  //
        old_checkpoint = checkpoint;

        // Find the nearest streak of states in the same knot basis. 
        // This looks for the most recent grid update, BEFORE the current step.
        // if It occurred ON the current step, then we are fine, unless there was another knot change in the last few steps. (Which shouldn't happen but I haven't completely verified I've caught all cases.)
        size_t checkpoint_n = n;
        while (Distribution::most_recent_knot_change_idx(checkpoint_n-1) >= n - this->order + 1){
            checkpoint_n--;
            assert(checkpoint_n > this->order);
            assert(int(checkpoint_n) > int(n) - this->order-1); // may trigger if knot changes too close together.
        }           
        std::vector<state_type> check_states;
        std::vector<double> check_times;
        {
            std::vector<state_type>::const_iterator start_vect_idx = y.begin() - order + 1 + checkpoint_n;  
            std::vector<state_type>::const_iterator end_vect_idx = y.begin() + checkpoint_n;  
            check_states = std::vector<state_type>(start_vect_idx, end_vect_idx+1);
        }
        {
            std::vector<double>::const_iterator start_vect_idx = t.begin() - order + 1 + checkpoint_n;  
            std::vector<double>::const_iterator end_vect_idx = t.begin() + checkpoint_n;  
            check_times = std::vector<double>(start_vect_idx, end_vect_idx+1);
        }        
        assert(check_states.size() == order);
        assert(check_times.size() == order);
        assert(check_states.front().F.container_size() == check_states.back().F.container_size());
        checkpoint = {checkpoint_n, Distribution::get_knot_energies(),this->regimes,check_states,check_times};
    }
    if (euler_exceeded || !good_state){     
        bool store_dt_increase = true;   
        if (euler_exceeded){    
            if (!increasing_dt){
                Display::popup_stream <<"\nMax euler iterations exceeded, reloading checkpoint at "<<this->t[old_checkpoint.n]*Constant::fs_per_au<< " fs and decreasing dt.\n\r";
            }
            else{
            Display::popup_stream << "\nReloading checkpoint at "<<t[old_checkpoint.n]*Constant::fs_per_au
            <<" fs, as increasing step size led to euler iterations being exceeded. Will attempt again after "<<2*checkpoint_gap<< "steps.\n\r";
            times_to_increase_dt.back() = t[n] + 2*checkpoint_gap*dt;  
            store_dt_increase = false;
            increasing_dt = 0;                
            }
        }
        else{
            Display::popup_stream <<"\n NaN encountered in ODE iteration. Reloading checkpoint at "<<t[old_checkpoint.n]*Constant::fs_per_au<< " fs and decreasing dt.\n\r";
        }
        Display::show(Display::display_stream,Display::popup_stream);
        // Reload at checkpoint's n, updating the n step, and decreasing time step length
        n = load_checkpoint_and_decrease_dt(_log,n,old_checkpoint);  // virtual function overridden by ElectronRateSolver
        
        if (!store_dt_increase){
            // remove the extra time added.   
            times_to_increase_dt.pop_back();           
        }
        good_state = true;
        euler_exceeded = false;
        checkpoint = old_checkpoint;
    }
    dyn_dt_time += std::chrono::high_resolution_clock::now() - t_start_dyn_dt;
    
    
    /// Save data periodically in case of crash (currently need to manually move files to a folder in output if want to resume simulation; as well as copy mol file to it from output/log/ if want to plot). 
    if(std::chrono::duration_cast<std::chrono::minutes>(std::chrono::high_resolution_clock::now() - time_of_last_save) > minutes_per_save){
        auto t_start_backup = std::chrono::high_resolution_clock::now();
          _log << "Saving output backup" << endl;
          Display::popup_stream <<"\nSaving output backup, simulation will resume soon...\n\r";
          Display::show(Display::display_stream,Display::popup_stream);        
          time_of_last_save = std::chrono::high_resolution_clock::now();
          size_t old_size = y.size();
          y.resize(n+1); t.resize(n+1);
          save(data_backup_folder);
          y.resize(old_size); t.resize(old_size);
          // re-set the future t values       
          for (int i=n+1; i<old_size; i++){   // note potential inconsistency(?) with hybrid's iterate(): npoints = (t_final - t_initial)/this->dt + 1
              this->t[i] = this->t[i-1] + this->dt;
          }          
        backup_time += std::chrono::high_resolution_clock::now() - t_start_backup;
    }
    
    
    pre_ode_time += std::chrono::high_resolution_clock::now() - t_start;  
}

int ElectronRateSolver::post_ode_step(ofstream& _log, size_t& n){
    auto t_start = std::chrono::high_resolution_clock::now();

    //////  Dynamic grid updater ////// 
    #ifndef SWITCH_OFF_ALL_DYNAMIC_UPDATES
    auto t_start_grid = std::chrono::high_resolution_clock::now();
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic && (n-this->order+1)%steps_per_grid_transform == 0){ // TODO would be good to have a variable that this is equal to that is modified to account for changes in time step size. If a dt decreases you push back the grid update. If you increase dt you could miss it.
        Display::popup_stream << "\n\rUpdating grid... \n\r"; 
        _log << "[ Dynamic Grid ] Updating grid" << endl;
        Display::show(Display::display_stream,Display::popup_stream);   
        update_grid(_log,n+1,false);       
    }   
    // move from initial grid to dynamic grid shortly after a fresh simulation's start.
    /*
    else if (n-this->order == max(2,(int)(steps_per_grid_transform/10)) && (input_params.Load_Folder() == "") && !grid_initialised){  // TODO make this an input param
        Display::popup_stream << "\n\rPerforming initial grid update... \n\r"; 
        _log << "[ Dynamic Grid ] Performing initial grid update..." << endl;
        Display::show(Display::display_stream,Display::popup_stream); 
        update_grid(_log,n+1,true); 
        //////// 
        // TODO restart simulation with this better grid.
        ////////
    }
    */
    dyn_grid_time += std::chrono::high_resolution_clock::now() - t_start_grid;  
    #endif //SWITCH_OFF_ALL_DYNAMIC_UPDATES
    
    #ifdef NCURSES
    //////  Check if user wants to end simulation early ////// 
    auto t_start_usr = std::chrono::high_resolution_clock::now();
    auto ch = wgetch(Display::win);
    if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){   
        flushinp(); // flush buffered inputs
        Display::popup_stream <<"\n\rExiting early... press backspace/del again to confirm or any other key to cancel and resume the simulation \n\r";
        Display::show(Display::display_stream,Display::popup_stream);
        nodelay(Display::win, false);
        ch = wgetch(Display::win);  // note implicitly refreshes screen
        if (ch == KEY_BACKSPACE || ch == KEY_DC || ch == 127){
            y.resize(n);
            t.resize(n);    
            return 1;
        }
        nodelay(Display::win, true);
        Display::show(Display::display_stream);
    }     
    user_input_time += std::chrono::high_resolution_clock::now() - t_start_usr;  
    #endif //NCURSES
    
    
    if (Display::popup_stream.rdbuf()->in_avail() != 0){
        // clear stream
        Display::popup_stream.str(std::string());
    }   
        
    post_ode_time += std::chrono::high_resolution_clock::now() - t_start;    
    return 0;
}

#else
void ElectronRateSolver::pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update){}

int ElectronRateSolver::post_ode_step(ofstream& _log, size_t& n){}
#endif

void ElectronRateSolver::set_zero_y(){
    // Makes a zero std::vector in a mildly spooky way 
    zero_y = get_initial_state(); // do this to make the underlying structure large enough
    zero_y *= 0.; // set it to Z E R O
}
