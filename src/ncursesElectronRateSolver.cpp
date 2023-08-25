#include "ncursesElectronRateSolver.h"

//IOFunctions found in IOFunctions.cpp
void ncursesElectronRateSolver::pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update){
    auto t_start = std::chrono::high_resolution_clock::now();
    

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

int ncursesElectronRateSolver::post_ode_step(ofstream& _log, size_t& n){
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
    #endif #SWITCH_OFF_ALL_DYNAMIC_UPDATES
    
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
