/**
 * @file ElectronRateSolverIO.cpp
 * @authors Spencer Passmore & Alaric Sanders 
 * @brief 
 * @note may want to change to have the raw files save with not much more fineness than the load ones.
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
#include <fstream>
#include <algorithm>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <math.h>
#include <omp.h>
// #include "config.h"
#include <filesystem>

// IO functions
void ElectronRateSolver::save(const std::string& _dir) {
    string dir = _dir; // make a copy of the const value
    dir = (dir.back() == '/') ? dir : dir + "/";

    std::cout << "[ Output ] \033[95mSaving to output folder \033[94m'"<<dir<<"'\033[95m..."<< std::endl;
    
    // output at least min_outputted_points, and above that output with spacing that is unaffected if simulation was cut off early.
    // (assumes t and y have been truncated to last point calculated when save() is called.)
    num_steps_out = max(
        input_params.Out_T_size(),
        (int)(0.5 + min_outputted_points * timespan_au/(t.back()-t.front()))
        );
    saveFree(dir+"freeDist.csv");
    saveFreeRaw(dir+"freeDistRaw.csv");
    saveKnots(dir + "knotHistory.csv");
    saveBound(dir);
    std::cout <<"\033[0m"<<std::endl;

    // Save intensity
    std::vector<double> times;
    double t_fineness = timespan_au  / num_steps_out;
    double previous_t = t[0]-t_fineness;
    int i = -1;
    while (i < static_cast<int>(t.size())-1){  //TODO make this some constructed function or something
        i++;
        if(t[i] < previous_t + t_fineness && i<= t.size()-extra_fine_steps_out){
            continue;
        }        
        times.push_back(t[i]);
        previous_t = t[i];
    }
    pf.save(times,dir+"intensity.csv");
}

void ElectronRateSolver::file_delete_check(const std::string& fname){
    if(std::filesystem::exists(fname)){
        if( remove( fname.c_str() ) != 0 ){
            perror( "[ Output ] Error saving file: could not delete existing file" );
            return;
            }
        else
            puts( "File successfully deleted" );
    }    
}

/// Save the free electron density distribution
void ElectronRateSolver::saveFree(const std::string& fname) {
    file_delete_check(fname);

    // Saves a table of free-electron dynamics to file fname
    ofstream f;
    cout << "Free: \033[94m'"<<fname<<"'\033[95m | ";
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Time (fs) | Density @ energy (eV):" <<endl;
    // We need to output to the same energies, so we choose the final knots for reference.
    std::vector<double> reference_knots = Distribution::load_knots_from_history(t.size()); 
    /* Alternative: Constant spacing reference knots:
    std::vector<double> reference_knots;
    std::vector<double>::iterator x;
    double max_e = Distribution::get_knots_from_history(0).back(); 
    double spacing = max_e/100;
    double val;
    for (x = reference_knots.begin(), val = 0; x != reference_knots.end(); ++x, val += spacing) {
        *x = val;
    }  
    */  
    f << "#           | "<<Distribution::output_energies_eV(this->input_params.Out_F_size())<<endl;
    #ifdef DEBUG
    cout << "[ Dynamic Grid ], writing densities to reference knot energies: \n";
    for (double elem : reference_knots) cout << elem *Constant::eV_per_Ha<< ' ';
    cout << endl;
    #endif  

    assert(y.size() == t.size());
    
    double t_fineness = timespan_au  / num_steps_out;
    double previous_t = t[0]-t_fineness;
    int i = -1; 
    size_t next_knot_update = 0;
    while (i <  static_cast<int>(t.size())-1){
        i++;
        if (i + order == next_knot_update or i == 0){  // EDIT: No this seems to be broken Original: The knot update indices correspond to an update at that index that also transformed the last 'order' of states (as in the solver's order).
            Distribution::load_knots_from_history(i+order);
            next_knot_update = Distribution::next_knot_change_idx(i+order);
        } 
        if(t[i] < previous_t + t_fineness && i<= t.size()-extra_fine_steps_out){
            continue;
        }
        f<<round_time(t[i]*Constant::fs_per_au)<<" "<<y[i].F.output_densities(this->input_params.Out_F_size(),reference_knots)<<endl;
        previous_t = t[i];
        
    }
    f.close();
    Distribution::load_knots_from_history(t.size()); // back to original state
}

/**
 * @brief Saves each time and corresponding B-spline coefficients. 
 * @param fname 
 */
void ElectronRateSolver::saveFreeRaw(const std::string& fname) {
    file_delete_check(fname);


    ofstream f;
    cout << "Free Raw: \033[94m'"<<fname<<"'\033[95m | ";
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Energy Knot: "<< Distribution::output_knots_eV() << endl;
    f << "# Time (fs) | Expansion Coeffs (not density)"  << endl;
    
    assert(y.size() == t.size());
    double t_fineness = timespan_au  / num_steps_out;
    double previous_t = t[0]-t_fineness;
    int i = -1; 
    size_t next_knot_update = 0;
    while (i <  static_cast<int>(t.size())-1){
        i++;
        if (i + order == next_knot_update or i == 0){
            Distribution::load_knots_from_history(i+order);
            next_knot_update = Distribution::next_knot_change_idx(i+order);
        } 
        if(t[i] < previous_t + t_fineness && i<= t.size()-extra_fine_steps_out){
            continue;
        }
        f<<round_time(t[i]*Constant::fs_per_au)<<" "<<y[i].F<<endl;  // Note that the << operator divides the factors by Constant::eV_per_Ha.
        previous_t = t[i];
    }
    f.close();
    Distribution::load_knots_from_history(t.size()); // back to original state
}


void ElectronRateSolver::saveBound(const std::string& dir) {
    // saves a table of bound-electron dynamics , split by atom, to folder dir.
    assert(y.size() == t.size());
    // Iterate over atom types
    ofstream f;    
    for (size_t a=0; a<input_params.Store.size(); a++) {
        string fname = dir+"dist_"+input_params.Store[a].name+".csv";
        file_delete_check(fname);
        
        cout << "Bound: \033[94m'"<<fname<<"'\033[95m | ";
        f.open(fname);
        f << "# Ionic electron dynamics"<<endl;
        f << "# Time (fs) | State occupancy (Probability times number of atoms)" <<endl;
        f << "#           | ";
        // Index, Max_occ inherited from MolInp
        for (auto& cfgname : input_params.Store[a].index_names) {
            f << cfgname << " ";
        }
        f<<endl;
        // Iterate over time.
        double t_fineness = timespan_au  / num_steps_out;
        double previous_t = t[0]-t_fineness;
        int i = -1;
        while (i <  static_cast<int>(t.size())-1){
            i++;
            if(t[i] < previous_t + t_fineness && i<= t.size()-extra_fine_steps_out){ 
                continue;
            }            
            // Make sure all "natom-dimensioned" objects are the size expected
            assert(input_params.Store.size() == y[i].atomP.size());
            
            f<<round_time(t[i]*Constant::fs_per_au) << ' ' << y[i].atomP[a]<<endl;   // prob. multiplied by 1./Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au
            previous_t = t[i];
        }
        f.close();           
    }
    // save rates. // DISABLED as rates not integrated with solver properly yet.
    /*
    string fname = dir+"rates.csv";
    f.open(fname);
    f << "# Densities transferred over each time period (previous_time;time]" <<endl;
    f << "# Time (fs) | photo|fluor|auger|transport|eii|tbr" <<endl;
    double t_fineness = timespan_au  / num_steps_out;
    double previous_t = t[0]-t_fineness;
    std::vector<double> density = {0,0,0,0,0,0};
    int i = -1;
    i++;
    while (i <  static_cast<int>(t.size())-1){
        i++;      
        // sum up the densities over the time gap. Rate at step i-1 corresponds to change added to step i.
        std::vector<double> tmp {photo_rate[i],fluor_rate[i],auger_rate[i],bound_transport_rate[i],eii_rate[i],tbr_rate[i]};
        for (size_t j = 0; j < density.size();j++)
            density[j] += tmp[j]*(t[i]-t[i-1]);            
        if(t[i] < previous_t + t_fineness && i<= t.size()-extra_fine_steps_out){ 
            continue;
        }
        f<<round_time(t[i]*Constant::fs_per_au) << ' ' << density<<endl;
        density = {0,0,0,0,0,0};       
    }
    f.close(); 
    */
}


void ElectronRateSolver::saveKnots(const std::string& fname) {
    file_delete_check(fname);

    ofstream f;
    cout << "Knot History: \033[94m'"<<fname<<"'\033[95m | ";
    f.open(fname);
    f << "# History of free electron grid's B-spline knot updates"<<endl;
    f << "# Time set (fs) | Energies"  << endl;

    assert(y.size() == t.size());
    size_t next_knot_update = 0;
    for (size_t i=0; i<t.size(); i++) {
        if (i + order == next_knot_update or i == 0){
            Distribution::load_knots_from_history(i+order);
            next_knot_update = Distribution::next_knot_change_idx(i+order);
            f<<t[i]*Constant::fs_per_au<<" "<<Distribution::output_knots_eV()<<endl;
        } 
    }
    f.close();
    Distribution::load_knots_from_history(t.size()); // back to original state
}


/**
 * @brief 
 * @details destructive
 */
void ElectronRateSolver::tokenise(std::string str, std::vector<double> &out, const size_t start_idx,const char delim)
{
    out.resize(0); // clear 
    std::istringstream ss(str);

    std::string s;
    size_t i = 0;
    while (std::getline(ss, s, delim)) {
        i++;
        if (i <= start_idx) continue;
        if (s.size() == 0 || (s.size() == 1 && s[0] == '\r')){continue;}  // This or part is just to deal with excel being a pain.
        double d = stold(s); 
        out.push_back(d);
    }
}

/**
 * @brief Loads the times and corresponding spline factors from a simulation's raw output into this simulation's y[i].F. 
 * The final time step is made as close to the load_time specified by input_params, without being over it or obviously divergent.   
 * @details This code loads the spline scale factors and associated times from the *raw* distribution output file. 
 * The code parasitically overrides the distribution object with the original spline basis coefficients, 
 * after which it transforms to the grid point basis submitted for *this* run via transform_basis().
 * It should be noted that transform_basis() is an *approximation* of the non-polynomial behaviourpyt, albeit it is very good due to using 
 * order 64 gaussian quadrature.
 * @todo there's no need to convert the basis of old lines, if we associate distribution at each time step to a basis.
 */
void ElectronRateSolver::loadFreeRaw_and_times() {
    vector<string> time_and_BS_factors;
    const std::string& fname = input_params.Load_Folder() + "freeDistRaw.csv";

    cout << "\n[ Free ] Applying free distribution from file with order 64 gaussian quadrature: "<<fname<<"..."<<endl;
    cout << "[ Caution ] Ensure same atomic input files are used!"<<endl;
    
    ifstream infile(fname);
	if (infile.good())
        std::cout<<"Opened successfully!"<<endl;
    else {
        std::cerr<<"Opening failed."<<endl;
		exit(EXIT_FAILURE); // Quit with a huff.
        return;
    }    
    
    // READ  
    string comment = "#";
    string grid_point_flag = "# Energy Knot:";
        
     // Get indices of lines to load
    std::string line;
    double previous_t;
    double t_fineness = pow(10,-loading_t_precision);// Maximum fineness allowed, since saved simulations restrict output step fineness up until last few steps, this is effectively turned off unless last few steps are incredibly fine. 
    vector<int> step_indices;
    int i = - bspline_order - 1;
    while (std::getline(infile, line)){
        i++;
        if(i >= 0){
            std::istringstream s(line);
            string str_time;
            s >> str_time;    
            double t = convert_str_time(str_time);     
            if(i >= 1 && t < previous_t + t_fineness){ // && t < input_params.Load_Time_Max()*Constant::fs_per_au-20*t_fineness){
                continue;
            }        
            step_indices.push_back(i);
            previous_t = t;
            if (t > input_params.Load_Time_Max())
                break;
        }
    }     

    if (input_params.Using_Input_Timestep() == false){   // Get dt at end of simulation.
        infile.clear();  
        infile.seekg(0, std::ios::beg);        
        // last_idx = i
        int j = -4;
        double dt = INFINITY;
        while (std::getline(infile, line)){
            std::istringstream s(line);
            string str_time;
            s >> str_time;
            j++;
            if (j == i-1){
                dt = stod(str_time);     // prev_time
            }
            if (j == i){
                dt = (stod(str_time) - dt)/Constant::fs_per_au; // last_time - prev_time
            }
        }
        assert(dt < timespan_au);
        input_params.num_time_steps = std::round(this->timespan_au/dt);
    }
    this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, IVP_step_tolerance); // Set up so we can load, but in load_bound we set zero_y to the correct basis.
    steps_per_time_update = max(1 , (int)(input_params.time_update_gap/(timespan_au/input_params.num_time_steps))); 


    int num_steps = step_indices.size();    
    // get each raw line
    infile.clear();  
    infile.seekg(0, std::ios::beg);
    i = -4;
    std::vector<double> saved_knots;
	while (!infile.eof())
	{    
        i++;
        string line;
        getline(infile, line);
         // KNOTS
        if (!line.compare(0,grid_point_flag.length(),grid_point_flag)){
            // Remove boilerplate
            line.erase(0,grid_point_flag.size()+1);
            // Get grid points (knots)
            this->tokenise(line,saved_knots);
            // Convert to Atomic units
            for (size_t k = 0; k < saved_knots.size();k++){
                saved_knots[k]/=Constant::eV_per_Ha; 
            }
            continue;
        }
        else if (!line.compare(0, 1, comment)){
            continue;
        }
        if (!line.compare(0, 1, "")) continue;
        // continue if i not in step_indices
        if (std::find(step_indices.begin(), step_indices.end(), i) == step_indices.end()) continue;
        time_and_BS_factors.push_back(line);
    }

    // Populate solver with distributions at each time step
    t.resize(num_steps,0);  
    y.resize(num_steps,y[0]);    
    bound_t saved_time(num_steps,0);  // this is a mere stand-in for t at best.  // WHaT doES THIS MeaN?
    i = 0;
    std::vector<double> saved_f;  // Defined outside of scope so that saved_f will be from the last time step)
    for(const string &elem : time_and_BS_factors){
        // TIME
        std::istringstream s(elem);
        string str_time;
        s >> str_time;
        saved_time[i] = convert_str_time(str_time);                

        if(saved_time[i] > input_params.Load_Time_Max() || i >= y.size()){
            // time is past the maximum
            y.resize(i); // (Resized later by integrator for full sim.)
            t.resize(i);
            break;
        }
        // SPLINE FACTORS
        this->tokenise(elem,saved_f);
        saved_f.erase(saved_f.begin());    
        // Convert to old scale.
        for(size_t j = 0; j < saved_f.size();j++){
            saved_f[j] *= Constant::eV_per_Ha;
        }
        // FREE DISTRIBUTION    
        std::vector<double> new_knots;
        if (input_params.elec_grid_type.mode != GridSpacing::dynamic){
            new_knots =  y[0].F.get_knot_energies();  
            y[i].F.set_distribution_STATIC_ONLY(saved_knots,saved_f);  
            // Set knots manually
            // To ensure compatibility, transform old distribution to new grid points.    
            // TODO should remove transforming basis unless last step (and check still works)
            y[i].F.transform_basis(new_knots); 
        }    
        
        else{
            y[i].F.set_spline_factors(saved_f);
            //Distribution::set_knot_history(0,saved_knots);
            //this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, IVP_step_tolerance);
            // Starting knots are given by loadKnots()
        }
        t[i] = saved_time[i];
        i++;
    }

    // Translate time to match input of THIS run. Disabled since it isn't implemented in the bound state loading.
    if ( !( simulation_start_time - pow(10,loading_t_precision) <= saved_time[0] && saved_time[0] <= simulation_start_time + pow(10,loading_t_precision) ) ){
        //#ifndef TIME_TRANSLATION_LOADING_ALLOWED
        std::runtime_error("Start time set for this simulation doesn't match first time in data of simulation state being loaded.");
        //#endif
        // for(size_t i = 0; i < saved_time.size(); i++){
        //     t[i] += simulation_start_time - saved_time[0];
        // }
    }
    // Shave off end until get a starting point that isn't obviously divergent.
    auto too_large = [](vector<state_type>& y){return y.end()[-1].F(0) > y.end()[-2].F(0);};
    auto too_small = [](vector<state_type>& y){return y.end()[-1].F(0) > 0;};
    while (too_large(y) || too_small(y)){
        assert(y.size()>0);
        y.resize(y.size()-1);
        t.resize(t.size()-1);
    }
    // // Clear obviously divergent F
    // for(size_t i = 0;i++;i < y.size()){
    //     if (too_large(y) || too_small(y)){
    //         y[i].F = 0;
    //     }
    // }  

    if (input_params.elec_grid_type.mode != GridSpacing::dynamic){ //TODO for this to work for dynamic will need to integrate knot loading with this function rather than be separate.
        vector<double> starting_knots = Distribution::get_knot_energies(); // The knots this simulation will start with.
        string col = "\033[95m"; string clrline = "\033[0m\n";
        cout <<  col + "Loaded simulation. Param comparisons are as follows:" << clrline
        << "Grid points at simulation resume time (" + col << t.back()*Constant::fs_per_au <<"\033[0m):" << clrline
        << "gp i        | (energy , spline factor)  " << clrline;
        vector<size_t> out_idxs = {0,1, 10,100};
        for(size_t j : out_idxs){
            if (j >= y.size()) continue;  //use setw
            cout << "Source gp "<<j <<" | " + col
            << "(" << saved_knots[j] << " , " << saved_f[j] << ")" << clrline
            << "New gp "<<j <<"   | " + col                                              
            << "(" << starting_knots[j] << " , " << y.back().F[j] << ")" 
            << clrline << "------------------" << clrline;
        }
        cout << endl;
    }
    else{
        string col = "\033[95m"; string clrline = "\033[0m\n";
        cout <<  col + "Loaded simulation. Param comparisons are as follows:" << clrline
        << "Grid points at simulation resume time (" + col << t.back()*Constant::fs_per_au <<"\033[0m):" << clrline
        << "gp i        | (energy , spline factor)  " << clrline
        << "<Not yet implemented for dynamic grid>"
        << clrline << "------------------" << clrline 
        << endl;        
    }

}   
void ElectronRateSolver::loadKnots() {
    // Ensure clear knot history
    while (Distribution::knots_history.size() > 0){
        Distribution::knots_history.pop_back();
        cout<<"Cleared knot from history"<<endl;
    }

    const std::string& fname = input_params.Load_Folder() + "knotHistory.csv";

    cout << "\n[ Dynamic Grid ] Loading knots from file path: "<<fname<<"..."<<endl;
    
    ifstream infile(fname);
	if (infile.good())
        std::cout<<"Opened successfully!"<<endl;
    else {
        std::cerr<<"Opening failed."<<endl;
		exit(EXIT_FAILURE); // Quit with a huff.
        return;
    }    
    
    string comment = "#";
    // get each raw line
    infile.clear();  
    infile.seekg(0, std::ios::beg);
    int i = bspline_order - 1;
    std::vector<double> saved_knots;
    //std::vector<indexed_knot> Distribution::knots_history
	while (!infile.eof())
	{    
        i++;
        string line;
        getline(infile, line);
        if (!line.compare(0, 1, comment)) continue;
        if (!line.compare(0, 1, "")) continue;
        // TIME
        std::istringstream s(line);
        string str_time;
        s >> str_time;
        double time = convert_str_time(str_time);
        size_t step;
        bool found_step = false;
        for (size_t n=0; n < this->t.size(); n++){
            if (this->t[n] == time){
                step = n;
                found_step = true;
                break;
            }            
            if (this->t[n] > time && n!= 0){
                step = n-1;
                found_step = true;
                break;
            }
        }        
        if (!found_step)
            break;
        // KNOTS
        // Get grid points (knots)
        this->tokenise(line,saved_knots,1);
        // Convert to Atomic units
        for (size_t k = 0; k < saved_knots.size();k++){
            saved_knots[k]/=Constant::eV_per_Ha; 
        }
        Distribution::knots_history.push_back(indexed_knot{step,saved_knots});
    }
    y[i].F.load_knots_from_history(y.size()-1);
}

/// Ensuring we have consistent removal of decimal places (should need to be called for loadFreeRaw, loadBound, and loadKnots)
/// Necessary due to high precision of Constant::fs_per_au
double ElectronRateSolver::convert_str_time(string str_time){
    double time = stod(str_time);
    // Convert to right units (based on saveFreeRaw)
    time /= Constant::fs_per_au;
    // Dodge floating point errors
    return round_time(time);
}
double ElectronRateSolver::round_time(double time,const bool ceil, const bool floor){
    assert(!(ceil && floor));
    double P = pow(10,loading_t_precision);
    if (ceil)
        return std::ceil(time * P)/P;
    if (floor)
        return std::floor(time * P)/P;
    return std::round(time * P)/P;

}

/**
 * @brief Loads all times and free e densities from previous simulation's raw output, and uses that to populate y[i].F, the free distribution.
 * @attention Currently assumes same input states
 */
void ElectronRateSolver::loadBound() {
    cout << "Loading atoms' orbital states"<<endl; 
    cout << "[ Caution ] Ensure same atoms and corresponding .inp files are used!"<<endl; 


    for (size_t a=0; a<input_params.Store.size(); a++) {
        // (unimplemented) select atom's bound file  
        const std::string& fname = input_params.Load_Folder() + "dist_" + input_params.Store[a].name + ".csv"; 

        cout << "[ Free ] Loading bound states from file path: "<<fname<<"..."<<endl;
         
        ifstream infile(fname);
        // READ
        // get each raw line
        string comment = "#";
        vector<string> saved_occupancies;
        while (!infile.eof())
        {    
            string line;
            getline(infile, line);
            if (!line.compare(0, 1, comment)) continue;
            if (!line.compare(0, 1, "")) continue;
            saved_occupancies.push_back(line);
        }        
        

        int num_steps = y.size();
        if (y.size()==0){
            cout << y.size() << " <- y.size()" << endl;
            cout << "[[Dev warning]] It seems loadBound was run before loadFreeRaw_and_times, but this means loadBound won't know what times to use." << endl;
        }
        
        //  initialise - fill with initial state
        for(size_t count = 1; count < num_steps; count++){
            this->y[count].atomP[a] = this->y[0].atomP[a];
        }
        
        
        // Iterate through and find each time that matches. 
        int matching_idx;
        for(string elem : saved_occupancies){
            // TIME
            std::stringstream s(elem);
            double elem_time;
            string str_time;
            s >> str_time;            
            elem_time = convert_str_time(str_time);
            if(elem_time > t.back()){
                break;
            }            
            matching_idx = find(t.begin(),t.end(),elem_time) - t.begin(); 
            if (matching_idx >= t.size()){
                std::cerr << "Warning, mismatch in points between bound and free files!" << std::endl;
                continue;  //Error! Couldn't find a corresponding point...
            }
            else{
                std::vector<double> occ_density;
                tokenise(elem,occ_density);
                occ_density.erase(occ_density.begin()); // remove time element
                // Convert to correct units
                const double units = 1./Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;  
                for(size_t k = 0; k < occ_density.size();k++){
                    occ_density[k] /= units;
                }                 
                y[matching_idx].atomP[a] = occ_density;
            }
        }
        // // Shave time and state containers to last matching state. (Disabled, this shouldn't happen now.)
        //y.resize(matching_idx + 1);
        //t.resize(matching_idx + 1);
        if(t.size() != matching_idx + 1){
            throw std::runtime_error("No bound state found for the final loaded step or times mismatched in files."); 
        }
    }

    size_t n = t.size()-1;
    
    if (input_params.elec_grid_type.mode == GridSpacing::dynamic){
        // Detect transition energy (in lieu of it not currently being in output file)
        // todo make a func
        double dirac_peak_cutoff_density = 0; // a peak's density has to be above this to count
        dirac_energy_bounds(n,regimes.dirac_maximums,regimes.dirac_minimums,regimes.dirac_peaks,true,regimes.num_dirac_peaks,dirac_peak_cutoff_density);
        mb_energy_bounds(n,regimes.mb_max,regimes.mb_min,regimes.mb_peak,false);
        transition_energy(n, param_cutoffs.transition_e);    
        param_cutoffs.transition_e = max(param_cutoffs.transition_e,2*regimes.mb_max); // mainly for case that transition region continues to dip into negative (in which case the transition region doesn't update).   
        // Set basis
        Distribution::set_basis(n, param_cutoffs, regimes,  Distribution::get_knots_from_history(n));
        state_type::set_P_shape(input_params.Store);
        this->zero_y = get_ground_state();
        this->zero_y *= 0.; // set it to Z E R O          
    }    
}

void ElectronRateSolver::log_config_settings(ofstream& _log){
    #ifdef NO_TBR
    _log << "[ Config ] Three Body Recombination disabled in config.h" << endl;
    #endif
    #ifdef NO_EE
    _log << "[ Config ] Electron-Electron interactions disabled in config.h" << endl;
    #endif
    #ifdef NO_EII
    _log << "[ Config ] Electron-Impact ionisation disabled in config.h" << endl;
    #endif
    #ifdef BOUND_GD_HACK
    _log << "[ Config ] hacky bound transport enabled" << endl;
    #endif
//    #ifdef NO_MINISTEPS
//    _log << "[ Config ] Stiff solver intermediate steps disabled in config.h" << endl;
//    #elif defined NO_MINISTEP_UPDATING
//    _log << "[ Config ] Stiff solver's (experimental) intermediate step size updating disabled" << endl;
//    #endif
}
