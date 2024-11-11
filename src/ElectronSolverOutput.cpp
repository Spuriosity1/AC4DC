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
#include "config.h"
#include <filesystem>

// IO functions
void ElectronRateSolver::save(const std::string& _dir) {
    string dir = _dir; // make a copy of the const value
    dir = (dir.back() == '/') ? dir : dir + "/";

    std::cout << "[ Output ] \033[95mSaving to output folder \033[94m'"<<dir<<"'\033[95m..."<< std::endl;
    
    // output at least min_outputted_points, and above that output with spacing that is unaffected if simulation was cut off early.
    // (assumes t and y have been truncated to last point calculated when save() is called.)
    num_steps_out = max(input_params.Out_T_size(),(int)(0.5 + min_outputted_points * timespan_au/(t.back()-t.front())));
    saveFree(dir+"freeDist");
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
        if(t[i] < previous_t + t_fineness && i<= int(t.size())-extra_fine_steps_out){
            continue;
        }        
        times.push_back(t[i]);
        previous_t = t[i];
    }
    pf.save(times,dir+"intensity.csv");
}

void ElectronRateSolver::file_delete_check(const std::string& fname){
    // Commented out messages due to not playing nice with ncurses
    if(std::filesystem::exists(fname)){
        if( remove( fname.c_str() ) != 0 ){
            //perror( "[ Output ] Error saving file: could not delete existing file" );
            return;
            }
        else
           ;// puts( "File successfully deleted" );
    }    
}

/// Save the free electron density distribution
void ElectronRateSolver::saveFree(const std::string& base_fname) {
    for (size_t _c = 0; _c < Distribution::num_continuums; _c++){
        std::string fname = base_fname;
        if (_c >= 1){
            fname.append("_");
            fname.append(input_params.Store[_c-1].name); 
        }
        fname.append(".csv");
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
            if (i == static_cast<int>(next_knot_update) or i == 0){
                Distribution::load_knots_from_history(i);
                next_knot_update = Distribution::next_knot_change_idx(i);
            } 
            if(t[i] < previous_t + t_fineness && i<= int(t.size())-extra_fine_steps_out){
                continue;
            }
            f<<round_time(t[i]*Constant::fs_per_au)<<" "<<y[i].F.output_densities(_c,this->input_params.Out_F_size(),reference_knots)<<endl;
            previous_t = t[i];
            
        }
        f.close();
        Distribution::load_knots_from_history(t.size()); // back to original state
    }
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
        if (i == static_cast<int>(next_knot_update) or i == 0){
            Distribution::load_knots_from_history(i);
            next_knot_update = Distribution::next_knot_change_idx(i);
        } 
        if(t[i] < previous_t + t_fineness && i<=static_cast<int>(t.size())-extra_fine_steps_out){
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
    ofstream f;    
    // Iterate over save file type { 0: bound | 1: photoionisation | }
    for (size_t mode=0; mode < 2; mode++)
        // Iterate over atom types
        for (size_t a=0; a<input_params.Store.size(); a++) {
            string fname;
            string header;
            switch(mode){
                case 0:
                    fname = dir+"dist_"+input_params.Store[a].name+".csv";
                    header = string("# Ionic electron dynamics\n") 
                    + string("# Time (fs) | State occupancy (Probability times number of atoms)\n");
                    std::cout << "Bound: \033[94m'"<<fname<<"'\033[95m | "<<std::endl;
                break;
                case 1:
                    fname = dir+"photo_"+input_params.Store[a].name+".csv";
                    header = string("# Cumulative photoionisation\n")
                    + string("# Time (fs) | Total photoionised electron density\n");
                    std::cout << "Rates: \033[94m'"<<fname<<"'\033[95m | "<<std::endl;
                break;
                default:
                    continue;
            }
             
            file_delete_check(fname);
                        
            f.open(fname);
            f << header<<std::flush;
            if (mode == 0){
                f << "#           | ";
                // Index, Max_occ inherited from MolInp
                for (auto& cfgname : input_params.Store[a].index_names) {
                    f << cfgname << " ";
                }
                f<<endl;
            }
            // Iterate over time.
            double t_fineness = timespan_au  / num_steps_out;
            double previous_t = t[0]-t_fineness;
            int i = -1;
            while (i <  static_cast<int>(t.size())-1){
                i++;
                if(t[i] < previous_t + t_fineness && i<= int(t.size())-extra_fine_steps_out){ 
                    continue;
                }            
                switch(mode){
                    case 0:
                        // Make sure all "natom-dimensioned" objects are the size expected //TODO failsafe?
                        assert(input_params.Store.size() == y[i].atomP.size());
                        
                        f<<round_time(t[i]*Constant::fs_per_au) << ' ' << y[i].atomP[a]<<endl;   // Multiplied by 1./Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au                    
                    break;
                    case 1:                 
                        f<<round_time(t[i]*Constant::fs_per_au) << ' ' << y[i].cumulative_photo[a]<<endl;
                    break;
                    default:
                        continue;
                }
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
        if(t[i] < previous_t + t_fineness && i<= int(t.size())-extra_fine_steps_out){ 
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
        if (i == next_knot_update or i == 0){
            Distribution::load_knots_from_history(i);
            next_knot_update = Distribution::next_knot_change_idx(i);
            f<<t[i]*Constant::fs_per_au<<" "<<Distribution::output_knots_eV()<<endl;
        } 
    }
    f.close();
    Distribution::load_knots_from_history(t.size()); // back to original state
}


void ElectronRateSolver::log_config_settings(ofstream& _log){
    #ifdef NO_TBR
    _log << "[ Config ] Three Body Recombination disabled in config.h" << endl;
    cout <<  "\033[101m[ Config ] Three Body Recombination disabled in config.h\033[0m"<<endl; 
    #endif
    #ifdef NO_EE
    _log << "[ Config ] Electron-Electron interactions disabled in config.h" << endl;
    cout <<  "\033[101m[ Config ] Electron-Electron interactions disabled in config.h\033[0m"<<endl; 
    #endif
    #ifdef NO_EII
    _log << "[ Config ] Electron-Impact ionisation disabled in config.h" << endl;
    cout <<  "\033[101m[ Config ] Electron-Impact ionisation disabled in config.h\033[0m"<<endl; 
    #endif
    #ifdef BOUND_GD_HACK
    _log << "[ Config ] hacky bound transport enabled" << endl;
    #endif
    #ifdef NO_MINISTEPS
    _log << "[ Config ] Stiff solver intermediate steps disabled in config.h" << endl;
    #elif defined NO_MINISTEP_UPDATING
    _log << "[ Config ] Stiff solver's (experimental) intermediate step size updating disabled" << endl;
    #endif
}
