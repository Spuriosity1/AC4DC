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

#include "ComputeRateParam.h"
#include "ElectronRateSolver.h"
#include "HFInputParam.h"
#include "Constant.h"
// #include "config.h"

#include "utility.h"

using namespace std;

////// The below paragraph isn't currently implemented; for now, the rates are calculated here (initialise_grid_with_computed_cross_sections()).
////// Note the computationally expensive part of the code is the solving of equations, not the rates, so 
////// that should be taken into account when considering the priority of refactoring.   - S.P.
    // Rate system solver.
    // Uses precomputed rates from AC4DC for all atomic cross-section data.
    // KEEP IN MIND:
    // - For every atom X listed in the .mol file, AC4DC must be run for the file X.inp
    // - AC4DC has input parameters for pulse width, energy and fluence.
    // - Only photon energy affects the rate calculations.
    // Let scripts/run.py handle all of these details.


int main(int argc, const char *argv[]) {
    CmdParser runsettings(argc, argv);
    if (!runsettings.valid_input) {
        return 1;
    }

    cout<<"\033[1m";
    print_banner("config/banner.txt");
    cout<<"\033[34m";
    print_banner("config/version.txt");
    cout<<"\033[0m"<<endl<<endl;

    string name, logpath, tmp_molfile, outdir;

    cout<<"Copyright (C) 2020  Alaric Sanders and Alexander Kozlov"<<endl;
    cout<<"This program comes with ABSOLUTELY NO WARRANTY; for details run `ac4dc -w'."<<endl;
    cout<<"This is free software, and you are welcome to redistribute it"<<endl;
    cout<<"under certain conditions; run `ac4dc -c' for details."<<endl;

    // Temporarily convert to string, so we can add .mol for ease of use.
    string input_file_path = string(argv[1]); 
    if (get_file_names(input_file_path, name, logpath, tmp_molfile, outdir) == 1)
        return 1;

    save_mol_file(input_file_path,tmp_molfile);

    cout<<"Running simulation for target "<<name<<endl;
    cout << "logfile name: " << logpath <<endl;
    ofstream log(logpath); 
    cout << "\033[1;32mInitialising... \033[0m" <<endl;
    const char* const_path = input_file_path.c_str();
    #ifdef PYBIND
    pybind11::initialize_interpreter();  
    #endif //PYBIND
    ElectronRateSolver S(const_path, log); // Contains all of the collision parameters.                                       
    cout << "\033[1;32mComputing cross sections... \033[0m" <<endl;
    S.set_up_grid_and_compute_cross_sections(log, true);
    if (runsettings.solve_rate_eq) {
        cout << "\033[1;32mSolving rate equations..." << "\033[35m\033[1mTarget: " << name << "\033[0m" <<endl;       
        
        string backup_dir = "output/backup_data";
        try_mkdir(backup_dir);
        S.solve(log, backup_dir);
        try_mkdir(outdir);
        S.save(outdir);    
        //pybind11::finalize_interpreter(); Commented out so is never called (via this or scoped_interpreter) since it crashes due to a missing pointer for whatever reason.
    }
    move_mol_file(tmp_molfile,outdir,name); 
    move_log_file(logpath,outdir,name);
    cout << "\033[38;5;47mDone! \033[0m" <<endl;
    return 0;
    
}
