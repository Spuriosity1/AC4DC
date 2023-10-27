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
// #include "ElectronRateSolver.h"
#include "RadialWF.h"
#include "HFInputParam.h"
#include "Constant.h"
#include <nlohmann/json.hpp>

// #include "utility.h"
#include <filesystem>

using namespace std;
using json = nlohmann::json;
using namespace InputData;

int main(int argc, const char *argv[]) {
if (argc < 3){ std::cout << "Usage: ac4dc /path/to/infile.toml /path/to/outdir/ [overrides...]" <<std::endl; return 1; }
    ifstream ifs(argv[1]);

    // ifs >> input_data;
    json input_data = json::parse(ifs,
                    /* callback */ nullptr,
                    /* allow exceptions */ true,
                    /* ignore_comments */ true);

    // Ensure that the input parameters are sensible
//    nlohmann::json_schema::json_validator validator;
//    validator.set_root_schema(HFInputParam::schema);
 //   validator.validate(input_data);

    HFInputParam in = input_data.get<HFInputParam>();

    filesystem::path outdir(argv[2]);

    

    filesystem::path outfile = outdir / (in.name + ".atomic_data");
    filesystem::path logfile = outdir / (in.name + ".log");

    cout<<"Computing atomic parameters for atom "<<in.name<<endl;
    cout << "logfile name: " << logfile <<endl;
    ofstream _log(logfile); 
    assert (_log.good());
    cout << "\033[1;32mComputing cross sections... \033[0m" <<endl;

    Grid lattice;
    lattice.logspace_from_nsteps(in.grid_min, in.grid_max, in.num_grid_points);
    Potential potential(lattice, in.nuclear_Z, in.charge_model);

   // set up the orbitals 
    std::vector<RadialWF> orbitals;
    for (auto& orb : in.orbitals){
        orbitals.push_back(RadialWF(in.num_grid_points));
        orbitals.back().set_N(orb.n);
        orbitals.back().set_L(orb.l);
        orbitals.back().set_occupancy(orb.occ);
        orbitals.back().flag_shell(); // Remeber that this is a shell (what does this mean?)
        
	    if (orb.n == 0 || orb.n > 10){
		    cerr << "[ Atomic ] \033[31;1m Incorrect number of orbital types specified in input file \033[0m" << endl;
		    cerr << "Note that a shell approximation for orbitals only counts for 1" << endl;
		    throw runtime_error("Bad atomic input");
	    }
    }
    

    HartreeFock HF(lattice, orbitals, potential, in, _log);

    ComputeRateParam Dynamics(lattice, orbitals, potential, in, _log);
    Dynamics.configure_calc(
            /* Auger */ false,
            /* Fluorescence */ false,
            /* Photo */ true,
            /* EII */ false,
            /* Fourier t'form */ true,
            /* bound */ false);

    vector<int> final_occ(orbitals.size(), 0);
    vector<int> max_occ(orbitals.size(), 0);
    vector<bool> shell_check(orbitals.size(),0); // Store the indices of shell-approximated orbitals
    for (size_t i = 0; i < orbitals.size(); i++) {
        // locks in electrons that are in a potential deeper than the (mean) photon energy
        if (fabs(orbitals[i].Energy) > in.omega) final_occ[i] = orbitals[i].occupancy();  // TODO what about continuum lowering/IPD or electrons with energies above the photon energy? -S.P.  
        max_occ[i] = orbitals[i].occupancy();
        shell_check[i] = orbitals[i].is_shell();
    }

    RateData::Atom atomic_data = Dynamics.SolvePlasmaBEB(max_occ, final_occ, shell_check);
    
    vector<vector<unsigned>> index = Dynamics.Get_Indexes();
    
    // Write the atomic rate data to file
    json j = atomic_data;
    std::ofstream ofs(outdir / (in.name + ".json") );
    assert (ofs.good());
    ofs << j;
    ofs.close();

    cout << "\033[38;5;47mDone! \033[0m" <<endl;
    return 0;
    
}
