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

#include "MolInp.h"
#include "Constant.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <map>
#include <cmath>
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include <assert.h>


MolInp::MolInp(const char* filename, ofstream & _log)
{
	// Input file for molecular ionization calculation.
	map<string, vector<string>> FileContent;

	std::cout<<"Opening molecular file "<<filename<<"... ";
	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	ifstream infile(filename);
	if (infile.good())
        std::cout<<"Success!"<<endl;
    else {
        std::cerr<<"Failed."<<endl;
		exit(EXIT_FAILURE); // chuck a hissy fit and quit.
        return;
    }
	string comment = "//";
	string curr_key = "";
	string end_file = "####END####";  // Do not parse text past a line beginning with this string 

	// Store the parameter-holding file content.
	// Note that end-line comments aren't actually ignored, just we only stream first term of each line.
	while (!infile.eof())
	{
		string line;
		getline(infile, line);
		// Skip lines with no information.
		if (!line.compare(0, 2, comment)) continue;   
		if (!line.compare(0, 1, "")) continue;
		if (!line.compare(0,end_file.length(),end_file)) break; 
		// A header, defining a category of inputs
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = vector<string>(0);
			}
			while (line.back() == ' '){
				line.resize(line.size()-1);
			}
			curr_key = line;
			
		// A line below the header, defining an input
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	size_t num_atoms = FileContent["#ATOMS"].size();
	size_t num_bound_free_exclusions = FileContent["#BOUND_FREE_EXCLUSIONS"].size();

	Orbits.clear();
	Orbits.resize(num_atoms);
	Latts.clear();
	Latts.resize(num_atoms, Grid(0));
	Pots.clear();
	Pots.resize(num_atoms);
	Atomic.clear();
	Store.clear();
	Store.resize(num_atoms);
	Index.clear();
	Index.resize(num_atoms);

	for (size_t n = 0; n < FileContent["#VOLUME"].size(); n++) {
		stringstream stream(FileContent["#VOLUME"][n]);

		if (n == 0) stream >> unit_V;
		if (n == 1) stream >> loss_geometry.L0;
		if (n == 2) stream >> loss_geometry;
	}

	for (size_t n = 0; n < FileContent["#OUTPUT"].size(); n++) {
		stringstream stream(FileContent["#OUTPUT"][n]);
		if (n == 0) stream >> out_T_size;
		if (n == 1) stream >> out_F_size;
		// if (n == 2) {
		// 	stream >> tmp;
		// 	if (tmp == 'Y') write_<placeholder> = true;
		// }
	}
	
	for (size_t n = 0; n < FileContent["#PULSE"].size(); n++) {  
		stringstream stream(FileContent["#PULSE"][n]);

		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> pulse_shape;  // (Note user-defined operators)
		if (n == 3) stream >> timespan_factor;
		if (n == 4) stream >> negative_timespan_factor;
	}

	string tmp = "";
	for (size_t n = 0; n < FileContent["#USE_COUNT"].size(); n++) {  
		stringstream stream(FileContent["#USE_COUNT"][n]);

		if (n == 0) stream >> tmp;
		if (tmp[0] == 't' || tmp[0] == 'T'){
			use_count = true; 
			if (n == 1) stream >> photon_count;
		}	
	}
	tmp = "";
	for (size_t n = 0; n < FileContent["#USE_INTENSITY"].size(); n++) {  
		stringstream stream(FileContent["#USE_INTENSITY"][n]);

		if (n == 0) stream >> tmp;
		if (tmp[0] == 't' || tmp[0] == 'T'){
			use_intensity = true; 
			if (n == 1) stream >> peak_intensity;
		}		
		
	}
	tmp = "";
	for (size_t n = 0; n < FileContent["#USE_FLUENCE"].size(); n++) {  
		stringstream stream(FileContent["#USE_FLUENCE"][n]);

		if (n == 0) stream >> tmp;
		if (tmp[0] == 't' || tmp[0] == 'T'){
			use_fluence = true; 
			if (n == 1) stream >> fluence;
		}
	}


	for (size_t n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		stringstream stream(FileContent["#NUMERICAL"][n]);

		if (n == 0) stream >> num_time_steps;
		if (n == 1) stream >> omp_threads;
		if (n == 2) stream >> param_cutoffs.min_coulomb_density;
	}

	for (size_t n = 0; n < FileContent["#MANUAL_GRID"].size(); n++) {
		stringstream stream(FileContent["#MANUAL_GRID"][n]);
		if (n == 0) stream >> elec_grid_type; // "true" for manual.		
		//#GRID ,
		if (n == 1) stream >> param_cutoffs.transition_e; 
		if (n == 2) stream >> elec_grid_regions; //elec_grid_regions.bndry_idx
		if (n == 3) stream >> elec_grid_regions; //elec_grid_regions.bndry_E
		if (n == 4) stream >> elec_grid_regions; //elec_grid_regions.powers
	}
	for (size_t n = 0; n < FileContent["#DYNAMIC_GRID"].size(); n++) {
		stringstream stream(FileContent["#DYNAMIC_GRID"][n]);
		if (n == 0) stream >> elec_grid_preset;
		if (n == 1) stream >> grid_update_period;
	}	
	for (size_t n = 0; n < FileContent["#ELECTRON_SOURCE"].size(); n++) {
		stringstream stream(FileContent["#ELECTRON_SOURCE"][n]);
		if (n == 0) stream >> electron_source_fraction;
		if (n == 1) stream >> electron_source_energy;
		if (n == 2) stream >> electron_source_duration;
		if (n == 3) stream >> electron_source_type;
	}	
	for (size_t n = 0; n < FileContent["#FILTRATION"].size(); n++) {
		stringstream stream(FileContent["#FILTRATION"][n]);
		if (n == 0){ stream >> filtration_file;
			filtration_file = "output/__Molecular/" + filtration_file + "/freeDistRaw.csv";
		}
	}

	for (size_t n = 0; n < FileContent["#LOAD"].size(); n++) {
		stringstream stream(FileContent["#LOAD"][n]);

		if (n == 0){stream >> load_folder;
			load_folder = "output/__Molecular/" + load_folder + "/";
		}
		if (n == 1) stream >> simulation_resume_time_max;
		if (n == 2 && (stream.get() =='t'||stream.get() =='T')){
			loading_uses_input_timestep = true;
		} 
	}

	for (size_t n = 0; n < FileContent["#DEBUG"].size(); n++) {
		stringstream stream(FileContent["#DEBUG"][n]);

		if (n == 0){ stream >> simulation_cutoff_time; cutoff_flag = true;}
		if (n == 1) stream >> time_update_gap;
		if (n == 2) stream >> steps_per_live_plot_update;

	}	

	// Thread capping
	#if defined THREAD_MAX_1
		omp_threads = 1
	#elif defined THREAD_MAX_4
		omp_threads = min(4,omp_threads);
	#elif defined THREAD_MAX_8
		omp_threads = min(8,omp_threads);
	#elif defined THREAD_MAX_12
		omp_threads = min(12,omp_threads);
	#elif defined THREAD_MAX_16
		omp_threads = min(16,omp_threads);
	#elif defined THREAD_MAX_20
		omp_threads = min(20,omp_threads);
	#elif defined THREAD_MAX_24
		omp_threads = min(24,omp_threads);
	#elif defined THREAD_MAX_28
		omp_threads = min(28,omp_threads);		
	#endif

	// Get fluence (in 10^4 * J.cm^-2).
	if (use_count){  // [10^12 ph.µm^-2]
		fluence = photon_count*omega*Constant::J_per_eV*1e12; // [J µm^-2]
		fluence *= 1e4;  // [10^4 J cm^-2]
	}	
	if (use_intensity){
		// TODO Need to double check that this was working as expected  
		// fluence = I0 * fwhm. I0 = I_avg*timespan/(fwhm). NB: var width = fwhm
		// if square: Ipeak = I0.  
		// if gaussian, I_peak = I0/norm, where norm = sqrt(pi/4/log(2)). 
		
		// Convert units, Peak intensity is units of 10^19 W/cm^2, fluence is in 10^4 * J/cm^2.
		double t_units = pow(10,-15);   
		double I_units = pow(10,15); 
		switch (pulse_shape)
		{
		case PulseShape::gaussian:{
			const double norm = sqrt(Constant::Pi/4/log(2)); // from  Pulse::operator() 
			fluence = peak_intensity*I_units*norm*(width*t_units);
			break;
			}
		case PulseShape::square:{
			fluence = peak_intensity*I_units*(width*t_units);
			break;
			}
		default:
			throw runtime_error("Pulse shape has not been set.");
			break;
		}    
	}
	if (use_fluence){
		// already set.
	}


	
	// Give dynamic grid regions pulse energy in eV 
	elec_grid_preset.pulse_omega = omega;
	// Give energy for electron source in eV
	#ifndef NO_ELECTRON_SOURCE
	if (electron_source_energy > 0){
		elec_grid_preset.electron_source_energy = electron_source_energy;
	}
	#endif //NO_ELECTRON_SOURCE


	// Hardcode the boundary conditions
	elec_grid_type.zero_degree_inf = 3;
	elec_grid_type.zero_degree_0 = 0;

	// TODO update this.
	const string bc = "\033[33m"; // begin colour escape code
	const string clr = "\033[0m"; // clear escape code
	// 80 equals signs
	const string banner = "================================================================================";
	cout<<banner<<endl;
	cout<<bc<<"Unit cell size: "<<clr<<unit_V<<" A^3"<<endl;
	cout<<bc<<"Droplet L0:     "<<clr<<loss_geometry.L0<<" A"<<endl;
	cout<<bc<<"Droplet shape:  "<<clr<<loss_geometry<<endl<<endl;

	cout<<bc<<"Photon energy:  "<<clr<<omega<<" eV"<<endl;
	cout<<bc<<"Pulse fluence:  "<<clr<<fluence*10000<<" J/cm^2 = "<<10000*fluence/omega/Constant::J_per_eV<<"ph cm^-2"<<endl;
	cout<<bc<<"Pulse FWHM:     "<<clr<<width<<" fs"<<endl;
	cout<<bc<<"Pulse shape:    "<<clr<<pulse_shape<<endl<<endl;

	// cout<<bc<<"Electron grid:  "<<clr<<min_elec_e<<" ... "<<max_elec_e<<" eV"<<endl;
	// cout<<    "                "<<num_elec_points<<" points"<<endl;
	cout<<bc<<"Grid type:      "<<clr<<elec_grid_type<<endl;
	cout<<bc<<"Low energy cutoff for Coulomb logarithm estimation: "<<clr<<param_cutoffs.transition_e<<"eV"<<endl;
	cout<<bc<<"Minimum num electrons per unit cell for Coulomb logarithm to be considered: "<<clr<<param_cutoffs.min_coulomb_density<<endl;
	cout<<endl;

	cout<<bc<<"ODE Iteration:  "<<clr<<num_time_steps<<" timesteps"<<endl;
	if(cutoff_flag){cout<<bc<<"Simulation set to cut off early at: "<<clr<<simulation_cutoff_time<<" fs"<<endl;}
	else{cout<<bc<<"Simulation cutoff is inactive."<<clr<<endl; }	
	cout <<endl;

	cout<<bc<<"Output:         "<<clr<<out_T_size<<" time grid points"<<endl;
	cout<<    "                "<<out_F_size<<" energy grid points"<<endl << endl;
	
	cout<<bc<<"OMP threads:  "<<clr<<omp_threads<<" threads"<<endl;

	cout<<banner<<endl;
	
	// Convert to number of photon flux.
	omega /= Constant::eV_per_Ha;
	electron_source_energy/= Constant::eV_per_Ha;
	fluence *= 10000/Constant::Jcm2_per_Haa02/omega;

	// Convert to atomic units.
	width /= Constant::fs_per_au;
	simulation_cutoff_time /= Constant::fs_per_au;
	time_update_gap /= Constant::fs_per_au;
	grid_update_period /= Constant::fs_per_au;
	loss_geometry.L0 /= Constant::Angs_per_au;
	unit_V /= Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;

	param_cutoffs.min_coulomb_density /= unit_V;

	param_cutoffs.transition_e /= Constant::eV_per_Ha;
	for (size_t i = 0; i < elec_grid_regions.bndry_E.size(); i++){
		elec_grid_regions.bndry_E[i] /= Constant::eV_per_Ha;
	}

	// Reads the very top of the file, expecting input of the form
	// H 2
	// O 1
	// Then scans for atomic files of the form 
	// input/atoms/H.inp
	// input/atoms/O.inp
	// Store is then populated with the atomic data read in below.
	for (size_t i = 0; i < num_atoms; i++) {
		string at_name;
		double at_num;

		stringstream stream(FileContent["#ATOMS"][i]);
		stream >> at_name >> at_num;

		Store[i].nAtoms = at_num/unit_V;
		Store[i].name = at_name;
		// Store[i].R = radius;

		at_name = "input/atoms/" + at_name + ".inp";

		Atomic.push_back(Input((char*)at_name.c_str(), Orbits[i], Latts[i], _log));
		// Overrides pulses found in .inp files
		Atomic.back().Set_Pulse(omega, fluence, width);
		Atomic.back().Set_Num_Threads(omp_threads);

		Potential U(&Latts[i], Atomic[i].Nuclear_Z(), Atomic[i].Pot_Model());

		Pots[i] = U;
	}
		// For specified atoms turn off secondary ionisation, i.e. EII and TBR. (Useful if have only small number of a species present and photoionisation dominates their contribution to the dynamics).
		bound_free_exclusions = std::vector<bool>(num_atoms,false);
		for (size_t i = 0; i < num_bound_free_exclusions; i++){
			stringstream stream(FileContent["#BOUND_FREE_EXCLUSIONS"][i]);
			string at_name;
			stream >> at_name;
			bool found_atom = false;
			for (size_t j = 0; j < num_atoms; j++){
				if(at_name == Store[j].name){
					bound_free_exclusions[j] = true;
					found_atom = true;
					break;
				}
			}
			assert(found_atom&&"Could not find matching atom present in #ATOMS that was specified in #BOUND_FREE_EXCLUSIONS");  
		}

	if (!validate_inputs()) {
		cerr<<endl<<endl<<endl<<"Exiting..."<<endl;
		throw runtime_error(".mol input file is invalid");
	}
}

bool MolInp::validate_inputs() { // TODO need to add checks probably -S.P. TODO if it breaks it assert that num time steps greater than loaded.
	bool is_valid=true;
	cerr<<"\033[31;1m";
	if (omega <= 0 ) { cerr<<"ERROR: pulse omega must be positive"; is_valid=false; }
	if (width <= 0 ) { cerr<<"ERROR: pulse width must be positive"; is_valid=false; }
	if (fluence <= 0 ) { cerr<<"ERROR: pulse fluence must be positive"; is_valid=false; }
	if (num_time_steps <= 0 ) { cerr<<"ERROR: got negative number of timesteps"; is_valid=false; }
	if (out_T_size <= 0) { cerr<<"ERROR: system set to output zero timesteps"; is_valid=false; }
	if (out_F_size <= 0) { cerr<<"ERROR: system set to output zero energy grid points"; is_valid=false; }
	if (loss_geometry.L0 <= 0) { cerr<<"ERROR: radius must be positive"; is_valid=false; }
	if (timespan_factor < 0 || negative_timespan_factor < 0) {cerr << "ERROR, timespan factors must be postive"; is_valid = false;}
	if (timespan_factor < negative_timespan_factor){cerr << "ERROR, the timespan factor for the negative times must be smaller than the full timespan factor";is_valid=false;}
	if (pulse_shape == PulseShape::square && (timespan_factor < 1 && timespan_factor !=0)){cerr << "ERROR, timespan too short to capture full square pulse";is_valid=false;}
	if (pulse_shape == PulseShape::square && negative_timespan_factor != 0){cerr << "ERROR, timespan for negative times cannot be specified with square pulse";is_valid=false;}
	if (use_fluence + use_count + use_intensity != 1) {cerr << "ERROR, require exactly one of #USE_FLUENCE, #USE_COUNT, and #USE_INTENSITY to be active ";is_valid = false;}
	if (omp_threads <= 0) { omp_threads = 4; cerr<<"Defaulting number of OMP threads to 4"; }
	if (steps_per_live_plot_update < 1){steps_per_live_plot_update = 1; cerr<<"Steps per live plot was raised to 1 from given value of "<<steps_per_live_plot_update;}
	if (electron_source_fraction != 0 && (electron_source_energy <= 0 || electron_source_fraction < 0)) {cerr<<"Invalid electron source parameters.";is_valid=false;}

	if (elec_grid_type.mode == GridSpacing::unknown) {
	cerr<<"ERROR: Grid type not recognised - param corresponding to use of manual must start with (t)rue or (f)alse";
	is_valid=false;
	}
	if (elec_grid_type.mode == GridSpacing::manual){
		// transition e.
		if (param_cutoffs.transition_e <= elec_grid_regions.bndry_E.front() || param_cutoffs.transition_e >= elec_grid_regions.bndry_E.back()) {
			std::cout <<"Transition energy "<<param_cutoffs.transition_e*Constant::eV_per_Ha<<" outside range: "
			<< elec_grid_regions.bndry_E.front()*Constant::eV_per_Ha<<" "<<param_cutoffs.transition_e*Constant::eV_per_Ha << " " << elec_grid_regions.bndry_E.back()*Constant::eV_per_Ha << std::endl;
			param_cutoffs.transition_e = elec_grid_regions.bndry_E.back()/4;
			cerr<<"Defaulting low-energy cutoff for coulomb _log calculations to "<<param_cutoffs.transition_e*Constant::eV_per_Ha; 
		}
		// Electron grid style
		if(elec_grid_regions.bndry_E.front() < 0 || elec_grid_regions.bndry_E.back() < 0 || elec_grid_regions.bndry_E.back() <= elec_grid_regions.bndry_E.front()) { cerr<<"ERROR: Electron grid specification invalid"; is_valid=false; }
		if (num_time_steps <= 0 ) { cerr<<"ERROR: got negative number of energy steps"; is_valid=false; }	
	}

	// unit cell volume.
	if (unit_V <= 0) { cerr<<"ERROR: unit xell volume must be positive"; is_valid=false; }

	cerr<<"\033[0m";
	return is_valid;
}

// if recalc is false, and atom's secondary ionization calcs are turned off, interpolates from stored primary ionization data (photoion., and fluor. + Auger decay rates).
// Note that if the primary ionization data is missing, it will be generated. This will take a long time, but never has to be performed again
void MolInp::calc_rates(ofstream &_log, bool recalc) {
	// Loop through atomic species.
	for (size_t a = 0; a < Atomic.size(); a++) {
		HartreeFock HF(Latts[a], Orbits[a], Pots[a], Atomic[a], _log);

		// This Computes the parameters for the rate equations to use, loading them into Init.
		ComputeRateParam Dynamics(Latts[a], Orbits[a], Pots[a], Atomic[a], recalc);
		vector<int> final_occ(Orbits[a].size(), 0);
		vector<int> max_occ(Orbits[a].size(), 0);
		vector<bool> shell_check(Orbits[a].size(),0); // Store the indices of shell-approximated orbitals
		for (size_t i = 0; i < Orbits[a].size(); i++) {
			// locks in electrons that are in a potential deeper than the (mean) photon energy
			if (fabs(Orbits[a][i].Energy) > Omega()) final_occ[i] = Orbits[a][i].occupancy();   
			max_occ[i] = Orbits[a][i].occupancy();
			shell_check[i] = Orbits[a][i].is_shell();
		}


		string name = Store[a].name;
		double nAtoms = Store[a].nAtoms;

		// We may need to call HartreeFock HF(Latts[a], Orbits[a], Pots[a], Atomic[a], _log);
		// However, calling it is expensive for heavier elements, so we only call it if haven't saved rates. // TODO Seems like this actually isn't expensive...
		// So we need to pass these args through
		Store[a] = Dynamics.SolveAtomicRatesAndPlasmaBEB(max_occ, final_occ, shell_check, !bound_free_exclusions[a],_log);
		Store[a].bound_free_excluded = bound_free_exclusions[a];
		Store[a].name = name;
		Store[a].nAtoms = nAtoms;
		
		// Store[a].R = dropl_R();
		Index[a] = Dynamics.Get_Indexes();
	}
}