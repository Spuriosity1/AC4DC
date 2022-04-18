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

#include "MolInp.hpp"
#include "Constant.hpp"
#include <fstream>
#include <sstream>
#include <iostream>
#include <stdexcept>
#include <map>
#include "HartreeFock.hpp"
#include "ComputeRateParam.hpp"


MolInp::MolInp(const char* filename, std::ofstream & log)
{
	// Input file for molecular ionization calculation.
	std::map<std::string, std::vector<std::string>> FileContent;

	std::cout<<"Opening molecular file "<<filename<<"... ";
	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	std::ifstream infile(filename);
	if (infile.good())
        std::cout<<"Success!"<<"\n";
    else {
        std::cerr<<"Failed."<<"\n";
		exit(EXIT_FAILURE); // chuck a hissy fit and quit.
        return;
    }
	std::string comment = "//";
	std::string curr_key = "";

	while (!infile.eof())
	{
		std::string line;
		getline(infile, line);
		if (!line.compare(0, 2, comment)) continue;
		if (!line.compare(0, 1, "")) continue;
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = std::vector<std::string>(0);
			}
			curr_key = line;
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	size_t num_atoms = FileContent["#ATOMS"].size();

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
		std::stringstream stream(FileContent["#VOLUME"][n]);

		if (n == 0) stream >> unit_V;
		if (n == 1) stream >> loss_geometry.L0;
		if (n == 2) stream >> loss_geometry;
	}

	for (size_t n = 0; n < FileContent["#OUTPUT"].size(); n++) {
		std::stringstream stream(FileContent["#OUTPUT"][n]);
		char tmp;

		if (n == 0) stream >> out_T_size;
		if (n == 1) stream >> out_F_size;
		if (n == 2) {
			stream >> tmp;
			if (tmp == 'Y') write_charges = true;
		}
		if (n == 3) {
			stream >> tmp;
			if (tmp == 'Y') write_intensity = true;
		}
		if (n == 4) {
			stream >> tmp;
			if (tmp != 'Y') write_md_data = false;
		}
	}

	for (size_t n = 0; n < FileContent["#PULSE"].size(); n++) {
		std::stringstream stream(FileContent["#PULSE"][n]);

		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
		if (n == 3) stream >> pulse_shape;
	}

	for (size_t n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		std::stringstream stream(FileContent["#NUMERICAL"][n]);

		if (n == 0) stream >> num_time_steps;
		if (n == 1) stream >> omp_threads;
		if (n == 2) stream >> min_elec_e;
		if (n == 3) stream >> max_elec_e;
		if (n == 4) stream >> num_elec_points;
		if (n == 5) stream >> elec_grid_type;
		if (n == 6) stream >> elec_grid_type.num_low;
		if (n == 7) stream >> elec_grid_type.transition_e;
		if (n == 8) stream >> elec_grid_type.min_coulomb_density;

	}

	// for (size_t n = 0; n < FileContent["#GRID"].size(); n++) {
	// 	std::stringstream stream(FileContent["#GRID"][n]);
	// 	if (n == 0) stream >> min_elec_e;
	// 	if (n == 1) stream >> max_elec_e;
	// 	if (n == 2) stream >> num_elec_points;
	// 	if (n == 3) stream >> elec_grid_type;
	// 	if (n == 4) stream >> elec_grid_type.num_low;
	// 	if (n == 5) stream >> elec_grid_type.transition_e;
	// }

	// Hardcode the boundary conditions
	elec_grid_type.zero_degree_inf = 3;
	elec_grid_type.zero_degree_0 = 0;

	const std::string bc = "\033[33m"; // begin colour escape code
	const std::string clr = "\033[0m"; // clear escape code
	// 80 equals signs
	const std::string banner = "================================================================================";
	std::cout<<banner<<"\n";
	std::cout<<bc<<"Unit cell size: "<<clr<<unit_V<<" A^3"<<"\n";
	std::cout<<bc<<"Droplet L0:     "<<clr<<loss_geometry.L0<<" A"<<"\n";
	std::cout<<bc<<"Droplet Shape:  "<<clr<<loss_geometry<<"\n"<<"\n";

	std::cout<<bc<<"Photon energy:  "<<clr<<omega<<" eV"<<"\n";
	std::cout<<bc<<"Pulse fluence:  "<<clr<<fluence*10000<<" J/cm^2 = "<<10000*fluence/omega/Constant::J_per_eV<<"ph cm^-2"<<"\n";
	std::cout<<bc<<"Pulse FWHM:     "<<clr<<width<<" fs"<<"\n";
	std::cout<<bc<<"Pulse shape:    "<<clr<<pulse_shape<<"\n"<<"\n";

	std::cout<<bc<<"Electron grid:  "<<clr<<min_elec_e<<" ... "<<max_elec_e<<" eV"<<"\n";
	std::cout<<    "                "<<num_elec_points<<" points"<<"\n";
	std::cout<<bc<<"Grid type:      "<<clr<<elec_grid_type<<"\n";
	std::cout<<bc<<"Low energy cutoff for Coulomb logarithm estimation: "<<clr<<elec_grid_type.transition_e<<"eV"<<"\n";
	std::cout<<bc<<"Minimum num electrons per unit cell for Coulomb logarithm to be considered: "<<clr<<elec_grid_type.min_coulomb_density<<"\n";
	std::cout<<"\n";

	std::cout<<bc<<"ODE Iteration:  "<<clr<<num_time_steps<<" timesteps"<<"\n"<<"\n";

	std::cout<<bc<<"Output:         "<<clr<<out_T_size<<" time grid points"<<"\n";
	std::cout<<    "                "<<out_F_size<<" energy grid points"<<"\n";
	std::cout<<banner<<"\n";

	// Convert to number of photon flux.
	omega /= Constant::eV_per_Ha;
	fluence *= 10000/Constant::Jcm2_per_Haa02/omega;

	// Convert to atomic units.
	width /= Constant::fs_per_au;
	loss_geometry.L0 /= Constant::Angs_per_au;
	unit_V /= Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;

	elec_grid_type.min_coulomb_density /= unit_V;

	min_elec_e /= Constant::eV_per_Ha;
	max_elec_e /= Constant::eV_per_Ha;

	elec_grid_type.transition_e /= Constant::eV_per_Ha;

	// Reads the very top of the file, expecting input of the form
	// H 2
	// O 1
	// Then scans for atomic files of the form 
	// input/H.inp
	// input/O.inp
	// Store is then populated with the atomic data read in below.
	for (size_t i = 0; i < num_atoms; i++) {
		std::string at_name;
		double at_num;

		std::stringstream stream(FileContent["#ATOMS"][i]);
		stream >> at_name >> at_num;

		Store[i].nAtoms = at_num/unit_V;
		Store[i].name = at_name;
		// Store[i].R = radius;

		at_name = "input/" + at_name + ".inp";

		Atomic.push_back(Input((char*)at_name.c_str(), Orbits[i], Latts[i], log));
		// Overrides pulses found in .inp files
		Atomic.back().Set_Pulse(omega, fluence, width);
		Atomic.back().Set_Num_Threads(omp_threads);

		Potential U(&Latts[i], Atomic[i].Nuclear_Z(), Atomic[i].Pot_Model());

		Pots[i] = U;
	}

	if (!validate_inputs()) {
		std::cerr<<"\n"<<"\n"<<"\n"<<"Exiting..."<<"\n";
		throw std::runtime_error(".mol input file is invalid");
	}
}

bool MolInp::validate_inputs() {
	bool is_valid=true;
	std::cerr<<"\033[31;1m";
	if (omega <= 0 ) { std::cerr<<"ERROR: pulse omega must be positive"; is_valid=false; }
	if (width <= 0 ) { std::cerr<<"ERROR: pulse width must be positive"; is_valid=false; }
	if (fluence <= 0 ) { std::cerr<<"ERROR: pulse fluence must be positive"; is_valid=false; }
	if (num_time_steps <= 0 ) { std::cerr<<"ERROR: got negative number of timesteps"; is_valid=false; }
	if (out_T_size <= 0) { std::cerr<<"ERROR: system set to output zero timesteps"; is_valid=false; }
	if (out_F_size <= 0) { std::cerr<<"ERROR: system set to output zero energy grid points"; is_valid=false; }
	if (loss_geometry.L0 <= 0) { std::cerr<<"ERROR: radius must be positive"; is_valid=false; }
	if (omp_threads <= 0) { omp_threads = 4; std::cerr<<"Defaulting number of OMP threads to 4"; }

	if (elec_grid_type.mode == GridSpacing::unknown) {
		std::cerr<<"ERROR: Grid spacing not recognised - must start with (l)inear, (q)uadratic,";
		std::cerr<<" (e)xponential, (h)ybrid or (p)owerlaw"; is_valid=false;
	}
	if (elec_grid_type.mode == GridSpacing::hybrid || elec_grid_type.mode == GridSpacing::powerlaw) {
		if (elec_grid_type.num_low <= 0 || elec_grid_type.num_low >= num_elec_points) { 
			std::cerr<<"Defaulting number of dense points to "<<num_elec_points/2;
			elec_grid_type.num_low = num_elec_points/2;
		}
	}
	if (elec_grid_type.transition_e <= min_elec_e || elec_grid_type.transition_e >= max_elec_e) {
		std::cerr<<"Defaulting low-energy cutoff to "<<max_elec_e/4;
		elec_grid_type.transition_e = max_elec_e/4;
	}
	

	// unit cell volume.
	if (unit_V <= 0) { std::cerr<<"ERROR: unit xell volume must be positive"; is_valid=false; }

	// Electron grid style
	if(min_elec_e < 0 || max_elec_e < 0 || max_elec_e <= min_elec_e) { std::cerr<<"ERROR: Electron grid specification invalid"; is_valid=false; }
	if (num_time_steps <= 0 ) { std::cerr<<"ERROR: got negative number of energy steps"; is_valid=false; }
	std::cerr<<"\033[0m";
	return is_valid;
}



// TODO: Delete this spaghetti code
void MolInp::calc_rates(std::ofstream &_log, bool recalc) {
	// Loop through atomic species.
	for (size_t a = 0; a < Atomic.size(); a++) {
		// HartreeFock HF(Latts[a], Orbits[a], Pots[a], Atomic[a], _log);

		// This Computes the parameters for the rate equations to use, loading them into Init.
		ComputeRateParam Dynamics(Latts[a], Orbits[a], Pots[a], Atomic[a], recalc);
		std::vector<int> final_occ(Orbits[a].size(), 0);
		std::vector<int> max_occ(Orbits[a].size(), 0);
		for (size_t i = 0; i < max_occ.size(); i++) {
			if (fabs(Orbits[a][i].Energy) > Omega()) final_occ[i] = Orbits[a][i].occupancy();
			max_occ[i] = Orbits[a][i].occupancy();
		}

		std::string name = Store[a].name;
		double nAtoms = Store[a].nAtoms;

		Store[a] = Dynamics.SolvePlasmaBEB(max_occ, final_occ, _log);
		Store[a].name = name;
		Store[a].nAtoms = nAtoms;
		// Store[a].R = dropl_R();
		Index[a] = Dynamics.Get_Indexes();
	}
}
