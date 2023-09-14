// /*===========================================================================
// This file is part of AC4DC.

//     AC4DC is free software: you can redistribute it and/or modify
//     it under the terms of the GNU General Public License as published by
//     the Free Software Foundation, either version 3 of the License, or
//     (at your option) any later version.

//     AC4DC is distributed in the hope that it will be useful,
//     but WITHOUT ANY WARRANTY; without even the implied warranty of
//     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//     GNU General Public License for more details.

//     You should have received a copy of the GNU General Public License
//     along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
// ===========================================================================*/
// #include "HFInputParam.h"
// #include "Constant.h"
// #include <fstream>
// #include <sstream>
// #include <iostream>
// #include <map>

// HFInputParam::Input(char *filename, vector<RadialWF> &Orbitals, Grid &Lattice, ofstream & log)
// {
// 	log << "[ Atomic ] Input file: " << filename << endl;

// 	name = filename;
// 	size_t lastdot = name.find_last_of(".");
// 	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
// 	size_t lastslash = name.find_last_of("/");
// 	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

// 	cout << "[ Atomic ] Opening atomic file "<< filename << "... ";
// 	ifstream infile(filename);
// 	if (infile.good()){
//         std::cout<<"Success!"<<endl;
// 	} else {
//         std::cerr<<"\033[31;1mFailed!\033[0m"<<endl;
// 		throw runtime_error("Could not find atomic input file."); // chuck a hissy fit and quit.
//         return;
//     }

// 	map<string, vector<string>> FileContent;
// 	string comment = "//";
// 	string curr_key;

// 	// Put the information into FileContent
// 	while (!infile.eof() && infile.is_open()) {
// 		string line;
// 		getline(infile, line);

// 		if (!line.compare(0, 2, comment)) continue;
// 		if (!line.compare(0, 2, "")) continue;
// 		if (!line.compare(0, 1, "#")) {
// 			if ( FileContent.find(line) == FileContent.end() ) {
// 				FileContent[line] = vector<string>(0);
// 			}
// 			curr_key = line;
// 		} else {
// 			FileContent[curr_key].push_back(line);
// 		}
// 	}
// 	// Seems like this is old code that hasn't been removed, presumably this part was moved to MolInp.cpp. Though the #OUTPUT section is still present in the files. TODO try turning this off. - S.P.
// 	// Update: 
// 	/////////////////////////////////////////////////////
// 	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
// 		stringstream stream(FileContent["#PULSE"][n]);
// 		if (n == 0) stream >> omega;
// 		if (n == 1) stream >> width;
// 		if (n == 2) stream >> fluence;
// 	}

// 	for (int n = 0; n < FileContent["#OUTPUT"].size(); n++) {
// 		stringstream stream(FileContent["#OUTPUT"][n]);
// 		if (n == 0) stream >> out_time_steps;
// 		if (n == 1) {
// 			char tmp;
// 			stream >> tmp;
// 			if (tmp == 'Y') write_charges = true;
// 		}
// 		if (n == 2) {
// 			char tmp;
// 			stream >> tmp;
// 			if (tmp == 'Y') write_intensity = true;
// 		}
// 	}
// 	/////////////////////////////////////////////////////////

// 	// This part almost definitely has old code that hasn't been removed, some of the parameters seems to pertain to plasma hyperparameters that the .mol files already do, so may be redundant
// 	/////////////////////////////////////////////////////////
// 	int line_indentifier = 0, num_grid_pts, num_orbitals, N, L, current_orbital = 0, occupancy;
// 	double r_min, r_box;
// 	map<char, int> subshell_to_angular = {{'s', 0}, {'p', 1}, {'d', 2}, {'f', 3},{'N',-10}};   // 'N' flags usage of the shell energy instead of the orbital.

// 	for (int n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
// 		stringstream stream(FileContent["#NUMERICAL"][n]);
// 		if (n == 0) stream >> num_grid_pts;
// 		if (n == 1) stream >> r_min;
// 		if (n == 2) stream >> r_box;
// 		if (n == 3) stream >> num_time_steps;
// 		if (n == 4) stream >> omp_threads;  // Overridden by .mol
// 		if (n == 5) {
// 			int tmp = 0;
// 			stream >> tmp;
// 			HF_tollerance = pow(10, tmp);
// 		}
// 		if (n == 6) {
// 			int tmp = 0;
// 			stream >> tmp;
// 			Master_tollerance = pow(10, tmp);
// 		}
// 		if (n == 7) {
// 			int tmp = 0;
// 			stream >> tmp;
// 			No_exchange_tollerance = pow(10, tmp);
// 		}
// 		if (n == 8) stream >> max_HF_iterations;
// 	}
// 	////////////////////////////////////////////////////////

// 	// Assign a default value to avoid undefiend comparisons
// 	num_orbitals = -10;
// 	for (int n = 0; n < FileContent["#ATOM"].size(); n++) {
// 		stringstream stream(FileContent["#ATOM"][n]);
// 		if (n == 0) stream >> Z;
// 		if (n == 1) stream >> model;
// 		if (n == 2) stream >> hamiltonian;
// 		if (n == 3) stream >> num_orbitals;
// 		//( Orbitals )\\		// TODO: Refactor this to be less of a kludge
// 		if (n > 3 && n < num_orbitals + 4) {
// 			char tmp;  // the subshell (spdf)
// 			stream >> N >> tmp >> occupancy;
// 			if (occupancy == 0){
// 				cerr << "[ Atomic ] \033[31;1mOrbital with N=" << N << ", L=" << subshell_to_angular[tmp] << " has occupancy=0 \033[0m" << endl;
// 			}
// 			// It seems like this effectively throws if the wrong number of orbitals is specified. As virtual orbitals are allowed. -S.P.
// 			if (N == 0){
// 				cerr << "[ Atomic ] \033[31;1m Incorrect number of orbital types specified in input file: "<<filename<<"\033[0m" << endl;
// 				cerr << "Note that a shell approximation for orbitals only counts for 1" << endl;
// 				throw runtime_error("Bad atomic input");
// 			}
// 			Orbitals.push_back(RadialWF(num_grid_pts)); // Adds a radial wavefunction
// 			Orbitals.back().set_N(N);
// 			Orbitals.back().set_L(subshell_to_angular[tmp]);//setting L overwrites occupancy with 4L+2 (which is simpler for DecayRates.cpp), so occupancy must be set after this.
// 			Orbitals.back().set_occupancy(occupancy);
// 			Orbitals.back().flag_shell();  // remember if this is a shell
// 		}
// 		if (n == num_orbitals + 4) stream >> potential;
// 		if (n == num_orbitals + 5) stream >> me_gauge;
// 	}

// 	Grid lattice(num_grid_pts, r_min/Z, r_box, 4);
// 	Lattice = lattice;
// 	omega /= Constant::eV_per_Ha;
// 	fluence *= 10000/omega/Constant::Jcm2_per_Haa02;
// }

// HFInputParam::HFInputParam(const HFInputParam & Other)
// {
// 	name = Other.name;
// 	model = Other.charge_model;
// 	potential = Other.potential;
// 	me_gauge = Other.me_gauge;
// 	hamiltonian = Other.hamiltonian;
// 	omega = Other.omega;
// 	width = Other.width;
// 	fluence = Other.fluence;
// 	num_time_steps = Other.num_time_steps;
// 	omp_threads = Other.omp_threads;
// 	Z = Other.Z;
// 	write_charges = Other.write_charges;
// 	write_intensity = Other.write_intensity;
// 	out_time_steps = Other.out_time_steps;
// 	Master_tollerance = Other.Master_tollerance;
// 	No_exchange_tollerance = Other.No_exchange_tollerance;
// 	HF_tollerance = Other.HF_tollerance;
// 	max_HF_iterations = Other.max_HF_iterations;
// }

// int HFInputParam::Hamiltonian()
// {
// 	if (hamiltonian == "HF") return 0;
// 	else return 1;
// }

// HFInputParam::~HFInputParam()
// {
// }
//
//
#include "HFInputParam.h"
const json HFInputParam::schema = R"(
		{
			"$schema": "http://json-schema.org/draft-07/schema#",
			"title": "HFInputParam",
			"properties": 
			{
				"name": {"type": "string"},
				"nuclear_Z": {"type": "integer", "minimum": 1},
				"orbitals": {"type": "string", "uniqueValues": "true"},
				"charge_model": {"type": "string", "oneOf": ["coulomb", "sphericalWell"]},
				"hamiltonian": {"type": "string", "oneOf": ["HF", "LDA"]},
				"potential": {"type": "string", "oneOf": ["V_N","V_Nm1","V_Nm1_no"] },
				"me_gauge": {"type":"string", "oneOf": ["length", "velocity"] },
				"num_grid_points": {"type": "integer", "minimum": 1},
				"grid_min": {"type": "number", "minimum": 0},
				"grid_max": {"type": "number", "minimum": 0},
				"omp_threads": {"type": "integer", "minimum": 1},
				"Master_tollerance": {"type": "number", "minimum": 0},
				"No_exchange_tollerance": {"type": "number", "minimum": 0},
				"HF_tollerance": {"type": "number", "minimum": 0},
				"max_HF_iterations": {"type": "integer", "minimum": 1},
			},
			"additionalProperties": false,
			"minProperties": 15
		}

	)"_json;

