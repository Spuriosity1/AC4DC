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
#include "HFInput.hpp"
#include "Constant.hpp"

#include <toml++/toml.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>
#include <map>
#include <limits>

const std::map<char, unsigned> Orbital::chemists_l = {
	{'s', 0},
	{'p', 1},
	{'d', 2},
	{'f', 3}
};

const bool Orbital::is_chemist_l(char symbol){
	return Orbital::chemists_l.find(symbol) != Orbital::chemists_l.end();
}

// weak string matching
bool string_similar(const std::string& s1, const std::string& s2, 
		const size_t num_match = std::numeric_limits<size_t>::max())
{
	size_t num_compared = std::min(std::min(s1.size(), s2.size()), num_match);
	if (num_compared == 0) {
		return false;
	}
	for (size_t i=0; i<std::min(s1.size(), num_match); i++){
		if (std::tolower(s1[i]) != std::tolower(s2[i])){
			return false;
		}
	}
	return true;
}

// void throw_if_ndef(toml::node_view t, const char* key){
// 	if (!t[key].has_value()){
// 		std::cerr<<"No configuration found for "<<key<<"\n";
// 		throw std::runtime_error("Bad config file");
// 	}
// }

void HFInput::interpret_npot(const std::string& s){
	if ( string_similar(s, "coulomb")){
		this->nuclear_potential = n_pot_t::coulomb;
	} else if ( string_similar(s, "sphere")){
		this->nuclear_potential = n_pot_t::sphere;
	} else {
		std::cerr << "Bad specifier: "<< s <<std::endl;
		throw std::runtime_error("Bad specifier");
	}
}


void HFInput::interpret_gauge(const std::string& s){
	if ( string_similar(s, "length") ){
		this->gauge = gauge_t::length;
	} else if ( string_similar(s, "velocity") ){
		this->gauge = gauge_t::velocity;
	} else {
		std::cerr << "Bad specifier: "<< s <<std::endl;
		throw std::runtime_error("Bad specifier");
	}
}

void HFInput::interpret_hmodel(const std::string& s){
	if ( string_similar(s, "HF", 1) ){
		this->ham_model = ham_mod_t::HF;
	} else if ( string_similar(s, "LDA", 1) ){
		this->ham_model = ham_mod_t::LDA;
	} else {
		std::cerr << "Bad specifier: "<< s <<std::endl;
		throw std::runtime_error("Bad specifier");
	}
}

void HFInput::interpret_epot(const std::string& s){
	if ( string_similar(s, "V_N") ){
		this->orbital_potential = e_pot_t::V_N;
	} else if ( string_similar(s, "V_N-1no") ){
		this->orbital_potential = e_pot_t::V_Nm1no;
	} else if ( string_similar(s, "V_N-1") ){
		this->orbital_potential = e_pot_t::V_Nm1;
	} else {
		std::cerr << "Bad specifier: "<< s <<std::endl;
		throw std::runtime_error("Bad specifier");
	}
}

void HFInput::read_orbital(const std::string& entry){
	std::stringstream ss(entry);
	Orbital orb;
	char ang_label;
	ss >> orb.n >> ang_label >> orb.max_occ;
	if (orb.n < 1 || !Orbital::is_chemist_l(ang_label) || orb.max_occ < 0){
		std::cerr<<"Bad orbital configuration\n";
		throw std::runtime_error("Bad specifier");
	} 
	orb.l = Orbital::chemists_l.at(ang_label);
	orbitals.push_back(orb);
}

/**
 * @brief Parses the stream 'ifs', expected to be TOML-formatted
 * 
 * @param ifs 
 */
void HFInput::from_toml(std::ifstream& ifs){
	toml::table config = toml::parse(ifs);
	this->name = *config["name"].value<std::string>();
	
	auto gs = config["electronic_gs"];

	this->Z = *gs["nuclear_charge"].value<int64_t>();
	
	// cursed flag interpreters

	// throw_if_ndef(gs, "nuclear_potential");
	interpret_npot(*gs["nuclear_potential"].value<std::string>());
	// throw_if_ndef(gs, "hamiltonian_form");
	interpret_hmodel(*gs["hamiltonian_form"].value<std::string>());
	// throw_if_ndef(gs, "potential_model");
	interpret_epot(*gs["potential_model"].value<std::string>());
	// throw_if_ndef(gs, "photon_gauge");
	interpret_gauge(*gs["photon_gauge"].value<std::string>());

	
	toml::array& orbs = *gs["electron_config"].as_array();

	orbitals.resize(0);
	for (auto&& o : orbs){
		read_orbital(*o.value<std::string>());
	}

	auto grid = config["radialgrid"].as_array();

	auto tol = config["tolerance"].as_array();

}

void HFInput::into_toml(std::ofstream& ifs){
	throw "Not Implemented";
}





/*

//// DEPRECEATED

Input::Input(char *filename, std::vector<RadialWF> &Orbitals, Grid &Lattice, std::ofstream & log)
{
	log << "[ Atomic ] Input file: " << filename << "\n";

	name = filename;
	size_t lastdot = name.find_last_of(".");
	if (lastdot != std::string::npos) name = name.substr(0, lastdot);
	size_t lastslash = name.find_last_of("/");
	if (lastdot != std::string::npos) name = name.substr(lastslash+1);

	std::cout << "[ Atomic ] Opening atomic file "<< filename << "... ";
	std::ifstream infile(filename);
	if (infile.good())
        std::cout<<"Success!\n";
    else {
        std::cerr<<"\033[31;1mFailed!\033[0m\n";
		throw std::runtime_error("Could not find atomic input file."); // chuck a hissy fit and quit.
        return;
    }

	std::map<std::string, std::vector<std::string>> FileContent;
	std::string comment = "//";
	std::string curr_key;

	while (!infile.eof() && infile.is_open()) {
		std::string line;
		getline(infile, line);

		if (!line.compare(0, 2, comment)) continue;
		if (!line.compare(0, 2, "")) continue;
		if (!line.compare(0, 1, "#")) {
			if ( FileContent.find(line) == FileContent.end() ) {
				FileContent[line] = std::vector<std::string>(0);
			}
			curr_key = line;
		} else {
			FileContent[curr_key].push_back(line);
		}
	}

	for (int n = 0; n < FileContent["#PULSE"].size(); n++) {
		std::stringstream stream(FileContent["#PULSE"][n]);
		if (n == 0) stream >> omega;
		if (n == 1) stream >> width;
		if (n == 2) stream >> fluence;
	}

	for (int n = 0; n < FileContent["#OUTPUT"].size(); n++) {
		std::stringstream stream(FileContent["#OUTPUT"][n]);
		if (n == 0) stream >> out_time_steps;
		if (n == 1) {
			char tmp;
			stream >> tmp;
			if (tmp == 'Y') write_charges = true;
		}
		if (n == 2) {
			char tmp;
			stream >> tmp;
			if (tmp == 'Y') write_intensity = true;
		}
	}

	int line_indentifier = 0, num_grid_pts, num_orbitals, N, L, current_orbital = 0, occupancy;
	double r_min, r_box;
	std::map<char, int> angular = {{'s', 0}, {'p', 1}, {'d', 2}, {'f', 3}};

	for (int n = 0; n < FileContent["#NUMERICAL"].size(); n++) {
		std::stringstream stream(FileContent["#NUMERICAL"][n]);
		if (n == 0) stream >> num_grid_pts;
		if (n == 1) stream >> r_min;
		if (n == 2) stream >> r_box;
		if (n == 3) stream >> num_time_steps;
		if (n == 4) stream >> omp_threads;
		if (n == 5) {
			int tmp = 0;
			stream >> tmp;
			HF_tollerance = pow(10, tmp);
		}
		if (n == 6) {
			int tmp = 0;
			stream >> tmp;
			Master_tollerance = pow(10, tmp);
		}
		if (n == 7) {
			int tmp = 0;
			stream >> tmp;
			No_exchange_tollerance = pow(10, tmp);
		}
		if (n == 8) stream >> max_HF_iterations;
	}

	// Assign a default value to avoid undefiend comparisons
	num_orbitals = -10;
	for (int n = 0; n < FileContent["#ATOM"].size(); n++) {
		std::stringstream stream(FileContent["#ATOM"][n]);
		if (n == 0) stream >> Z;
		if (n == 1) stream >> model;
		if (n == 2) stream >> hamiltonian;
		if (n == 3) stream >> num_orbitals;
		// TODO: Refactor this to be less of a kludge
		if (n > 3 && n < num_orbitals + 4) {
			char tmp;
			stream >> N >> tmp >> occupancy;
			if (occupancy == 0){
				std::cerr << "[ Atomic ] \033[31;1mOrbital with N=" << N << ", L=" << angular[tmp] << " has occupancy=0 \033[0m\n";
			}
			if (N == 0){
				std::cerr << "[ Atomic ] \033[31;1m Orbital with N=0 encountered: "<<filename<<"\033[0m\n";
				std::cerr << "[ Atomic ] Did you specify the right number of orbitals?" << std::endl;
				throw std::runtime_error("Bad atomic input");
			}
			Orbitals.push_back(RadialWF(num_grid_pts));
			Orbitals.back().set_N(N);
			Orbitals.back().set_L(angular[tmp]);//setting L overwrites occupancy with 4L+2
			Orbitals.back().set_occupancy(occupancy);
		}
		if (n == num_orbitals + 4) stream >> potential;
		if (n == num_orbitals + 5) stream >> me_gauge;
	}

	Grid lattice(num_grid_pts, r_min/Z, r_box, 4);
	Lattice = lattice;
	omega /= Constant::eV_per_Ha;
	fluence *= 10000/omega/Constant::Jcm2_per_Haa02;
}

Input::Input(const Input & Other)
{
	name = Other.name;
	model = Other.model;
	potential = Other.potential;
	me_gauge = Other.me_gauge;
	hamiltonian = Other.hamiltonian;
	omega = Other.omega;
	width = Other.width;
	fluence = Other.fluence;
	num_time_steps = Other.num_time_steps;
	omp_threads = Other.omp_threads;
	Z = Other.Z;
	write_charges = Other.write_charges;
	write_intensity = Other.write_intensity;
	out_time_steps = Other.out_time_steps;
	Master_tollerance = Other.Master_tollerance;
	No_exchange_tollerance = Other.No_exchange_tollerance;
	HF_tollerance = Other.HF_tollerance;
	max_HF_iterations = Other.max_HF_iterations;
}

int Input::Hamiltonian()
{
	if (hamiltonian == "HF") return 0;
	else return 1;
}

Input::~Input()
{
}


*/