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
#ifndef AC4DC_INPUT_CXX_H
#define AC4DC_INPUT_CXX_H

#include "RadialWF.hpp"
#include "Potential.hpp"
#include "Grid.hpp"
#include <vector>
#include <string>
#include <cmath>



struct HFInput
{
	enum class n_pot_t{coulomb, sphere};
	enum class gauge_t{length, velocity};
	enum class ham_mod_t{HF, LDA};
	enum class e_pot_t{V_N, V_Nm1no, V_Nm1};
	
	std::string name = ""; // Name of the atom
	// std::string model;
	// std::string potential = "V_N";
	// std::string me_gauge = "length";
	// std::string hamiltonian = "LDA";
	double omega = 5000;// XFEL field frequency, Ha
	
	int num_time_steps = 0; // Guess number of time steps for time dynamics.
	int omp_threads = 1;
	int Z;

	bool write_charges = false;
	bool write_intensity = false;
	int out_time_steps = 500; // Guess number of time steps for time dynamics.

	double master_tolerance = pow(10, -10);
	double no_exch_tolerance = pow(10, -3);
	double HF_tolerance = pow(10, -6);
	int max_HF_iterations = 500;

	n_pot_t nuclear_potential = n_pot_t::coulomb;
	gauge_t gauge = gauge_t::length;
	ham_mod_t ham_model = ham_mod_t::LDA;
	e_pot_t orbital_potential = e_pot_t::V_N;
	
	// serialisers
	void from_toml(ifstream& ifs);
	void into_toml(ofstream& ifs);
};


// Original input function (deprecated)
class Input
{
public:
	//Read or set the configuration. Assign all the quantum numbers and trial energies, lattice and potential
	Input(char* filename, std::vector<RadialWF> &Orbitals, Grid &Lattice, std::ofstream & log);
	//Input(char* filename, std::vector<RadialWF> &Orbitals,  std::vector<RadialWF> &Virtual, Grid &Lattice, std::ofstream & log);
	Input(const Input & Other);

	// this is a really weird function
	Input& operator=(Input other) {
		if (&other == this) {
			return *this;
		}

		throw "Equals operator called";

		std::swap(name, other.name);
		std::swap(model, other.model);
		std::swap(potential, other.potential);
		std::swap(hamiltonian, other.hamiltonian);
		std::swap(me_gauge, other.me_gauge);
		std::swap(omega, other.omega);
		std::swap(width, other.width);
		std::swap(fluence, other.fluence);
		std::swap(num_time_steps, other.num_time_steps);
		std::swap(omp_threads, other.omp_threads);
		std::swap(Z, other.Z);
		std::swap(write_charges, other.write_charges);
		std::swap(write_intensity, other.write_intensity);
		std::swap(out_time_steps, other.out_time_steps);
		std::swap(Master_tollerance, other.Master_tollerance);
		std::swap(No_exchange_tollerance, other.No_exchange_tollerance);
		std::swap(HF_tollerance, other.HF_tollerance);
		std::swap(max_HF_iterations, other.max_HF_iterations);

		return *this;
	}

	std::string Pot_Model() { return model; }
	std::string Exited_Pot_Model() { return potential; }
	std::string Gauge() { return me_gauge; }
	std::string Name() { return name; }
	double Omega() { return omega; }
	double Width() { return width; }
	double Fluence() { return fluence; }
	void Set_Width(double ext_width) { width = ext_width; }
	void Set_Fluence(double ext_fluence) { fluence = ext_fluence; }

	void Set_Pulse(double ext_omega, double ext_fluence, double ext_width, bool write_ch = false, bool write_int = false, int ext_T_size = 0) {
		omega = ext_omega;
		fluence = ext_fluence;
		width = ext_width;
		num_time_steps = ext_T_size;
	}
	int TimePts() { return num_time_steps; }
	int Hamiltonian();
	int Nuclear_Z() { return Z; }

	int Num_Threads() {return omp_threads; }
	void Set_Num_Threads(int new_num_threads) {omp_threads = new_num_threads; }
	int max_HF_iters() {return max_HF_iterations; }
	double Master_toll() {return Master_tollerance; }
	double No_Exch_toll() {return No_exchange_tollerance; }
	double HF_toll() {return HF_tollerance; }

	bool Write_Charges() {return write_charges; }
	bool Write_Intensity() {return write_intensity; }
	int Out_T_size() {return out_time_steps; }

	~Input();
private:
	std::string name = "";// A name of argv[1] suppilied as filename without extension. Is added to output files.
	std::string model;
	std::string potential = "V_N";
	std::string me_gauge = "length";
	std::string hamiltonian = "LDA";
	double omega = 5000;// XFEL field frequency.
	double width = 5; // XFEL pulse width. Gaussian profile hardcoded.
	double fluence = 0; // XFEL pulse fluence.
	int num_time_steps = 0; // Guess number of time steps for time dynamics.
	int omp_threads = 1;
	int Z;

	bool write_charges = false;
	bool write_intensity = false;
	int out_time_steps = 500; // Guess number of time steps for time dynamics.

	double Master_tollerance = pow(10, -10);
	double No_exchange_tollerance = pow(10, -3);
	double HF_tollerance = pow(10, -6);
	int max_HF_iterations = 500;
};


#endif /* end of include guard: AC4DC_INPUT_CXX_H */
