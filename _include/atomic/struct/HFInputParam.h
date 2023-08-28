/**
 * @file Input.h
 * @brief  Defines the Input class, which holds the inputs from .inp atomic files. (Don't invite Input.h to your party, they're a wet blanket).
 * @details The Input class reads and stores the configuration-defining information from the input files or an Input object.  
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
#ifndef AC4DC_INPUT_CXX_H
#define AC4DC_INPUT_CXX_H

#include "RadialWF.h"
#include "Potential.h"
#include "Grid.h"
#include <vector>
#include <string>
#include <cmath>

using namespace std;

class Input
{
public:
	//Read or set the configuration. Assign all the quantum numbers and trial energies, lattice and potential
	Input(char* filename, vector<RadialWF> &Orbitals, Grid &Lattice, ofstream & log);
	//Input(char* filename, vector<RadialWF> &Orbitals,  vector<RadialWF> &Virtual, Grid &Lattice, ofstream & log);
	Input(const Input & Other);

	Input& operator=(Input other) {
		if (&other == this) {
			return *this;
		}

		swap(name, other.name);
		swap(model, other.model);
		swap(potential, other.potential);
		swap(hamiltonian, other.hamiltonian);
		swap(me_gauge, other.me_gauge);
		swap(omega, other.omega);
		swap(width, other.width);
		swap(fluence, other.fluence);
		swap(num_time_steps, other.num_time_steps);
		swap(omp_threads, other.omp_threads);
		swap(Z, other.Z);
		swap(write_charges, other.write_charges);
		swap(write_intensity, other.write_intensity);
		swap(out_time_steps, other.out_time_steps);
		swap(Master_tollerance, other.Master_tollerance);
		swap(No_exchange_tollerance, other.No_exchange_tollerance);
		swap(HF_tollerance, other.HF_tollerance);
		swap(max_HF_iterations, other.max_HF_iterations);

		return *this;
	}

	string Pot_Model() { return model; }
	string Exited_Pot_Model() { return potential; }
	string Gauge() { return me_gauge; }
	string Name() { return name; }
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
	string name = "";// A name of argv[1] suppilied as filename without extension. Is added to output files.
	string model;
	string potential = "V_N";
	string me_gauge = "length";
	string hamiltonian = "LDA";
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
