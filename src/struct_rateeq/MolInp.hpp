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

#ifndef AC4DC_CXX_MOLINP_H
#define AC4DC_CXX_MOLINP_H

#include "Input.hpp"
#include "Constant.hpp"
#include "SplineBasis.hpp"
#include "GridSpacing.hpp"
#include "LossGeometry.hpp"
#include "Pulse.hpp"

class MolInp
{
	// Molecular input for coupled atom/electron plasma calcualtions.
public:
	MolInp(const char* filename, std::ofstream & log);
	~MolInp() {}

	
	std::vector<Input> Atomic; // Vector of atomic input objects
	std::vector<RateData::Atom> Store; // Stores all atomic parameters: EII, photoionisation, fluorescence, 

	std::vector<Potential> Pots;
	std::vector<std::vector<RadialWF>> Orbits;
	std::vector<Grid> Latts;
	std::vector<std::vector<std::vector<int>>> Index;

	double Omega() {return omega;}
	double Width() {return width;}
	double Fluence() {return fluence;}
	int Num_Time_Steps() {return num_time_steps;}
	double dropl_R() {return radius;}

  	// void Set_Fluence(double new_fluence) {fluence = new_fluence;}
	bool Write_Charges() {return write_charges; }
	bool Write_Intensity() {return write_intensity; }
	bool Write_MD_data() {return write_md_data; }

	int Out_T_size() {return out_T_size; }
	int Out_F_size() {return out_F_size; }

	double Max_Elec_E() {return max_elec_e;}
	double Min_Elec_E() {return min_elec_e;}
	size_t Num_Elec_Points() {return num_elec_points;}


	std::string name = "";

	// Scans for atomic process rates in the folder output/[atom name]/Xsections and recalculates if absent.
	// 'Recalculate' flag skips scanning and forces recomputation.
	// The result is available as Store.
    void calc_rates(std::ofstream &_log, bool recalc=true);

	GridSpacing elec_grid_type;
	LossGeometry loss_geometry;
	PulseShape pulse_shape = PulseShape::none;

	int num_time_steps = -1; // Number of time steps for time dynamics

protected:

	bool validate_inputs();

	double omega = -1;// XFEL photon energy, au.
	double width = -1; // XFEL pulse width in au. Gaussian profile hardcoded.
	double fluence = -1; // XFEL pulse fluence, au.
	
	int out_T_size = 0; // Unlike atomic input, causes to output all points.
	int out_F_size = 100; // Number of F grid points to use when outputting.
	double radius = -1; // Droplet radius
	int omp_threads = 1;

	// Flags for outputting
	bool write_charges = false;
	bool write_intensity = false;
	bool write_md_data = true;

	// unit volume.
	double unit_V = -1.;

	// Electron grid style
	double min_elec_e = -1;
	double max_elec_e = -1;
	size_t num_elec_points = -1; // Number of cells in the free-electron distribution expansion

};


#endif
