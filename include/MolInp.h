/**
 * @file MolInp.h
 * @brief Molecular input for coupled atom/electron plasma calculations.
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

#ifndef AC4DC_CXX_MOLINP_H
#define AC4DC_CXX_MOLINP_H

#include "Input.h"
#include "Constant.h"
#include "SplineBasis.h"
#include "GridSpacing.hpp"
#include "LossGeometry.hpp"
#include "Pulse.h"


class MolInp
{
	// Molecular input for coupled atom/electron plasma calculations.
public:
	MolInp(const char* filename, ofstream & _log);
	~MolInp() {}

	/// Vector of atomic input objects
	vector<Input> Atomic; 
	/// Stores all atomic parameters: EII, photoionisation, fluorescence, 
	vector<RateData::Atom> Store; 

	vector<Potential> Pots;
	vector<vector<RadialWF>> Orbits;
	vector<Grid> Latts;
	vector<vector<vector<int>>> Index;

	double Omega() {return omega;}
	double Width() {return width;}
	double Cutoff_Inputted() {return cutoff_flag;}
	double Simulation_Cutoff() {return simulation_cutoff_time;}
	double Fluence() {return fluence;}
	int Num_Time_Steps() {return num_time_steps;}
	double dropl_R() {return radius;}

	int Out_T_size() {return out_T_size; }
	int Out_F_size() {return out_F_size; }

	int Plasma_Threads(){return omp_threads;}

	string Filtration_File(){return filtration_file;}

	string Load_Folder(){return load_folder;}
	double Load_Time_Max(){return simulation_resume_time_max/Constant::fs_per_au;}
	bool Using_Input_Timestep(){return loading_uses_input_timestep;}

	double Grid_Update_Period(){return grid_update_period;}

	string name = "";

	// Scans for atomic process rates in the folder output/[atom name]/Xsections and recalculates if absent.
	// 'Recalculate' flag skips scanning for primary ionization data and forces computation for the specific photon energy used. As it only generates data for a single photon energy, it is quicker than the initial generation of data for when recalc is False. However, the data is not (currently) stored for use in later runs.
	// The result is available as Store.
    void calc_rates(ofstream &_log, bool recalc=true);

	GridSpacing elec_grid_type;
	DynamicGridPreset elec_grid_preset;
	Cutoffs param_cutoffs;
	ManualGridBoundaries elec_grid_regions;
	LossGeometry loss_geometry;
	PulseShape pulse_shape = PulseShape::none;
	// Gaussian modifications
	double timespan_factor = 0;  // if set, is the total simulated timespan in units of the pulse width (fwhm) 
	double negative_timespan_factor = 0; // if set, the total simulated timespan of the pulse before the t=0 point.

	int num_time_steps = -1; // Number of time steps for time dynamics

	double time_update_gap = 0; // interval **in fs** between each cout'd time.	May be useful to set to a high number for HPC or when otherwise not using ncurses.
    int steps_per_live_plot_update = 20; // Interval  **in steps** between plotting of the free electron distribution to _live_plot.png. Setting to 1 (updating every step) has a negligible effect on speed outside of very fast high step count simulations. 

	double electron_source_fraction = 0;
	double electron_source_energy = -1;
	double electron_source_duration = 1; // As fraction of entire pulse
	char electron_source_type = 'c'; // (c)onstant: rate is ([intensity]/[initial intensity]) * [initial photoion. rate of atoms in target]  *  [electron source fraction]. | (p)roportional: rate is [source fraction] * [total photionisation rate of all atoms in target].   

protected:

	bool validate_inputs();

	double omega = -1;// XFEL photon energy, au.
	double width = -1; // XFEL pulse width in au. Gaussian profile hardcoded.

	bool use_fluence = false;
	bool use_count = false;
	bool use_intensity = false; // peak intensity
	double fluence = -1; // XFEL pulse fluence, au.
	double photon_count = -1;
	double peak_intensity = -1;
	
	// Simulation end time, inputted as fs. Note that simulation (w/ rate equations) temporal width is 4*FWHM.
	double simulation_cutoff_time = 1; 
	// Whether the cutoff was given
	bool cutoff_flag = false;

	int out_T_size = 0; // Unlike atomic input, causes to output all points.
	int out_F_size = 100; // Number of F grid points to use when outputting.
	double radius = -1; // Droplet radius
	int omp_threads = 1;

	// unit volume.
	double unit_V = -1.;

	// Simulation loading parameters
	double simulation_resume_time_max = INFINITY; // will attempt to load closest to this time but not after.
	string load_folder = ""; // If "" don't load anything. 
	string filtration_file = ""; // If "" don't load anything.
	bool loading_uses_input_timestep = false; // 

	// Dynamic grid
	double grid_update_period; // time period between dynamic grid updates, fs.

	// Rate calc exclusions
	std::vector<bool> bound_free_exclusions;


};


#endif
