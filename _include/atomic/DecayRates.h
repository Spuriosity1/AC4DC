/**
 * @file DecayRates.h
 * 
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
#pragma once

#include "RadialWF.h"
#include "Grid.h"
#include "Potential.h"
#include <vector>
#include "Constant.h"
#include "RateData.h"
#include <fstream>
#include "HFInputParam.h"

//Class calculates total photoionization crossection for a given orbital
class DecayRates
{
public:
	DecayRates(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, InputData::HFInputParam & Input);

	static void set_Max_occ(vector<int> & Max_occ, vector<RadialWF> orbitals);

	
	vector<RateData::photo> Photo_Ion(double omega, ofstream & log); // All photoinonization crossections. Position in vector indicates orbital
	vector<RateData::fluor> Fluor(); // Fluorescence rates for all channels.
	vector<RateData::auger> Auger(vector<int> Max_occ, ofstream & log); // Auger decay rates for all channels.

	vector<double> FT_density(double Q_min = 0, double Q_max = 2, int Q_size = 20);
	~DecayRates();

private:
	int IntegrateContinuum(const Grid &Lattice,  Potential &U, const vector<RadialWF> &Core, RadialWF& Current, int c = 0);
	Grid& lattice;
	vector<RadialWF>& orbitals;
	Potential& u;
	InputData::HFInputParam & input;

	// EII internal intepolated data.
	// Grid CntLattice = Grid(0);// Grid for continuum wave calculations.
	//Potential CntU;
	// vector<RadialWF> CntOrbitals;// Interpolated onto continuum grid "orbitals".
};
