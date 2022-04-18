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
#include <fstream>
#include "HFInput.h"


//Class calculates total photoionization crossection for a given orbital

class DecayRates
{
public:
	DecayRates(Grid &Lattice, std::vector<RadialWF> &Orbitals, Potential &U, HFInput & Input);

	std::vector<photo> Photo_Ion(double omega, std::ofstream & log); // All photoinonization crossections. Position in vector indicates orbital
	std::vector<fluor> Fluor(); // Fluorescence rates for all channels.
	std::vector<auger> Auger(std::vector<int> Max_occ, std::ofstream & log); // Auger decay rates for all channels.

	std::vector<double> FT_density(double Q_min = 0, double Q_max = 2, int Q_size = 20);
	~DecayRates();

private:
	int IntegrateContinuum(Grid &Lattice, Potential &U, std::vector<RadialWF> &Core, RadialWF* Current, int c = 0);
	Grid& lattice;
	std::vector<RadialWF>& orbitals;
	Potential& u;
	HFInput & input;

	// EII internal intepolated data.
	Grid CntLattice = Grid(0);// Grid for continuum wave calculations.
	Potential CntU;
	std::vector<RadialWF> CntOrbitals;// Interpolated onto continuum grid "orbitals".
};
