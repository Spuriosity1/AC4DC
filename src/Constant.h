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

/* This header contains some Physical constants,
functions for calculations of Wigner 3j, 6j symbols, and Clebsh-Gordan coefficients,
and some data containers used througout the code. */
#include <vector>
using namespace std;

namespace Constant
{
	const double Pi = 3.1415926535897932384626433832795;
	const double Alpha = 0.0072973525698;
	const double Alpha2 = Constant::Alpha * Constant::Alpha;
	const double Fm = 1.8897261246 / 100000;
	const double fs_in_au = 0.02418884;//femtosecond in au
	const double eV_in_au = 27.211385;//electron volts in atomic units
	const double Intensity_in_au = 6.434;// x10^15 W/cm^2 //3.50944758;//x10^2 W/nm^2
	const double Fluence_in_au = 0.155689291;//J/cm^2
	const double au2_in_barn = 5.2917721067*5.2917721067*1000000;//atomic units to Barns.
	const double au2_in_Mbarn = 5.2917721067*5.2917721067;//atomic units to Mega Barns.
	const double RiemannZeta3 = 1.202056903159594;
	const double au_in_Angs = 0.52917721067; // Bohr radius = 1 atomic unit in Angstrom.

	double Wigner3j(double, double, double, double, double, double);
}

namespace CustomDataType
{
	struct photo//photoionization rate
	{
		double val;//value of the rate in a.u.
		int hole;//orbital with hole
		double energy;
	};

	struct fluor//fluorescence rate
	{
		double val;//value of the rate in a.u.
		int fill;//filled orbital
		int hole;//orbital with hole
	};

	struct auger//auger decay rate
	{
		double val;
		int hole;
		int fill;
		int eject;
		double energy = 0;
	};

	struct ffactor//form factor for Q_mesh values in FormFactor class
	{
		int index;
		vector<double> val;
	};

	struct polarize
	{
		vector<int> reference;// Contains only occupancies of Orbitals, no Virtual included.
		double refEnergy;// Energy of the reference configuration.
		/* E1-selected excited configurations.
		excited[][0] - orbital in reference from which to excite.
		excited[][1] - orbital in reference to which to excite.
		*/
		vector<vector<int>> excited;// Includes both Orbitals and Virtual.
		vector<double> extEnergy;
		vector<double> Dipoles; // Reduced transition dipole matrix elements form 'reference' to 'excited'.
		int index = 0;
	};

	struct EIIdata
	{
		int init;
		vector<int> fin;
		vector<int> occ;
		vector<float> ionB;
		vector<float> kin;
	};
}

static const double Moulton_5[5] = { 251. / 720., 646. / 720., -264. / 720., 106. / 720., -19. / 720. }; //Adams-Moulton method
static const double Bashforth_5[5] = { 1901. / 720., -1378. / 360., 109. / 30., -637. / 360., 251. / 720. }; //Adams-Bashforth method
static const double gaussX_10[10] = {-0.973906528517, -0.865063366689, -0.679409568299, -0.433395394129, -0.148874338982,
									  0.148874338982, 0.433395394129, 0.679409568299, 0.865063366689, 0.973906528517};
static const double gaussW_10[10] = {0.066671344308688, 0.14945134915058, 0.21908636251598, 0.26926671931000, 0.29552422471475,
									 0.29552422471475, 0.26926671931000, 0.21908636251598, 0.14945134915058, 0.066671344308688};
static const double gaussX_13[13] = {-0.9841830547185881, -0.9175983992229779, -0.8015780907333099, -0.6423493394403401,
								-0.4484927510364467, -0.23045831595513477, 0., 0.23045831595513477, 0.4484927510364467,
								0.6423493394403401, 0.8015780907333099, 0.9175983992229779, 0.9841830547185881};
static const double gaussW_13[13] = {0.04048400476531614, 0.09212149983772834, 0.13887351021978736, 0.17814598076194582,
								0.20781604753688845, 0.2262831802628971, 0.23255155323087365, 0.2262831802628971, 0.20781604753688845,
								0.17814598076194582, 0.13887351021978736, 0.09212149983772834, 0.04048400476531614};
