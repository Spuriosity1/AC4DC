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
#include <string>
#include <cassert>


namespace Constant
{
	const double Pi = 3.1415926535897932384626433832795;
	const double Alpha = 0.0072973525698;
	const double Alpha2 = Constant::Alpha * Constant::Alpha;
	const double Fm = 1.8897261246 / 100000;

	// Conversions
	// Read these as "One au is 0.024 fs"
	const double fs_per_au = 0.024188843265857;//femtosecond in au
	const double eV_per_Ha = 27.211385;//electron volts in atomic units
	const double J_per_eV = 1.60217662e-19;
	// const double Intensity_in_au = 6.434;// x10^15 W/cm^2 //3.50944758;//x10^2 W/nm^2
	// const double Jcm2_per_Haa02 = 1./6.4230434293;//J/cm^2 
	const double au2_in_barn = 5.2917721067*5.2917721067*1000000;//atomic units to Barns.
	const double au2_in_Mbarn = 5.2917721067*5.2917721067;//atomic units to Mega Barns.
	const double RiemannZeta3 = 1.202056903159594;
	const double Angs_per_au = 0.52917721067; // Bohr radius = 1 atomic unit in Angstrom.
	const double cm_per_au = Angs_per_au*1e-8;
	const double kb_eV = 8.617333262145e-5; // Boltzmann constant, electronvolt per Kelvin
	const double kb_Ha = 8.617333262145e-5/eV_per_Ha; // Boltzmann constant, Ha per Kelvin

	const double Jcm2_per_Haa02 = J_per_eV/cm_per_au/cm_per_au * eV_per_Ha;

	double Wigner3j(double, double, double, double, double, double);
}

namespace PhysicalRate
{
	struct photo //photoionisation rate
	{
		double val;//value of the rate in a.u.
		int hole;//orbital with hole
		double energy; // Ha
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
		std::vector<double> val;
	};

	struct polarize
	{
		std::vector<int> reference;// Contains only occupancies of Orbitals, no Virtual included.
		double refEnergy;// Energy of the reference configuration.
		/* E1-selected excited configurations.
		excited[][0] - orbital in reference from which to excite.
		excited[][1] - orbital in reference to which to excite.
		*/
		std::vector<std::vector<int>> excited;// Includes both Orbitals and Virtual.
		std::vector<double> extEnergy;
		std::vector<double> ImpactCrossSections; // Reduced transition dipole matrix elements form 'reference' to 'excited'.
		int index = 0;
	};
}

// todo check if these are needed
// static const double Moulton_5[5] = { 251. / 720., 646. / 720., -264. / 720., 106. / 720., -19. / 720. }; //Adams-Moulton method
// static const double Bashforth_5[5] = { 1901. / 720., -1378. / 360., 109. / 30., -637. / 360., 251. / 720. }; //Adams-Bashforth method
