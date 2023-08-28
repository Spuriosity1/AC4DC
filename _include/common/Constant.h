/**
 * @file Constant.h
 * @brief This header contains some Physical constants,
 * @details note
functions for calculations of Wigner 3j, 6j symbols, and Clebsh-Gordan coefficients,
and some data containers used throughout the code. 
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

#include <vector>
#include <string>
#include <cassert>
using namespace std;

namespace Constant
{
	const double Pi = 3.1415926535897932384626433832795;
	const double Alpha = 0.0072973525698;
	const double Alpha2 = Constant::Alpha * Constant::Alpha;
	const double Fm = 1.8897261246 / 100000;

	//**// Conversions
	// Read these as "One au is 0.024 fs".
	// Note unit conversion direction is consistent, i.e. always divide to convert to atomic units.
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
	//**//
	const double kb_eV = 8.617333262145e-5; // Boltzmann constant, electronvolt per Kelvin
	const double kb_Ha = 8.617333262145e-5/eV_per_Ha; // Boltzmann constant, Ha per Kelvin

	const double Jcm2_per_Haa02 = J_per_eV/cm_per_au/cm_per_au * eV_per_Ha;

}

typedef std::vector<double> bound_t; // TODO I'm debating removing this since there are lots of std::vector<double> declarations that this makes confusing -S.P.
// even though it is just vector, the typedef is useful in case we want to change the implementation of bound_t later for speed - A.S.

/// Ensures we have consistent removal of decimal places which is important for loading sims.
/// Necessary due to high precision of Constant::fs_per_au
static const double loading_t_precision = 9;