/**
 * @file Dipole.h 
 * @brief 
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

#ifndef AC4DC_CXX_DIPOLE_H
#define AC4DC_CXX_DIPOLE_H


namespace Dipole
{
	// p - impactor electron momentum
	// p_s - secondary electron momentum
	// B - ionization potential
	// occ - orbitals occupancy
	double DsigmaBEB(double T, double W, double B, double u, int occ);
	double sigmaBEB(double T, double B, double u, int occ);
	// Int_0^(T - B) dW W d(sigmaBED)/dW - total energy absorbance
	double sigmaBEBw1(double T, double B, double u, int occ);
}


#endif /* end of include guard: AC4DC_CXX_DIPOLE_H */
