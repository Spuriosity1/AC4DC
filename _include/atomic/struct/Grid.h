/**
 * @file Grid.h
 * @brief Defines the Grid class
 * @details The Grid class creates the coordinate grid on the interval from [r_min, r_max]
	with "num_grid_pts" points. The grid's spacing follows a logarithmic relationship. See Grid.
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

using namespace std;

/** 
 * @brief The Grid class creates the coordinate grid on the interval from [r_min, r_max]
	with "num_grid_pts" points. The grid has not constant node spacing. It is constant
	for the following variable:
	S = R + beta*ln(R).
	According to Vladimir Dzuba the best value for beta=4.
*/
class Grid
{	
private:
	vector<double> r, dr;
	double ds;
	double beta;
	int NumPts;
public:
	Grid(int num_grid_pts, double r_min, double r_max, double Beta = 4);
	Grid(double r_min, double r_max, double dR_max);//same linear logarithm, but with maximum dR_max. 
	// Defines number of points on its own. Great for integrals with oscillating functions.
	Grid(int X) { NumPts = X; }//empty lattice to be defined elsewhere
	Grid(vector<double> & X, vector<double> & dX);

	// Exponential grid for integrals over Gaussian basis set and uniform for continuum states
	Grid(int num_grid_pts, double r_min, double r_max, std::string mode);
	//		Grid(const std::string& filename);
	~Grid(void);

	void Extend(double new_max_R);
	double R(int i);
	double dR(int i);
	double dR_dS(int i);
	double dS();

	const Grid& operator=(const Grid& lattice)
	{
		r = lattice.r;
		dr = lattice.dr;
		ds = lattice.ds;
		NumPts = lattice.NumPts;
		beta = lattice.beta;

		return *this;
	}

	int size() {
		return NumPts;
	}

  const double * ptr_R() {return r.data();}
  const double * ptr_dR() {return dr.data();}
	const decltype(r) & all_R() const {return r;}
};
