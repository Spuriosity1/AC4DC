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
#include "Grid.h"
#include "stdafx.h"

Grid::Grid(int num_grid_pts, double r_min, double r_max, double Beta)
{
	double r_tmp = r_min;
	double s_tmp = 0.0;
	double s_i = 0.0;
	unsigned int alert = 1;
	beta = Beta;

	r.clear();
	dr.clear();

	if (num_grid_pts <= 0) { std::cout << "number of radial lattice points is 0 or negative"; }
	else
	{
		r.resize(num_grid_pts);
		dr.resize(num_grid_pts);

		NumPts = num_grid_pts;

		r[0] = r_min;
		dr[0] = r_min;
		s_tmp = r_min + beta*log(r_min);
		ds = (r_max + beta*log(r_max) - s_tmp) / (num_grid_pts - 1);

		for (unsigned int i = 1; i < num_grid_pts; i++)
		{
			s_i = s_tmp + ds;
			alert = 1;
			do
			{
				r_tmp += r_tmp*(s_i - s_tmp) / (r_tmp + beta);
				s_tmp = r_tmp + beta*log(r_tmp);
				alert++;
				if (alert > 100) { break; printf("grid is not created! \n"); }
			} while (fabs(s_i - s_tmp) / ds > 1.0*pow(10.0, -15.0));

			r[i] = r_tmp;
			dr[i] = r_tmp / (r_tmp + beta)*ds;//dr[i] = (dr/ds)_i*ds_i
		}
		r[num_grid_pts - 1] = r_max;
	}
}

Grid::Grid(double r_min, double r_max, double dR_max)
{
	// self-adjustable Grid.
	// Box size r_max.
	// Maximum spacing between the points at infinity dR_max.
	double h = 0.05;
	double exp_h = exp(h);
	double r_tmp = r_min, dr_tmp = r_min*(1+h);
	int n = 0;

	// First part of the grid is fine to account for nuclear potential.
	//      r_n = r_min*exp(h*n) + r_min*n [= r_tmp + r_min*n];
	//     dr_n = h*r_min*exp(h*n) + r_min [= h*r_tmp + r_min];
	while (dr_tmp < dR_max) {
		r.push_back(r_tmp + n*r_min);
		dr.push_back(dr_tmp);
		n++;
		r_tmp = r_tmp*exp_h;
		dr_tmp = h*r_tmp + r_min;
	}
	r_tmp = r.back();
	dr_tmp = dr.back();
	// After an interval has reached dR_max, linear grid.
	while(r_tmp < r_max) {
		r_tmp += dr_tmp;
		r.push_back(r_tmp);
		dr.push_back(dr_tmp);
	}

	NumPts = r.size();
}

Grid::Grid(int num_grid_pts, double r_min, double r_max, std::string mode)
{
	double s = 0;
	double h = 0;
	double tmp = 0;

	r.clear();
	dr.clear();

	if (num_grid_pts <= 0) { std::cout << "number of radial lattice points is 0 or negative"; }
	else
	{
		r.resize(num_grid_pts);
		dr.resize(num_grid_pts);
		r[0] = r_min;
		NumPts = num_grid_pts;

		if (mode == "exponential+")
		{
			if (r_max/r_min > num_grid_pts) {
				s = pow(r_max / r_min - num_grid_pts + 1, 1. / (num_grid_pts - 1));
				h = log(s);
				tmp = s;
				for (int i = 1; i < r.size(); i++)
				{
					r[i] = (tmp + i)*r_min;
					dr[i] = r[i] * h + r_min*(1 - h*i);
					tmp *= s;
				}
			} else mode = "exponential";
		}
		if (mode == "exponential") {
				s = pow(r_max / r_min, 1. / (num_grid_pts - 1));
				h = log(s);
				for (int i = 1; i < r.size(); i++)
				{
					r[i] = r[i - 1] * s;
					dr[i] = r[i] * h;
				}
		}
		if (mode == "linear")
		{
			h = (r_max - r_min) / (num_grid_pts - 1);
			for (int i = 1; i < r.size(); i++)
			{
				r[i] = r[i - 1] + h;
				dr[i] = h;
			}
		}
	}
}

Grid::Grid(vector<double> & X, vector<double> & dX) : r(X), dr(dX)
{
	if (r.size() != dr.size()) NumPts = -1;
	else NumPts = r.size();
	beta = 0;
	ds = 0;
}

void Grid::Extend(double new_max_R)
{
	// Expand an existing mesh by adding more exponentially spaced points to the end.
	double exp_h = r[r.size() - 1]/r[r.size() - 2];
	double h = log(exp_h);
	/*if (h < 0.01) {
		h = 0.01;
		exp_h = exp(0.01);
	}*/
	dr.back() = h * r.back();

	double last_r = exp_h * r.back();
	double last_dr = h * last_r;
	while (r.back() < new_max_R) {
		r.push_back(last_r);
		dr.push_back(last_dr);
		last_r *= exp_h;
		last_dr = h * last_r;
	}

	NumPts = r.size();
}

double Grid::R(int i)
{
	return r[i];
}

double Grid::dR(int i)
{
	return dr[i];
}

double Grid::dS()
{
	return ds;
}

double Grid::dR_dS(int i)
{
	return (r[i] / (beta + r[i]));
}

Grid::~Grid(void)
{}
