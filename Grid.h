#pragma once

#include <vector>
#include <string>

using namespace std;

class Grid
{	/* The Grid class creates the coordinate grid on the interval from [r_min, r_max]
	with "num_grid_pts" points. The grid has not constant node spacing. It is constant
	for the following variable
	S = R + beta*ln(R)
	According to Vladimir Dzuba the best value for beth=4.
	*/
public:
	Grid(int num_grid_pts, double r_min, double r_max, double Beta = 4);
	Grid(double r_max, double dR_max, int Z);//same linear logarithm, but with maximum dR_max. Defines number of points itself
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

	int size() { return NumPts; }

private:
	std::vector<double> r, dr;
	double ds;
	double beta;
	int NumPts;
};
