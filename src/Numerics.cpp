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
#include "Numerics.h"
#include <vector>
#include "EigenSolver.h"
#include <algorithm>

static const double adams_10[10] = { 2082753.0 / 7257600.0, 9449717.0 / 7257600.0, -11271304.0 / 7257600.0, 16002320.0 / 7257600.0, -17283646.0 / 7257600.0,
13510082.0 / 7257600.0, -7394032.0 / 7257600.0, 2687864.0 / 7257600.0, -583435.0 / 7257600.0, 57281.0 / 7257600.0 };

static const double adams_5[5] = { 251. / 720., 646. / 720., -264. / 720., 106. / 720., -19. / 720. };

static const double adams_4[4] = { 1. / 24., -5. / 24., 19. / 24., 9. / 24. };

static const double Lagrange[5][10] =
{
	{ -7129. / 2520.,           9.,         -18.,        28., -63. / 2., 126. / 5.,      -14.,  36. / 7.,  -9. / 8.,    1. / 9. },
    {       -1. / 9., -481. / 280.,           4.,  -14. / 3.,  14. / 3.,  -7. / 2., 28. / 15.,  -2. / 3.,   1. / 7.,  -1. / 72. },
    {       1. / 72.,     -1. / 4., -153. / 140.,    7. / 3.,  -7. / 4.,   7. / 6., -7. / 12.,   1. / 5., -1. / 24.,  1. / 252. },
    {     -1. / 252.,     3. / 56.,     -3. / 7., -37. / 60.,   3. / 2.,  -3. / 4.,   1. / 3., -3. / 28., 3. / 140., -1. / 504. },
    {      1. / 504.,    -1. / 42.,      1. / 7.,   -2. / 3.,  -1. / 5.,        1.,  -1. / 3.,  2. / 21., -1. / 56.,  1. / 630. }
};

static const double Lagrange_4[2][4] =
{
	{-11./6.,     3., -3./2.,  1./3.},
	{ -1./3., -1./2.,     1.,  2./3.}
};

static const double gaussX_2[2] = {-0.57735026918963, 0.57735026918963};
static const double gaussW_2[2] = {1., 1.};

static const double gaussX_4[4] = {-0.8611363115941, -0.3399810435849, 0.3399810435849, 0.8611363115941};
static const double gaussW_4[4] = {0.34785484513745, 0.65214515486255, 0.65214515486255, 0.34785484513745};

static const double gaussX_5[5] = {-0.9061798459387, -0.5384693101057, 0., 0.5384693101057, 0.9061798459387};
static const double gaussW_5[5] = {0.23692688505619, 0.47862867049937, 0.56888888888889, 0.47862867049937, 0.23692688505619};

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

Adams::Adams(Grid &Latt, int AdamsOrder) : Lattice(Latt)
{
	Adams_Coeff.clear();

	A.clear();
	B.clear();
	C.clear();
	D.clear();
	X.clear();
	Y.clear();

	A.resize(Lattice.size());
	B.resize(Lattice.size());
	C.resize(Lattice.size());
	D.resize(Lattice.size());
	X.resize(Lattice.size());
	Y.resize(Lattice.size());

	if (AdamsOrder < 6)	{ Adams_N = 5; }
	else { Adams_N = 10; }
//	Adams_N = 10;

	Adams_Coeff.resize(Adams_N);

	if (Adams_N == 5)
	{
		for (int i = 0; i < Adams_N; i++)
		{
			Adams_Coeff[i] = adams_5[i];
		}
	}
	else
	{
		for (int i = 0; i < Adams_N; i++)
		{
			Adams_Coeff[i] = adams_10[i];
		}
	}
}

Adams::~Adams()
{}

void Adams::Integrate(RadialWF* Psi, int start_pt, int end_pt)
{
	int incr, start, end;
	double Det, F_tmp, G_tmp;
	FirstMaxima = 0;
	NumNodes = 0;
	vector<double> dF_dR(Lattice.size(), 0), dG_dR(Lattice.size(), 0);

	if (start_pt > end_pt)
	{
		incr = -1;
		start = start_pt - Adams_N;
		end = end_pt;
	}
	else
	{
		incr = 1;
		start = start_pt + Adams_N;
		end = end_pt ;
	}

	for (int j = start; incr*j >= start - incr*Adams_N; j-=incr)
	{
		dF_dR[j] = A[j] * Psi->F[j] + B[j] * Psi->G[j] + X[j];
		dG_dR[j] = C[j] * Psi->F[j] + D[j] * Psi->G[j] + Y[j];
	}

	for (int i = start; incr*i <= incr*end; i+=incr)
	{
		Det = (1.0 - incr * Adams_Coeff[0] * Lattice.dR(i) * A[i])*(1.0 - incr * Adams_Coeff[0] * Lattice.dR(i) * D[i])
			- Adams_Coeff[0] * Lattice.dR(i) * C[i] * Adams_Coeff[0] * Lattice.dR(i) * B[i];

		F_tmp = Psi->F[i - incr] + incr*Adams_Coeff[0] * Lattice.dR(i) * X[i] ;
		G_tmp = Psi->G[i - incr] + incr*Adams_Coeff[0] * Lattice.dR(i) * Y[i] ;
		if (std::isnan(F_tmp) || std::isnan(G_tmp)){throw std::runtime_error("F_tmp or G_tmp is nan!");}
		for (int j = 1; j < Adams_N; j++)
		{
			F_tmp += incr * Adams_Coeff[j] * dF_dR[i - incr*j] * Lattice.dR(i - incr*j);
			G_tmp += incr * Adams_Coeff[j] * dG_dR[i - incr*j] * Lattice.dR(i - incr*j);
		}
		if (std::isnan(Psi->F[i])){

			throw std::runtime_error("Psi->F[i] is nan!0");

		}
		
		Psi->F[i] = (F_tmp * (1.0 - incr * Lattice.dR(i) * Adams_Coeff[0] * D[i]) + G_tmp * incr * Lattice.dR(i) * Adams_Coeff[0] * B[i]) / Det;
		Psi->G[i] = (F_tmp * incr * Lattice.dR(i) * Adams_Coeff[0] * C[i] + G_tmp * (1.0 - incr * Lattice.dR(i) * Adams_Coeff[0] * A[i])) / Det;
		if (std::isnan(Psi->F[i])){	
			
			throw std::runtime_error("Psi->F[i] is nan!1");

		}
		
		dF_dR[i] = A[i] * Psi->F[i] + B[i] * Psi->G[i] + X[i];
		dG_dR[i] = C[i] * Psi->F[i] + D[i] * Psi->G[i] + Y[i];
		if (std::isnan(Psi->F[i])){throw std::runtime_error("Psi->F[i] is nan!2");}
		
		if (dF_dR[i-incr] * dF_dR[i] < 0)
		{
			if (FirstMaxima == 0)
			{
				FirstMaxima = abs(i);
			}
			else
			{
				NumNodes++;
				if (fabs(Psi->F[FirstMaxima]) < fabs(Psi->F[i])) { FirstMaxima = abs(i); }
			}
		}
		if (std::isnan(Psi->F[i])){throw std::runtime_error("Psi->F[i] is nan!3 End of func: Adams::Integrate()");}
	}
}

void Adams::Integrate(std::vector<double> &Func, std::vector<double> &Result, int start_pt, int end_pt)
{
	int incr, start, end;
	double Func_tmp=0.0;

	if (start_pt > end_pt)
	{
		incr = -1;
		start = start_pt - Adams_N;
		end = end_pt;
	}
	else
	{
		incr = 1;
		start = start_pt + Adams_N;
		end = end_pt;
	}

	//check if the first Adams_N points are required to be calculated or they are given by the user.

	for (int i = 1; i < Adams_N; i++)
	{
		Func_tmp += fabs(Result[start_pt + incr*i]);
	}

	if (Func_tmp == 0)
	{
		std::vector< std::vector<double>> LeftMatr;
		std::vector<double> RightVect;

		int Lagrange_N = 9;//Always calculating first 3 points. If Adams_N < 10 it will overwrite an extra points

		LeftMatr.resize(Lagrange_N);
		for (int i = 0; i < Lagrange_N; i++)
		{
			LeftMatr[i].resize(Lagrange_N);
		}

		RightVect.resize(Lagrange_N);

		// Set up the LeftMatr. Simpler than start
		for (int i = 0; i < Lagrange_N; i++)
		{
			for (int j = 0; j < Lagrange_N; j++)
			{
				if (i < Lagrange_N / 2) { LeftMatr[i][j] = Lagrange[i + 1][j + 1]; }
				else { LeftMatr[i][j] = -Lagrange[Lagrange_N - i - 1][Lagrange_N - j - 1]; }
			}
		}
		//set up RightVect
		for (int i = 0; i < Lagrange_N; i++)
		{
			if (i < Lagrange_N / 2)
			{
				RightVect[i] = -Lagrange[i + 1][0] * Result[start_pt];
			}
			else
			{
				RightVect[i] = Lagrange[Lagrange_N - i - 1][Lagrange_N] * Result[start_pt];
			}

			RightVect[i] += Func[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
		}

		EigenSolver W;
		W.SolveSystem(LeftMatr, RightVect, Lagrange_N);

		for (int i = 0; i < Lagrange_N; i++)
		{
			Result[start_pt + incr*(i + 1)] = RightVect[i];
		}
	}

//	Main loop with Adams_N order
	for (int i = start; incr*i <= incr*end; i += incr)
	{
		Func_tmp = Result[i - incr] + incr * Adams_Coeff[0] * Lattice.dR(i) * Func[i];

		for (int j = 1; j < Adams_N; j++)
		{
			Func_tmp += incr * Adams_Coeff[j] * Func[i - incr*j] * Lattice.dR(i - incr*j);
		}

		Result[i] = Func_tmp;
	}
}

double Adams::Integrate(std::vector<double>* Func, int start_pt, int end_pt)
{
	/*
	vector<double> Result = vector<double>(Func->size(), 0.);
	Integrate(*Func, Result, start_pt, end_pt);
	return Result.back();
	*/
	int incr, start, end;
	vector<double> Result(Adams_N);

	if (start_pt > end_pt)
	{
		incr = -1;
		start = start_pt - Adams_N;
		end = end_pt;
	}
	else
	{
		incr = 1;
		start = start_pt + Adams_N;
		end = end_pt;
	}

	//check if the first Adams_N points are required to be calculated or they are given by the user.
	vector<vector<double>> LeftMatr;
	vector<double> RightVect;

	int Lagrange_N = 10;//Always calculating first 10 points. If Adams_N < 10 it will overwrite an extra points

	LeftMatr.resize(Lagrange_N);
	for (int i = 0; i < Lagrange_N; i++)
	{
		LeftMatr[i].resize(Lagrange_N);
	}

	RightVect.resize(Lagrange_N);

	// Set up the LeftMatr. Simpler than start
	for (int i = 0; i < Lagrange_N; i++)
	{
		for (int j = 0; j < Lagrange_N; j++)
		{
			if (i < Lagrange_N / 2) { LeftMatr[i][j] = Lagrange[i + 1][j + 1]; }
			else { LeftMatr[i][j] = -Lagrange[Lagrange_N - i - 1][Lagrange_N - j - 1]; }
		}
	}
	//set up RightVect
	for (int i = 0; i < Lagrange_N; i++)
	{
		if (i < Lagrange_N / 2)
		{
			RightVect[i] = -Lagrange[i + 1][0] * Result[0];
		}
		else
		{
			RightVect[i] = Lagrange[Lagrange_N - i - 1][Lagrange_N] * Result[0];
		}

		RightVect[i] += Func->at(start_pt + incr*(i + 1)) * Lattice.dR(start_pt + incr*(i + 1));  // Func->at(i) equiv. to density[i]
	}

	EigenSolver W;
	W.SolveSystem(LeftMatr, RightVect, Lagrange_N);

	assert(Result.size() <= RightVect.size());

	for (int i = 0; i < Result.size(); i++)
	{
		if (std::isnan(Result[i])){throw std::runtime_error("Adams::Integrate returned nan!");}
		Result[i] = RightVect[i];
	}

	for (int i = start; incr*i <= incr*end; i += incr)
	{
		if (std::isnan(Result[Adams_N - 1])){throw std::runtime_error("Adams::Integrate returned nan!");}
		Result[Adams_N-1] = Result[Adams_N - 2] + incr * Adams_Coeff[0] * Lattice.dR(i) * Func->at(i);

		for (int j = 1; j < Adams_N; j++)
		{
			Result[Adams_N-1] += incr * Adams_Coeff[j] * Func->at(i - incr*j) * Lattice.dR(i - incr*j);
		}

		for (int j = 0; j < (Adams_N-1); j++)
		{
			Result[j] = Result[j+1];
		}
	}
	return Result[Adams_N - 1];
}

void Adams::StartAdams(RadialWF* Psi, int start_pt, bool forward)
{
	std::vector< std::vector<double>> LeftMatr;
	std::vector<double> RightVect;

	int Lagrange_N = 9;//Always calculating first 10 points. If Adams_N < 10 it will overwrite an extra points

	int incr;
	if (forward) { incr = 1; }
	else { incr = -1; }

	LeftMatr.resize(2 * Lagrange_N);
	for (int i = 0; i < 2 * Lagrange_N; i++)
	{
		LeftMatr[i].resize(2 * Lagrange_N);
	}

	RightVect.resize(2 * Lagrange_N);

	//set up upper left
	for (int i = 0; i < Lagrange_N; i++)
	{
		for (int j = 0; j < Lagrange_N; j++)
		{
			if (i < Lagrange_N / 2) { LeftMatr[i][j] = Lagrange[i + 1][j + 1]; }
			else { LeftMatr[i][j] = - Lagrange[Lagrange_N - i - 1][Lagrange_N - j - 1]; }

			if (i == j)
			{
				LeftMatr[i][j] -= incr*A[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			}
		}
	}

	//set up lower right
	for (int i = 0; i < Lagrange_N; i++)
	{
		for (int j = 0; j < Lagrange_N; j++)
		{
			if (i < Lagrange_N / 2) { LeftMatr[Lagrange_N + i][Lagrange_N + j] = Lagrange[i + 1][j + 1]; }
			else { LeftMatr[Lagrange_N + i][Lagrange_N + j] = - Lagrange[Lagrange_N - i - 1][Lagrange_N - j - 1]; }

			if (i == j)
			{
				LeftMatr[Lagrange_N + i][Lagrange_N + j] -= incr*D[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			}
		}
	}

	//set up upper right
	for (int i = 0; i < Lagrange_N; i++)
	{
		for (int j = 0; j < Lagrange_N; j++)
		{
			LeftMatr[i][Lagrange_N + j] = 0.0;
			if (i == j)
			{
				LeftMatr[i][Lagrange_N + j] -= incr*B[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			}
		}
	}

	//set up lower left
	for (int i = 0; i < Lagrange_N; i++)
	{
		for (int j = 0; j < Lagrange_N; j++)
		{
			LeftMatr[Lagrange_N + i][j] = 0.0;
			if (i == j)
			{
				LeftMatr[Lagrange_N + i][j] -= incr*C[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			}
		}
	}

	//set up RightVect
	for (int i = 0; i < Lagrange_N; i++)
	{
		if (i < Lagrange_N / 2)
		{
			RightVect[i] = -Lagrange[i + 1][0] * Psi->F[start_pt] + incr*X[start_pt + incr*(i+1)] * Lattice.dR(start_pt + incr*(i+1));
			RightVect[i + Lagrange_N] = -Lagrange[i + 1][0] * Psi->G[start_pt] + incr*Y[start_pt + incr*(i+1)] * Lattice.dR(start_pt + incr*(i+1));
		}
		else
		{
			RightVect[i] = Lagrange[Lagrange_N - i - 1][Lagrange_N] * Psi->F[start_pt] + incr*X[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			RightVect[i + Lagrange_N] = Lagrange[Lagrange_N - i - 1][Lagrange_N] * Psi->G[start_pt] + incr*Y[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
		}
		if (std::isinf(RightVect[i])){throw std::runtime_error("RightVect element is inf!");}
		if (std::isnan(RightVect[i])){throw std::runtime_error("RightVect element is nan!");}
		if (std::isnan(RightVect[i+Lagrange_N])){throw std::runtime_error("RightVect element is nan!");}
	}
		
	//Solve the system of inhomogenious equations
	// LaftMatr * X = RightVect
	// the result is returned to RightVect

	EigenSolver W;
	W.SolveSystem(LeftMatr, RightVect, (2 * Lagrange_N));

	for (int i = 0; i < Lagrange_N; i++)
	{		
		Psi->F[start_pt + incr*(i + 1)] = RightVect[i];
		Psi->G[start_pt + incr*(i + 1)] = RightVect[Lagrange_N + i];
		if (std::isinf(Psi->F[start_pt + incr*(i + 1)])){throw std::runtime_error("Psi->F[] is inf!");}
		if (std::isnan(Psi->F[start_pt + incr*(i + 1)])){throw std::runtime_error("Psi->F[] is nan!");}
		if (std::isnan(Psi->G[start_pt + incr*(i + 1)])){throw std::runtime_error("Psi->G[] is nan!");}
		if (std::isinf(Psi->G[start_pt + incr*(i + 1)])){throw std::runtime_error("Psi->G[] is inf!");}
	}
}

std::vector<double> Adams::GreenOrigin(RadialWF* Psi)
{
	//this function returns the following integral
	//int_(0)^(Lattice.R(end_pt)) Psi.F*Y*Lattice.dR,
	// where Y includes the Wronskian.

	std::vector<double> Result;
	double Func_tmp;

	Result.clear();
	Result.resize( Lattice.size() );

	Result[0] = 0.5*Lattice.dR(0) * Psi->F[0] * Y[0];//trapezoid rule for first interval [0...Lattice.dR(0)]

	for (int i = 1; i < Adams_N; i++)
	{
		Result[i] = Result[i - 1] + 0.5*Lattice.dR(i) * Psi->F[i] * Y[i] +
					0.5*Lattice.dR(i - 1) * Psi->F[i - 1] * Y[i - 1];
	}

	for (int i = Adams_N; i <= Psi->pract_infinity(); i++)
	{
		Func_tmp = Result[i - 1] + Adams_Coeff[0] * Lattice.dR(i) * Psi->F[i]*Y[i];

		for (int j = 1; j < Adams_N; j++)
		{
			Func_tmp += Adams_Coeff[j] * Psi->F[i - j]*Y[i - j] * Lattice.dR(i - j);
		}

		Result[i] = Func_tmp;
	}

	return Result;
}

std::vector<double> Adams::GreenInfinity(RadialWF* Psi)
{
	//this function returns the following integral
	//int_(Lattice.R(start_pt))^(Lattice.R(end_pt)) Psi.F*Y*Lattice.dR,
	// where Y is inhomogenous term in the original ODE. 1/W multiplier is ignored, since it is constant and resultant function will be normalized. It should be constant, but just in case it is recalculated in each grid point

	std::vector<double> Result(Lattice.size(), 0.);
	double Func_tmp;
	int Infty = Psi->pract_infinity();

	Result[Infty] = 0.5*Lattice.dR(Infty) * Psi->F[Infty] * Y[Infty];

	for (int i = Infty - 1; i > (Infty - Adams_N); i--)
	{
		Result[i] = Result[i + 1] + 0.5*Lattice.dR(i) * Psi->F[i] * Y[i] +
			0.5*Lattice.dR(i + 1) * Psi->F[i + 1] * Y[i + 1];
	}

	for (int i = (Infty - Adams_N); i >= 0; i--)
	{
		Func_tmp = Result[i + 1] + Adams_Coeff[0] * Lattice.dR(i) * Psi->F[i] * Y[i];

		for (int j = 1; j < Adams_N; j++)
		{
			Func_tmp += Adams_Coeff[j] * Psi->F[i + j] * Y[i + j] * Lattice.dR(i + j);
		}

		Result[i] = Func_tmp;
	}

	return Result;
}

void Adams::Integrate_ODE(std::vector<double> &f, int start_pt, int end_pt)
{
	//solves equation df/dr = A*f + X;
	std::vector<double> df_dR(f.size(), 0.);
	std::vector< std::vector<double>> LeftMatr;
	std::vector<double> RightVect;

	int Lagrange_N = 9;//Always calculating first 10 points. If Adams_N < 10 it will overwrite an extra points
	int incr, start, end;
	double Det=0;

	if (start_pt > end_pt)
	{
		incr = -1;
		start = start_pt - Adams_N;
		end = end_pt;
	}
	else
	{
		incr = 1;
		start = start_pt + Adams_N;
		end = end_pt ;
	}


	//Check if first Lagrange_N points have being precalculated
	for (int i = 1; i <= Lagrange_N; i++)
	{
		Det += fabs(f[i]);
	}

	if (Det == 0)//Only first point is given. Calculate remaining.
	{
		LeftMatr.resize(Lagrange_N);
		for (int i = 0; i < Lagrange_N; i++)
		{
			LeftMatr[i].resize(Lagrange_N);
		}

		RightVect.resize(Lagrange_N);

		//set up upper left
		for (int i = 0; i < Lagrange_N; i++)
		{
			for (int j = 0; j < Lagrange_N; j++)
			{
				if (i < Lagrange_N / 2) { LeftMatr[i][j] = Lagrange[i + 1][j + 1]; }
				else { LeftMatr[i][j] = -Lagrange[Lagrange_N - i - 1][Lagrange_N - j - 1]; }

				if (i == j)
				{
					LeftMatr[i][j] -= incr*A[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
				}
			}
		}


		//set up RightVect
		for (int i = 0; i < Lagrange_N; i++)
		{
			if (i < Lagrange_N / 2)
			{
				RightVect[i] = -Lagrange[i + 1][0] * f[start_pt] + incr*X[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			}
			else
			{
				RightVect[i] = Lagrange[Lagrange_N - i - 1][Lagrange_N] * f[start_pt] + incr*X[start_pt + incr*(i + 1)] * Lattice.dR(start_pt + incr*(i + 1));
			}
		}

		//Solve simple inhomogenious equation
		// LaftMatr * X = RightVect
		// the result is returned to RightVect

		EigenSolver W;
		W.SolveSystem(LeftMatr, RightVect, Lagrange_N);

		for (int i = 0; i < Lagrange_N; i++)
		{
			f[start_pt + incr*(i + 1)] = RightVect[i];
		}
	}

	//Calculate derivatives in first Lagrange_N+1 points
	for (int i = 0; i <= Lagrange_N; i++)
	{
		df_dR[start_pt + incr*i] = A[start_pt + incr*i] * f[start_pt + incr*i] + X[start_pt + incr*i];
	}


	for (int i = start; incr*i <= incr*end; i += incr)
	{
		Det = (1.0 - incr * Adams_Coeff[0] * Lattice.dR(i) * A[i]);

		f[i] = f[i - incr] + incr*Adams_Coeff[0] * Lattice.dR(i) * X[i];

		for (int j = 1; j < Adams_N; j++)
		{
			f[i] += incr * Adams_Coeff[j] * df_dR[i - incr*j] * Lattice.dR(i - incr*j);
		}

		f[i] /= Det;

		df_dR[i] = A[i] * f[i] + X[i];
	}
}

Interpolation::Interpolation(int Order)
{
	order = Order;//fixed order, not too high, to avoid oscillations
}

vector<double> Interpolation::get_value(const vector<double> &f, const vector<double> &x_ini, double X)
{
	int close_left = 0;
	vector<double> Result(2, 0);
	vector<double> P(order + 1, 0);
	vector<double> x(order + 1, 0);
	vector<double> tmp(order, 0);
	vector<double> d_tmp(order, 0);
	vector<double> dP(order + 1, 0);

	if (!x_ini.empty() && !f.empty() && f.size() == x_ini.size())
	{
		if (X >= x_ini.front() || X < x_ini.back())
		{
			while (X > x_ini[close_left + 1]) { close_left++; }

			if (close_left < x_ini.size() - order - 1)//Interpolate forwards.
			{
				for (int i = 0; i <= order; i++)
				{
					P[i] = f[close_left + i];
					x[i] = x_ini[close_left + i];
				}
			}
			else//X is too close to the end point. Interpolate backwards.
			{
				for (int i = 0; i <= order; i++)
				{
					P[i] = f[close_left + 1 - i];
					x[i] = x_ini[close_left + 1 - i];
				}
			}

			for (int i = 1; i <= order; i++)
			{
				for (int j = 0; j <= order - i; j++)
				{
					tmp[j] = (P[j] * (x[j + i] - X) + P[j + 1] * (X - x[j])) / (x[j + i] - x[j]);
					d_tmp[j] = (dP[j] * (x[j + i] - X) - P[j] + dP[j + 1] * (X - x[j]) + P[j + 1]) / (x[j + i] - x[j]);
				}
				P = tmp;
				dP = d_tmp;
			}

			Result[0] = P[0];
			Result[1] = dP[0];

		}
	}

	return Result;
}

vector<double> Interpolation::get_value(PairFunction &S_old, Grid &Lattice_old, double X)
{
	int close_left = 0;
	vector<double> Result(2, 0);
	vector<double> P(order + 1, 0);
	vector<double> x(order + 1, 0);
	vector<double> tmp(order, 0);
	vector<double> d_tmp(order, 0);
	vector<double> dP(order + 1, 0);

	if (Lattice_old.size() > 0 && Lattice_old.size() == S_old.size())
	{
		if (X >= Lattice_old.R(0) || X < Lattice_old.R(Lattice_old.size() - 1))
		{
			while (X > Lattice_old.R(close_left + 1)) { close_left++; }

			if (close_left < Lattice_old.size() - order - 1)//Interpolate forwards.
			{
				for (int i = 0; i <= order; i++)
				{
					P[i] = S_old.F[close_left + i];
					dP[i] = 0;
					x[i] = Lattice_old.R(close_left + i);
				}
			}
			else//X is too close to the end point. Interpolate backwards.
			{
				for (int i = 0; i <= order; i++)
				{
					P[i] = S_old.F[close_left + 1 - i];
					dP[i] = 0;
					x[i] = Lattice_old.R(close_left + 1 - i);
				}
			}

			for (int i = 1; i <= order; i++)
			{
				for (int j = 0; j <= order - i; j++)
				{
					tmp[j] = (P[j] * (x[j + i] - X) + P[j + 1] * (X - x[j])) / (x[j + i] - x[j]);
					d_tmp[j] = (dP[j] * (x[j + i] - X) - P[j] + dP[j + 1] * (X - x[j]) + P[j + 1]) / (x[j + i] - x[j]);
				}
				P = tmp;
				dP = d_tmp;
			}

			Result[0] = P[0];
			Result[1] = dP[0];

		}
	}

	return Result;
}

int Interpolation::RecalcWF(RadialWF &S_old, Grid &Lattice_old, RadialWF &S_new, Grid &Lattice_new)
{
	if (S_old.size() == Lattice_old.size())
	{
		int j = 0;
		double Infinity = Lattice_old.R(S_old.pract_infinity());
		vector<double> Tmp(2, 0);
		S_new.clear();
		S_new.resize(Lattice_new.size());
		while (Lattice_new.R(j) < Infinity && j < Lattice_new.size() - 1)
		{
			Tmp = get_value(S_old, Lattice_old, Lattice_new.R(j));
			S_new.F[j] = Tmp[0];
			S_new.G[j] = Tmp[1];
			j++;
		}
		S_new.set_infinity(j);

		return 0;
	}
	else return 1;
}



void Interpolation::gaussian_sum(vector<double> & Vals, Grid & Lattice, int max_iter, double conv_toll)
{
  // Function 'Vals' (radial atomic density), defined on a mesh 'R' (with intervals 'dR) is interpolated with a sum of 'Order'
  // Gaussians as 
     
  // Vals(R) ~ S(R)
  // S(R) = sum_{n=1}^{Order} R^2 * b_i*exp(-a_i*R^2)
  // a_i >= 0,  b_i > 0 
  // is enforsed at every iteration.

  // The loss
  // Loss({a_i, b_i}) = int_{R(start)}^{R(end)} (Vals(R) - S(R))^2*dR 
  // is minimized iteratively.
  
  int start_pt = 0;
  int end_pt = Lattice.size()-1;
  while (Vals[end_pt] < 0.000001) end_pt--;

  vector<double> a(order, 0);
  vector<double> b(order, 0);
  // Gradients. First 'order' values are a[i], last are b[i].
  vector<double> grad_a(order, 0); 
  vector<double> grad_b(order, 0); 
  // Adam paramteres.
  vector<double> M(2*order, 0);
  vector<double> V(2*order, 0);
  double Beta1 = 0.9, Beta2 = 0.999, eta=1;
  double Beta1_pow_m = 1, Beta2_pow_m = 1, M_hat, V_hat, e = 1;
  
  double Loss = 0;
  Adams I(Lattice, 5);
  vector<double> aux_func(end_pt-start_pt + 1, 0);
  vector<double> Vals_pred(end_pt-start_pt + 1, 0);
  // Initialize parameters:
  double vals_norm = I.Integrate(&Vals, start_pt, end_pt);
  double tmp, exp_ar2, pred_norm, Lambda = 1, lr = 1e-3;

  // TODO: adjust Beta in geometric set from 4 to something 'order' dependent.
  for (int i = 0; i < order; i++) {
    if (i == 0) a[0] = 0.1;
    else a[i] = 3*a[i-1];

    for (int j = start_pt; j < end_pt; j++) {
      tmp = Lattice.R(j)*Lattice.R(j);
      aux_func[j] = tmp*exp(-a[i]*tmp);
    }
    tmp = I.Integrate(&aux_func, start_pt, end_pt);
    b[i] = vals_norm/tmp/order; // To normalize 'S(R)'.
  }
  
  double epsilon = 1; // relative error.
  
  for (int m = 0; m < max_iter; m++) {
    // Calculate current loss.
    Loss = 0;
    for (int j = start_pt; j < end_pt; ++j) {
      Vals_pred[j] = 0;
      tmp = Lattice.R(j);
      tmp *= tmp;
      for (int i = 0; i < order; i++) Vals_pred[j] += b[i]*tmp*exp(-a[i]*tmp);
      tmp = Vals_pred[j] - Vals[j];
      aux_func[j] = tmp*tmp;
    }
    Loss = I.Integrate(&aux_func, start_pt, end_pt);
    //Loss = I.Integrate(&aux_func, start_pt, end_pt);
    pred_norm = I.Integrate(&Vals_pred, start_pt, end_pt);

    //Loss += Lambda*fabs(pred_norm - vals_norm);

    // Calculate gradients.
    for (int i = 0; i < order; i++) {
      // Grad b.
      grad_b[i] = 0;

      for (int j = start_pt; j < end_pt; j++) {
        tmp = Lattice.R(j);
        tmp *= tmp;
        exp_ar2 = tmp*exp(-a[i]*tmp);
        aux_func[j] = (Vals_pred[j] - Vals[j])*exp_ar2;
      }
      if (pred_norm > vals_norm) tmp = 1;
      else tmp = -1;

      grad_b[i] = 2*I.Integrate(&aux_func, start_pt, end_pt) + tmp*Lambda*0.25*sqrt(M_PI/a[i])/a[i];

      // Grad_a.
      grad_a[i] = 0;
      for (int j = start_pt; j < end_pt; j++) {
        tmp = Lattice.R(j);
        tmp *= tmp;
        exp_ar2 = tmp*tmp*exp(-a[i]*tmp);
        aux_func[j] = (Vals_pred[j] - Vals[j])*exp_ar2;
      }

      grad_a[i] = b[i]*(-2*I.Integrate(&aux_func, start_pt, end_pt) - tmp*Lambda*0.375*sqrt(M_PI/a[i])/(a[i]*a[i]));
    }

    printf("Loss = %3.6f | Norm(S) = %3.6f | Norm(V) = %3.6f\n", Loss, pred_norm, vals_norm);
    
    for (int i = 0; i < order; i++) {
      tmp = grad_a[i]*lr;
      if (fabs(tmp) > a[i]*0.5) tmp = 0.5*tmp/fabs(tmp)*a[i];
      a[i] -= tmp;

      tmp = grad_b[i]*lr;
      if (fabs(tmp) > b[i]*0.5) tmp = 0.5*tmp/fabs(tmp)*b[i];
      b[i] -= tmp;
    }

    lr *= 0.98;
    /*
    // Adam algorithm for paramter optimization.
    Beta1_pow_m *= Beta1;
    Beta2_pow_m *= Beta2;

    for (int i = 0; i < order; i++) {
      M[i] = Beta1*M[i] + (1-Beta1)*dLoss_dw[i];
      V[i] = Beta2*V[i] + (1-Beta2)*dLoss_dw[i]*dLoss_dw[i];

      M_hat = M[i]/(1 - Beta1_pow_m);
      V_hat = V[i]/(1 - Beta2_pow_m);

      if (i < order) {
        tmp = a[i] - eta*M_hat/(sqrt(V_hat) + e);
        if (tmp < 0.01) a[i] /= 3;
        else a[i] = tmp;
      } else {
        tmp = b[i] - eta*M_hat/(sqrt(V_hat) + e);
        if (tmp < 0) b[i] /= 3;
        else b[i] = tmp;
      }
    }
    */
    // TODO: Assemble new utility function, print out it's value, test and denug the code.
    if (epsilon < conv_toll) break;
  }

}

GaussQuad::GaussQuad(int Order)
{
	GaussX.clear();
	GaussW.clear();
	order = Order;
	switch(order) {
	case 2:
		GaussW.assign(gaussW_2, gaussW_2 + 2);
		GaussX.assign(gaussX_2, gaussX_2 + 2);
		break;
	case 4:
		GaussW.assign(gaussW_4, gaussW_4 + 4);
		GaussX.assign(gaussX_4, gaussX_4 + 4);
		break;
	case 5:
		GaussW.assign(gaussW_5, gaussW_5 + 5);
		GaussX.assign(gaussX_5, gaussX_5 + 5);
		break;
	default:
		GaussW.assign(gaussW_10, gaussW_10 + 10);
		GaussX.assign(gaussX_10, gaussX_10 + 10);
		break;
	}
}

void GaussQuad::set_order(int Order)
{
	GaussX.clear();
	GaussW.clear();
	order = Order;
	switch(order) {
	case 2:
		GaussW.assign(gaussW_2, gaussW_2 + 2);
		GaussX.assign(gaussX_2, gaussX_2 + 2);
		break;
	case 4:
		GaussW.assign(gaussW_4, gaussW_4 + 4);
		GaussX.assign(gaussX_4, gaussX_4 + 4);
		break;
	case 5:
		GaussW.assign(gaussW_5, gaussW_5 + 5);
		GaussX.assign(gaussX_5, gaussX_5 + 5);
		break;
  case 10:
  	GaussW.assign(gaussW_10, gaussW_10 + 10);
		GaussX.assign(gaussX_10, gaussX_10 + 10);
		break;
	default:
		GaussW.assign(gaussW_13, gaussW_13 + 13);
		GaussX.assign(gaussX_13, gaussX_13 + 13);
		break;
	}
}

vector<double> GaussQuad::get_Gauss_X(double a, double b)
{
	// Linear mapping: X[i] = 0.5(b-a)*GaussX[i] + 0.5*(a+b).
	double l = 0.5*(a + b);
	vector<double> Result(GaussX.size(), l);
	for (int i = 0; i < GaussX.size(); i++) {
		Result[i] += 0.5*(b - a)*GaussX[i];
	}

	return Result;
}

double GaussQuad::Integrate(vector<double> &F, vector<double> &x, double a, double b)
{
	// Gauss quadrature integration formula.
	double Result = 0;
	if (GaussX.size() != x.size()) return 0;
	for (int i = 0; i < x.size(); i++) {
		Result += F[i]*GaussW[i];
	}
	Result *= 0.5*(b - a);

	return Result;
}
