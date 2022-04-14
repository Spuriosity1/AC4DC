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
#include <vector>>
#include "Constant.h"
// #include "IntegrateRateEquation.h"
#include "Input.h"
#include "MolInp.h"
#include "HartreeFock.h"
#include "DecayRates.h"
#include "Numerics.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <sstream>
#include <string>
#include <omp.h>
#include <algorithm>
#include "EigenSolver.h"
// #include "Plasma.h"
#include <utility>




// 'Mothership' class for calculating rate coefficients
// Does photoionisation, fluorescence, Auger decay and others
//


class ComputeRateParam
{
public:
	//Orbitals are HF wavefunctions. This configuration is an initial state.
	//Assuming there are no unoccupied states in initial configuration!!!
	ComputeRateParam(Grid &Lattice, std::vector<RadialWF> &Orbitals, Potential &U, Input & Inp, bool recalc=true) :
	 	lattice(Lattice), input(Inp), orbitals(Orbitals), u(U), recalculate(recalc) {
		};
	~ComputeRateParam();

	// Halfwidth = 5/Constant::Time -> 5 fs half width.
	int SolveFrozen(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & log);
	RateData::Atom SolvePlasmaBEB(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & log);
	// // Atomic.
	// int SetupAndSolve(ofstream & log);
	// // Molecular.
	// int SetupAndSolve(MolInp & Input, std::ofstream & log);

	//string CompareRates(std::string RateFile1, std::string RateFile2, std::ofstream & log);// Find the difference in rate equation using two different rates.


	int Symbolic(const std::string & input, const std::string & output);//convertes configuration indexes in human readable format
	int Charge(int Iconf);
	std::vector<double> PerturbMe(std::vector<RadialWF> & Virtual, double Dist, double Einit);
	std::vector<double> Secular(std::vector<RadialWF> & Virtual, double Dist, double Einit);

	int NumPath() { return dimension; }
	std::vector<double> generate_G();
	std::vector<double> Times() { return T; }
	std::vector<double> dTimes() { return dT; }
	std::vector<double> Probs(int i) { return P[i]; }
	std::vector<std::vector<double>> AllProbs() {return P;}


	bool SetupIndex(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & log);
	std::vector<std::vector<int>> Get_Indexes() { return Index; }

  // Atomic data containers.
	std::vector<std::vector<double>> density = std::vector<std::vector<double>>(0);

  	Grid & Atom_Mesh() { return lattice; }

protected:
	Grid & lattice;
	Input & input;
	std::vector<RadialWF> & orbitals;
	Potential& u;
	bool recalculate; // Flag to determine whether or not to force-recompute everything

	std::vector<CustomDataType::polarize> MixMe;
	int dimension;//number of configurations
	std::vector<std::vector<double>> charge;
	std::vector<double> T;// Time grid points.
	std::vector<double> dT;// Accurate differentials.
	std::vector<std::vector<double>> P;// P[i][m] is the probabilities of having configurations "i" at time T[m].
	std::vector<std::vector<int> > Index;
	int mapOccInd(std::vector<RadialWF> & Orbitals);// Inverse of what Index returns.

	// Returns LaTeX formatted electron config referred to by index i
	std::string InterpretIndex(int i);

	RateData::Atom Store;

	std::vector<CustomDataType::ffactor> FF;
	std::vector<int> hole_posit;

	int extend_I(std::vector<double>& Intensity, double new_max_T, double step_T);
    std::vector<double> generate_I(std::vector<double>& T, double I_max, double HalfWidth);
	std::vector<double> generate_T(std::vector<double>& dT);
	std::vector<double> generate_dT(int num_elem);
    double T_avg_RMS(std::vector<std::pair<double, int>> conf_RMS);
	double T_avg_Charge();

	static bool sortEIIbyInd(RateData::EIIdata A, RateData::EIIdata B) { return (A.init < B.init); }
	static bool sortRatesFrom(RateData::Rate A, RateData::Rate B) { return (A.from < B.from); }
	static bool sortRatesTo(RateData::Rate A, RateData::Rate B) { return (A.to < B.to); }
	// Keys allow to quickly find the required element. See the GenerateFromKeys().
	std::vector<int> RatesFromKeys;
	void GenerateRateKeys(std::vector<RateData::Rate> & ToSort);
};
