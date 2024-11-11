/**
 * @file ComputeRateParam.h
 * @brief 
 * @details 
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
#include "RadialWF.h"
#include "Grid.h"
#include "Potential.h"
#include <vector>
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

using namespace std;



// 'Mothership' class for calculating rate coefficients
// Does photoionisation, fluorescence, Auger decay and others
//


class ComputeRateParam
{
public:
	//Orbitals are HF wavefunctions. This configuration is an initial state.
	//Assuming there are no unoccupied states in initial configuration!!!
	ComputeRateParam(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Inp, bool recalc=true, bool secondary_ion = true) :
	 	lattice(Lattice), input(Inp), orbitals(Orbitals), u(U), recalculate(recalc){
		};
	~ComputeRateParam();

	// Halfwidth = 5/Constant::Time -> 5 fs half width.
	RateData::Atom SolveAtomicRatesAndPlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, vector<bool> shell_check, bool calculate_secondary_ionisation, ofstream & log);
	// // Atomic.
	// int SetupAndSolve(ofstream & log);
	// // Molecular.
	// int SetupAndSolve(MolInp & Input, ofstream & log);

	//string CompareRates(string RateFile1, string RateFile2, ofstream & log);// Find the difference in rate equation using two different rates.


	int Symbolic(const string & input, const string & output);//convertes configuration indexes in human readable format
	int Charge(int Iconf);
	vector<double> PerturbMe(vector<RadialWF> & Virtual, double Dist, double Einit);
	vector<double> Secular(vector<RadialWF> & Virtual, double Dist, double Einit);

	int NumPath() { return dimension; }
	vector<double> generate_G();
	vector<double> Times() { return T; }
	vector<double> dTimes() { return dT; }
	vector<double> Probs(int i) { return P[i]; }
	vector<vector<double>> AllProbs() {return P;}

	bool SetupIndex(vector<int> Max_occ, vector<int> Final_occ, ofstream & log);
	vector<vector<int>> Get_Indexes() { return Index; }

  // Atomic data containers.
	vector<vector<double>> density = vector<vector<double>>(0);

  	Grid & Atom_Mesh() { return lattice; }

protected:
	Grid & lattice;
	Input & input;
	vector<RadialWF> & orbitals;
	Potential& u;
	bool recalculate; // Flag to determine whether or not to force-recompute everything
	vector<CustomDataType::polarize> MixMe;
	int dimension;//number of configurations
	vector<vector<double>> charge;
	vector<double> T;// Time grid points.
	vector<double> dT;// Accurate differentials.
	vector<vector<double>> P;// P[i][m] is the probabilities of having configurations "i" at time T[m].
	vector<vector<int> > Index;
	int mapOccInd(vector<RadialWF> & Orbitals);// Inverse of what Index returns.

	// Returns LaTeX formatted electron config referred to by index i
	string InterpretIndex(int i);

	RateData::Atom Store;

	vector<CustomDataType::ffactor> FF;
	vector<int> hole_posit;

	int extend_I(vector<double>& Intensity, double new_max_T, double step_T);
    vector<double> generate_I(vector<double>& T, double I_max, double HalfWidth);
	vector<double> generate_T(vector<double>& dT);
	vector<double> generate_dT(int num_elem);
    double T_avg_RMS(vector<pair<double, int>> conf_RMS);
	double T_avg_Charge();

	static bool sortEIIbyInd(RateData::EIIdata A, RateData::EIIdata B) { return (A.init < B.init); }
	static bool sortRatesFrom(RateData::Rate A, RateData::Rate B) { return (A.from < B.from); }
	static bool sortRatesTo(RateData::Rate A, RateData::Rate B) { return (A.to < B.to); }
	// Keys allow to quickly find the required element. See the GenerateFromKeys().
	vector<int> RatesFromKeys;
	void GenerateRateKeys(vector<RateData::Rate> & ToSort);
};
