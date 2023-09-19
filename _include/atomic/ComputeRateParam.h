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
#include "HFInputParam.h"
#include "RateData.h"
#include "HartreeFock.h"
#include "DecayRates.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <sstream>
#include <string>
#include <omp.h>
#include <algorithm>
#include "numerics/EigenSolver.h"
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
	ComputeRateParam(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, InputData::HFInputParam & Inp, ofstream& _log) :
	 	lattice(Lattice), input(Inp), orbitals(Orbitals), u(U), log(_log) {
		};
	~ComputeRateParam();

	void configure_calc(bool _calc_auger, bool _calc_fluor, bool _calc_photo, bool _calc_eii, bool _calc_FT, bool _calc_bound){
		this->calc_auger = _calc_auger;
		this->calc_fluor = _calc_fluor;
		this->calc_photo = _calc_photo;
		this->calc_eii = _calc_eii;
		this->calc_FT = _calc_FT;
		this->calc_bound_transport = _calc_bound;
	}

	// Halfwidth = 5/Constant::Time -> 5 fs half width.
	// int SolveFrozen(vector<int> Max_occ, vector<int> Final_occ, ofstream & log);
	RateData::Atom SolvePlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, vector<bool> shell_check);

	// std::vector<RateData::photo> calc_photo(const std::vector<int>& Max_occ, const std::vector<int> Final_occ, vector<bool> shell_check);
	// std::vector<RateData::fluor> calc_fluor(const std::vector<int>& Max_occ, const std::vector<int> Final_occ, vector<bool> shell_check);
	// std::vector<RateData::auger> calc_auger(const std::vector<int>& Max_occ, const std::vector<int> Final_occ, vector<bool> shell_check);
	// std::vector<RateData::EIIdata> calc_eii(const std::vector<int>& Max_occ, const std::vector<int> Final_occ, vector<bool> shell_check);
	


	// // Atomic.
	// int SetupAndSolve(ofstream & log);
	// // Molecular.
	// int SetupAndSolve(MolInp & Input, ofstream & log);

	//string CompareRates(string RateFile1, string RateFile2, ofstream & log);// Find the difference in rate equation using two different rates.


	// int Symbolic(const string & input, const string & output);//convertes configuration indexes in human readable format
	int Charge(int Iconf);
	// vector<double> PerturbMe(vector<RadialWF> & Virtual, double Dist, double Einit);
	// vector<double> Secular(vector<RadialWF> & Virtual, double Dist, double Einit);

	int NumPath() { return dimension; }
	// vector<double> generate_G();
	vector<double> Times() { return T; }
	vector<double> dTimes() { return dT; }
	vector<double> Probs(int i) { return P[i]; }
	vector<vector<double>> AllProbs() {return P;}

	bool SetupIndex(vector<int> Max_occ, vector<int> Final_occ);
	vector<vector<unsigned>> Get_Indexes() { return Index; }

  // Atomic data containers.
	vector<vector<double>> density = vector<vector<double>>(0);

  	Grid & Atom_Mesh() { return lattice; }

protected:

	bool calc_auger, calc_fluor, calc_photo, calc_eii, calc_FT, calc_bound_transport;

	Grid & lattice;
	InputData::HFInputParam & input;
	vector<RadialWF> & orbitals;
	Potential& u;
	ofstream& log;


	vector<RateData::polarize> MixMe;
	int dimension;//number of configurations
	vector<vector<double>> charge;
	vector<double> T;// Time grid points.
	vector<double> dT;// Accurate differentials.
	vector<vector<double>> P;// P[i][m] is the probabilities of having configurations "i" at time T[m].
	vector<vector<unsigned> > Index;
	int mapOccInd(vector<RadialWF> & Orbitals);// Inverse of what Index returns.

	// Returns LaTeX formatted electron config referred to by index i
	string InterpretIndex(int i);

	RateData::Atom Store;

	vector<RateData::ffactor> FF;
	vector<int> hole_posit;

	// int extend_I(vector<double>& Intensity, double new_max_T, double step_T);
    // vector<double> generate_I(vector<double>& T, double I_max, double HalfWidth);
	// vector<double> generate_T(vector<double>& dT);
	// vector<double> generate_dT(int num_elem);
    // double T_avg_RMS(vector<pair<double, int>> conf_RMS);
	// double T_avg_Charge();

	static bool sortEIIbyInd(RateData::EIIdata A, RateData::EIIdata B) { return (A.init < B.init); }
	static bool sortRatesFrom(RateData::Rate A, RateData::Rate B) { return (A.from < B.from); }
	static bool sortRatesTo(RateData::Rate A, RateData::Rate B) { return (A.to < B.to); }
	// Keys allow to quickly find the required element. See the GenerateFromKeys().
	vector<int> RatesFromKeys;
	void GenerateRateKeys(vector<RateData::Rate> & ToSort);
};
