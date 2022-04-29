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
#ifndef COMPUTE_RATE_PARAM_H
#define COMPUTE_RATE_PARAM_H



#include "RadialWF.hpp"
#include "Grid.hpp"
#include "Potential.hpp"
#include <vector>
#include "Constant.hpp"
#include "RateData.hpp"
// #include "IntegrateRateEquation.hpp"
#include "HFInput.hpp"
// #include "MolInp.hpp"
#include "HartreeFock.hpp"
#include "DecayRates.hpp"
#include "Numerics.hpp"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <sstream>
#include <string>
#include <omp.h>
#include <algorithm>
#include "EigenSolver.hpp"
// #include "Plasma.hpp"
#include <utility>
#include <filesystem>




// 'Mothership' class for calculating rate coefficients
// Does photoionisation, fluorescence, Auger decay and others
//

class ComputeRateParam
{
public:
	//Orbitals are HF wavefunctions. This configuration is an initial state.
	//Assuming there are no unoccupied states in initial configuration!!!
	ComputeRateParam(Grid &Lattice, std::vector<RadialWF> &Orbitals, Potential &U, HFInput & Inp, unsigned num_threads = 4) :
	 	lattice(Lattice), input(Inp), orbitals(Orbitals), u(U), num_threads(num_threads) {
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
	// std::vector<double> generate_G();
	// std::vector<double> Times() { return T; }
	// std::vector<double> dTimes() { return dT; }
	// std::vector<double> Probs(int i) { return P[i]; }
	// std::vector<std::vector<double>> AllProbs() {return P;}


	bool SetupIndex(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & log);
	std::vector<std::vector<int>> Get_Indexes() { return Index; }

  // Atomic data containers.
	std::vector<std::vector<double>> density = std::vector<std::vector<double>>(0);

  	const Grid & Atom_Mesh() const  { return lattice; }

	void save_as(std::filesystem::path p){
		RateData::save_csv(p, Store);
	}

protected:
	const Grid & lattice;
	const HFInput & input;
	std::vector<RadialWF> & orbitals;
	const Potential& u;


	unsigned num_threads;

	std::vector<PhysicalRate::polarize> MixMe;
	int dimension;//number of configurations
	// std::vector<std::vector<double>> charge;

	
	std::vector<std::vector<int> > Index;
	int mapOccInd(std::vector<RadialWF> & Orbitals);// Inverse of what Index returns.

	// Returns LaTeX formatted electron config referred to by index i
	std::string InterpretIndex(int i);

	RateData::Atom Store;

	std::vector<PhysicalRate::ffactor> FF;
	std::vector<int> hole_posit;

	static bool sortEIIbyInd(RateData::EIIdata A, RateData::EIIdata B) { return (A.init < B.init); }
	static bool sortRatesFrom(RateData::Rate A, RateData::Rate B) { return (A.from < B.from); }
	static bool sortRatesTo(RateData::Rate A, RateData::Rate B) { return (A.to < B.to); }
	// Keys allow to quickly find the required element. See the GenerateFromKeys().
	std::vector<int> RatesFromKeys;
	void GenerateRateKeys(std::vector<RateData::Rate> & ToSort);
};

#endif