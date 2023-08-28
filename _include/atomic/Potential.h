/**
 * @file Potential.h
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
#pragma once

#include "Grid.h"
#include "RadialWF.h"
#include <string>
#include <vector>
#include <memory>

using namespace std;

class Potential
{
public:
	Potential(Grid * Lattice, int Z, string mod = "coulomb", double Rad_well = 0);
	Potential(int i = 0) {}
	~Potential(void) {};

	// Direct and Current orbital exchange in HF approximation.
	int HF_upd_dir(RadialWF * Current, vector<RadialWF> & Orbitals);
	// HFS local potential.
	int LDA_upd_dir(vector<RadialWF> & Orbitals);
	// Exchange except self-interaction (included in Direct) of Current.
	int HF_upd_exc(RadialWF * Current, vector<RadialWF> & Orbitals);
	// Updates potential for virtual states in V_(N-1) approximation.
	int HF_V_N1(RadialWF * Current, vector<RadialWF> & Orbitals, int c, bool UpdDir, bool UpdExc);
	void Reset(void);
	void ScaleNucl(double Scl_dir);

	void NewLattice(Grid * Lattice);
	void GenerateTrial(vector<RadialWF> & Orbitals);

	//Radial Coulomb integral dr1 dr2 P_a(1) P_b(2) (r<)^k/(r>)^{k+1} P_c(1) P_d(2)
	double R_k(int k, RadialWF & A, RadialWF & B, RadialWF & C, RadialWF & D);
	std::vector<double> Y_k(int k, std::vector<double> density, int infinity, int L);
	double Overlap(std::vector<double> density, int infinity);
	vector<double> make_density(vector<RadialWF> & Orbitals);

	string Type();
	int NuclCharge() { return n_charge; }
	double R_well() { return r_well; }
	vector<double> V;//Nuclear + Direct
	vector<double> Exchange;
	vector<double> Trial;
	vector<double> LocExc;// Local density exchange approximation with tail correction.
	vector<double> Asympt;// Tail correction for LDA.
	vector<double> v0_N1;
	// Replace V[j] + LocExch[j] by Asympt[j] when first exceeds second manually.

	// Auxillary functions.
	vector<float> Get_Kinetic(vector<RadialWF> & Orbitals, int start_with = 0);

	const Potential& operator = (const Potential& Other) {
		V = Other.V;
		Exchange = Other.Exchange;
		Trial = Other.Trial;
		LocExc = Other.LocExc;
		Asympt = Other.Asympt;
		v0_N1 = Other.v0_N1;
		model = Other.model;
		n_charge = Other.n_charge;
		nuclear = Other.nuclear;
		r_well = Other.r_well;
		delete lattice;
		lattice = Other.lattice;
		return *this;
	}

protected:
	void GenerateNuclear(void);
	string model = "coulomb";
	int n_charge = 1;
	double r_well = 0.0000001;
	std::vector<double> nuclear;
	Grid * lattice = nullptr;
};


class MatrixElems
{
public:
	// Class for evaluation of various matrix elements.
	MatrixElems(Grid * Lattice);

	// Radial Coulomb integral dr1 dr2 P_a(1) P_b(2) (r<)^k/(r>)^{k+1} P_c(1) P_d(2)
//	double R_k(int k, RadialWF & A, RadialWF & B, RadialWF & C, RadialWF & D);
	// Reduced dipole dr P_a(r)P_b(r) r . Gauge can be either "length" or "velocity".
	double Dipole(RadialWF & A, RadialWF & B, string gauge);
	// Average over configuration dipole matrix element.
	double DipoleAvg(RadialWF & A, RadialWF & B, string gauge);
  // RMS radius of a slater determinant.
  double R_pow_k(vector<RadialWF> & Orbitals, int k);

	~MatrixElems() {};
private:
	double Msum(int La, int Lb, int k);
	Grid * lattice;
};
