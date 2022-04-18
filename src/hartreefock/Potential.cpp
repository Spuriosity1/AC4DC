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
#include "Potential.hpp"
#include "stdafx.hpp"
#include "Constant.hpp"
#include "Numerics.hpp"
#include <algorithm>
#include "Wigner/wignerSymbols.hpp"

using namespace std;

Potential::Potential(Grid * Lattice, int Z, std::string mod, double Rad_well) : lattice(Lattice)//Nuclear
{
	V.clear();
	V.resize(Lattice->size());
	nuclear.clear();
	nuclear.resize(Lattice->size());
	Exchange.clear();
	Exchange.resize(Lattice->size());
	Trial.clear();
	Trial.resize(Lattice->size());
	LocExc.clear();
	LocExc.resize(Lattice->size());
	v0_N1.clear();
	v0_N1.resize(Lattice->size());
	model = mod;
	n_charge = Z;

	r_well = Rad_well;

	GenerateNuclear();
}

void Potential::GenerateNuclear(void)
{
	if (model == "coulomb") {
		for (int i = 0; i < lattice->size(); i++) {
			nuclear[i] = -n_charge / lattice->R(i);
		}
	} else {
		for (int i = 0; i < lattice->size(); i++) {
			if (lattice->R(i) <= r_well) {
				nuclear[i] = -n_charge*(3.0 - pow((lattice->R(i) / r_well), 2)) / 2 / r_well;
			} else {
				nuclear[i] = -n_charge / lattice->R(i);
			}
		}
	}
}

void Potential::GenerateTrial(std::vector<RadialWF> & Orbitals)
{
	Trial.clear();
	Trial.resize(lattice->size());
	Asympt.clear();
	Asympt.resize(lattice->size());

	int n_elec = 0;
	for (auto& v : Orbitals)
	{
		n_elec += v.occupancy();
	}

	if (n_elec > 1)
	{
		double KOP = n_charge - n_elec;
		double AM = 3.1, DD = 0.4, POT = 2.235;
		double HH = exp(-AM / DD);

		POT = POT * pow(n_charge / 55., 1. / 3.);
		if (KOP == 0.) { KOP = 1.; }

		for (int i = 0; i < lattice->size(); i++)
		{
			double UR = lattice->R(i) / DD;
			double ZR;

			if (UR <= 20.)
			{
				ZR = (n_charge - KOP) / ((1. + POT*(lattice->R(i))) * (1. + POT*(lattice->R(i))) * (HH + 1.));
				ZR = ZR / (HH * exp(UR) + 1.) + KOP;
			}
			else { ZR = KOP; }

			V[i] = -ZR / lattice->R(i);
			Trial[i] = V[i];
		}
	}
	else
	{
		for (int i = 0; i < lattice->size(); i++) V[i] = nuclear[i];
	}
}

void Potential::Reset()//clear the Direct potential V and sets and equates it to the Nuclear
{
	GenerateNuclear();
}

void Potential::ScaleNucl(double Scl_dir)
{
	double Scaling = Scl_dir / n_charge;
	for (int i = 0; i < lattice->size(); i++) {
		V[i] += Scaling * nuclear[i];
	}
}

std::string Potential::Type()
{
	return model;
}

int Potential::HF_upd_dir(RadialWF* Current, std::vector<RadialWF> &Orbitals)
{
	// Hartree-Fock Direct + Orbital self-interaction exchange. If Current is that orbital, than both
	// exchange and direct are included. If Current is some other orbital, only the direct part of
	// the potential is evaluated.
	std::vector<double> y_0(lattice->size(), 0.);
	std::vector<double> y(lattice->size(), 0.);
	std::vector<double> density(lattice->size(), 0.);
	std::vector<double> density_current(lattice->size(), 0.);
	int infinity = 0;
	double Q = 1; // electron-electron interaction weight.
	int L_min = Orbitals[0].L();
	double angular = 0;
	int N_elec = 0;

	for (auto& Orb: Orbitals) N_elec += Orb.occupancy();

	if (N_elec < 2) {
		V = nuclear;
		return 0;
	}

	for (int i = 0; i < Orbitals.size(); i++)
	{
		if (Orbitals[i].occupancy() == 0) continue;
		else Q = Orbitals[i].occupancy();

		// Line is giving 'uninitialised value' when run under valgrind when called by upd_HF_dir.
		if (Current->L() == Orbitals[i].L() && Current->N() == Orbitals[i].N())	{
			if (Current->occupancy() > 1) {
				Q = Orbitals[i].occupancy() - 1;
				for (int j = 0; j < Orbitals[i].pract_infinity(); j++) {
					density_current[j] = Q * Orbitals[i].F[j] * Orbitals[i].F[j];
				}

				y = Y_k(0, density_current, Orbitals[i].pract_infinity(), 2 * Current->L());

				//for (std::vector<double>::iterator Y = y.begin(); Y != y.end(); ++Y) *Y /= Q;

				if (Current->L() > 0) {
					for (int k = 2; k <= 2 * Orbitals[i].L(); k += 2) {
						angular = Constant::Wigner3j(Orbitals[i].L(), k, Orbitals[i].L(), 0, 0, 0);
						angular = 0.5*angular*angular*(4*Current->L() + 2)/(4*Current->L() + 1);

						y_0 = Y_k(k, density_current, Orbitals[i].pract_infinity(), 2 * Current->L());
					 	for (int j = 0; j < lattice->size(); j++)	{
							y[j] -= angular*y_0[j];
						}
					}
				}
				if (infinity < Orbitals[i].pract_infinity()) { infinity = Orbitals[i].pract_infinity(); }
			}
			continue;//orbital self inteaction is calculated exactly
		}

		if (Orbitals[i].L() < L_min) {
			L_min = Orbitals[i].L();
		}

		for (int j = 0; j < Orbitals[i].pract_infinity(); j++) {
			density[j] += Q * Orbitals[i].F[j] * Orbitals[i].F[j];
		}

		if (infinity < Orbitals[i].pract_infinity()) {
			infinity = Orbitals[i].pract_infinity();
		}
	}

	if (infinity != 0)	y_0 = Y_k(0, density, infinity, 2 * L_min);

	for (int i = 0; i < lattice->size(); i++) {
		V[i] = nuclear[i] + (y_0[i] + y[i]) / lattice->R(i);
		if (density[i] != 0) LocExc[i] = -0.635348143*pow((density[i] / lattice->R(i) / lattice->R(i)), 1. / 3.);
	}

	return 0;
}

int Potential::LDA_upd_dir(std::vector<RadialWF> &Orbitals)
{
	//updates Hartree-Fock self-consistent field at every iteration. Function recalculates direct and exchange potential.
	//Local exchange is also updated here.
	std::vector<double> y_0(lattice->size(), 0.);
	std::vector<double> density(lattice->size(), 0.);
	int infinity = 0;
	int L_min = Orbitals[0].L();
	int N_elec = 0;
	double Q = 1;

	LocExc.clear();
	LocExc.resize(lattice->size());
	Asympt.clear();
	Asympt.resize(lattice->size());

	for (int i = 0; i < Orbitals.size(); i++)
	{
		if (Orbitals[i].occupancy() == 0) continue;
		else Q = Orbitals[i].occupancy();
		if (Orbitals[i].L() < L_min) { L_min = Orbitals[i].L(); }

		for (int j = 0; j < Orbitals[i].pract_infinity(); j++) {
			density[j] += Q * Orbitals[i].F[j] * Orbitals[i].F[j];
		}

		N_elec += Orbitals[i].occupancy();
		if (infinity < Orbitals[i].pract_infinity()) infinity = Orbitals[i].pract_infinity();
	}

	if (N_elec == 1) {
		V = nuclear;
		for (int i = 0; i < lattice->size(); i++) {
			LocExc[i] = 0;
			Asympt[i] = 0;
		}
		return 0;
	}
	if (infinity != 0)	y_0 = Y_k(0, density, infinity, 2 * L_min);
	int Z_eff = N_elec - n_charge - 1;

	double V_tmp = 0;
	if (model == "coulomb") {
		for (int i = 0; i < lattice->size(); i++) {
			Asympt[i] = Z_eff / lattice->R(i);
			LocExc[i] = -0.635348143*pow((density[i] / lattice->R(i) / lattice->R(i)), 1. / 3.);
			V[i] = nuclear[i] + y_0[i] / lattice->R(i) + LocExc[i];// Nuclear + Direct potential.
			if (V[i] > Asympt[i]) V[i] = Asympt[i];// HFS/LDA tail correction.
		}
	}
	else {
		for (int i = 0; i < lattice->size(); i++) {
			if (lattice->R(i) <= r_well) {
				Asympt[i] = -Z_eff*(3.0 - pow((lattice->R(i) / r_well), 2)) / 2 / r_well;
			} else {
				Asympt[i] = -Z_eff / lattice->R(i);
			}
			LocExc[i] = -0.635348143*pow((density[i] / lattice->R(i) / lattice->R(i)), 1. / 3.);
			V[i] = nuclear[i] + y_0[i] / lattice->R(i) + LocExc[i];// Nuclear + Direct potential.
			if (V[i] > Asympt[i]) V[i] = Asympt[i];// HFS/LDA tail correction.
		}
	}

	return 0;
}

int Potential::HF_upd_exc(RadialWF * Current, std::vector<RadialWF> &Orbitals)
{
	// Function recalculates Exchange potential for Current wavefunction in the field of Orbitals.
	std::vector<double> y_k(lattice->size(), 0.);
	std::vector<double> density_exchange(lattice->size(), 0.);
	int infinity = lattice->size()-1;
	double angular, Q = 1;
	int N_elec = 0;

	for (int i = 0; i < lattice->size(); i++) {
		Exchange[i] = 0;
	}
	for (auto& Orb: Orbitals) N_elec += Orb.occupancy();

	if (N_elec < 2) return 0;

	for (int i = 0; i < Orbitals.size(); i++)
	{
		if (Orbitals[i].occupancy() == 0) { continue; }
		if (Current->L() == Orbitals[i].L() && Current->N() == Orbitals[i].N())	continue;
		Q = Orbitals[i].occupancy();

		if (Current->pract_infinity() < Orbitals[i].pract_infinity()) infinity = Current->pract_infinity();
		else infinity = Orbitals[i].pract_infinity();

		for (int j = 0; j <= infinity; j++)	density_exchange[j] = Q * Current->F[j] * Orbitals[i].F[j];

		for (int k = abs(Current->L() - Orbitals[i].L()); k <= (Current->L() + Orbitals[i].L()); k++)
		{
			if ((Current->L() + k + Orbitals[i].L()) % 2 != 0) continue;
			angular = Constant::Wigner3j(Current->L(), k, Orbitals[i].L(), 0, 0, 0);
			angular = 0.5*angular*angular;

			y_k = Y_k(k, density_exchange, infinity, Current->L() + Orbitals[i].L());

			for (int j = 0; j < lattice->size(); j++) Exchange[j] -= angular * y_k[j] * Orbitals[i].F[j] / lattice->R(j);
		}
	}
	return 0;
}

int Potential::HF_V_N1(RadialWF * Current, std::vector<RadialWF> & Orbitals, int c, bool UpdDir, bool UpdExc)
{
	// Check if both switches are false. If so no action is taken.
	// Required to insure no multiple subtraction from V & Exchange.
	if (!UpdDir && !UpdExc) return 0;

	int N_elec = 0;
	for (auto& Orb: Orbitals) N_elec += Orb.occupancy();

	std::vector<double> density(lattice->size(), 0);
	// V_(N-1) approximation for accurate virtual states. Has no effect on core and shouldn't be calculated for core orbitals.
	// Follows W. Johnson p. 127
	if (v0_N1[0] == 0 && N_elec > 1) {
		// First use. Define v0 as Y_0(c, c, r) / r, where "c" is a core state orbital, so 0 <= c <= Orbtials.size().

		for (int i = 0; i <= Orbitals[c].pract_infinity(); i++) {
			density[i] = Orbitals[c].F[i] * Orbitals[c].F[i];
		}
		v0_N1 = Y_k(0, density, Orbitals[c].pract_infinity(), 2*Orbitals[c].L());
		for (int i = 0; i < lattice->size(); i++) v0_N1[i] /= lattice->R(i);
	}

	if (UpdDir) {
		for (int i = 0; i < lattice->size(); i++) V[i] -= v0_N1[i];
	}

	if (UpdExc) {
		Adams I(*lattice, 10);
		double MatrElem = 0;
		int infty = Current->pract_infinity();
		for (auto& Orb: Orbitals) {
			if (Orb.L() != Current->L()) continue; // Sperically symmetric potential in MatrElem.

			infty = Current->pract_infinity();
			if (Orb.pract_infinity() < infty) infty = Orb.pract_infinity();

			for (int i = 0; i < lattice->size(); i++) {
				if (i <= infty) density[i] = Current->F[i] * v0_N1[i] * Orb.F[i];
				else density[i] = 0;
			}
			MatrElem = I.Integrate(&density, 0, infty);

			for (int i = 0; i <= infty; i++) Exchange[i] += MatrElem*Orb.F[i];
		}
	}

	return 0;
}

double Potential::R_k(int k, RadialWF &A, RadialWF &B, RadialWF &C, RadialWF &D)
{
	//Radial Coulomb integral dr1 dr2 P_a(1) P_b(2) (r<)^k/(r>)^{k+1} P_c(1) P_d(2)
	double Result = 0;
	int infinity = min(B.pract_infinity(), D.pract_infinity());
	std::vector<double> density(lattice->size());

	for (int i = 0; i <= infinity; i++) { density[i] = B.F[i] * D.F[i]; }

	density = Y_k(k, density, infinity, B.L() + D.L());
	infinity = min(A.pract_infinity(), C.pract_infinity());
	for (int i = 0; i <= infinity; i++) { density[i] *= A.F[i] * C.F[i] / lattice->R(i); }

	Adams I(*lattice, 10);
	Result = I.Integrate(&density, 0, infinity);

	return Result;
}

std::vector<double> Potential::Y_k(int k, std::vector<double> density, int infinity, int L)
{
	//density.size() = lattice->size(); start_pt=0, end_pt=lattice->size()
	std::vector<double> Result;
	std::vector<double> Y_less(infinity+1, 0.);
	std::vector<double> Y_gtr(infinity+1, 0.);
	int adams_order = 10;

	Y_less[0] = density[0] * lattice->R(0)/(L+3);
	Y_gtr[infinity] = density[infinity] * lattice->dR(infinity);

	Adams W(*lattice, adams_order);
	for (int i = 0; i < density.size(); i++) {
		W.A[i] = - k / lattice->R(i);
		W.X[i] = density[i];
	}
	// Rough integration.
	for (int i = 0; i < adams_order; i++) {
		Y_less[i + 1] = (W.A[i] * Y_less[i] + W.X[i])*lattice->dR(i) + Y_less[i];
		Y_less[i + 1] = 0.5*((W.A[i+1] * Y_less[i+1] + W.X[i+1])*lattice->dR(i+1) + Y_less[i+1]);
	}

	W.Integrate_ODE(Y_less, 0, infinity);

	for (int i = 0; i < lattice->size(); i++) {
		W.A[i] = (k+1) / lattice->R(i);
		W.X[i] = -density[i];
	}
	// Rough integration.
	for (int i = infinity; i > infinity - adams_order; i--) {
		Y_gtr[i - 1] = -(W.A[i] * Y_gtr[i] + W.X[i])*lattice->dR(i) + Y_gtr[i];
		Y_gtr[i - 1] = -0.5*((W.A[i-1] * Y_gtr[i-1] + W.X[i-1])*lattice->dR(i-1) - Y_gtr[i-1]);
	}

	W.Integrate_ODE(Y_gtr, infinity, 0);
	Result.resize(lattice->size());

	for (int i = 0; i < density.size(); i++)
	{
		if (i < infinity) { Result[i] = Y_less[i] + Y_gtr[i]; }
		else { Result[i] = Y_less[infinity]; }

	}

	return Result;
}

std::vector<double> Potential::make_density(std::vector<RadialWF> & Orbitals)
{
	std::vector<double> Result(lattice->size(), 0.);
	int infty = 0;
	double occ = 0;
	for (int i = 0 ; i < Orbitals.size(); i++) {
		infty = Orbitals[i].pract_infinity();
		occ = Orbitals[i].occupancy();
		double* f = Orbitals[i].F.data();
		for (int j = 0; j < Result.size(); j++) {
			if (j > infty) break;
			Result[j] += occ * *(f + j) * *(f + j);
		}
	}

	return Result;
}

double Potential::Overlap(std::vector<double> density, int infinity)
{
	double Result;
	Adams W(*lattice, 10);

	Result = W.Integrate(&density, 0, infinity);

	return Result;
}

std::vector<float> Potential::Get_Kinetic(std::vector<RadialWF> & Orbitals, int start_with)
{
	int size = 0;
	for (int i = start_with; i < Orbitals.size(); i++) if (Orbitals[i].occupancy() != 0) size++;

	if (size == 0) return std::vector<float>(0);

	std::vector<float> Result(size, 0);
	// Get Kinetic energies.
	std::vector<double> density(lattice->size(), 0.);
	int infinity = 0;
	Adams I(*lattice, 5);

	size = 0;
	for (int i = start_with; i < Orbitals.size(); i++) {
		if (Orbitals[i].occupancy() == 0) continue;
		infinity = Orbitals[i].pract_infinity();
		for (int j = 0; j <= infinity; j++) {
			density[j] = Orbitals[i].F[j] * Orbitals[i].F[j] * V[j];
		}
		Result[size] = float(Orbitals[i].Energy - I.Integrate(&density, 0, infinity));
		size++;
	}

	return Result;
}

MatrixElems::MatrixElems(Grid * Lattice) : lattice(Lattice)
{
}

double MatrixElems::ImpactCrossSection(RadialWF &A, RadialWF &B, std::string gauge)
{
	double Result = 0;

	int infty = A.pract_infinity();
	if (B.pract_infinity() > infty) infty = B.pract_infinity();

    std::vector<double> density(infty+1, 0);

    if (gauge == "length") {
		for (int i = 0; i < density.size(); i++) {
			density[i] = lattice->R(i)*A.F[i]*B.F[i];
		}
    }
    else {
		double ang_coeff =  0.5*(A.L() - B.L())*(A.L() + B.L() + 1);
			for (int i = 0; i < density.size(); i++)	{
			density[i] = B.F[i] *(A.G[i] + ang_coeff * A.F[i]/lattice->R(i)) ;
		}
    }

	Adams I(*lattice, 10);
	Result = I.Integrate(&density, 0, density.size()-1);

	return Result;
}

double MatrixElems::DipoleAvg(RadialWF & A, RadialWF & B, std::string gauge)
{
	double Result = 0;
	if (A.L() > B.L()) Result = sqrt((double)A.L());
	else Result = sqrt((double)B.L());

	if ( (A.L() + (A.L() + B.L() + 1)/2) % 2 == 1) Result = -1.*Result;
	Result *= Msum(A.L(), B.L(), 1)*Dipole(A, B, gauge);

	return Result;
}


double MatrixElems::Msum(int La, int Lb, int k)
{
	// Calculates Sum_m (-1)^m / La Lb k \
	//                         \ -m m 0 /
	// Used for average over configuration calculations of matrix elements g(abcd) - g(abdc).
    double Result = 0, tmp = 0;
    // Check selection rules.
    if ((La + Lb + k) % 2 != 0) return 0;
    if (La + Lb < k || La + k < Lb || Lb + k < La) return 0;
    if (abs(La - Lb) < k || abs(La - k) < Lb || abs(k - Lb) < La) return 0;
	// Summation.
	int maxM = La;
	if (La > Lb) maxM = Lb;
	for (int m = -maxM; m <= maxM; m++) {
		tmp = WignerSymbols::wigner3j(La, k, Lb, -m, 0, m);
		if (m%2 == 0) Result += tmp;
		else Result -= tmp;
	}

	return Result;
}


double MatrixElems::R_pow_k(std::vector<RadialWF> & Orbitals, int k)
{
  double Result = 0;

  int infty = 0;
  for (auto & orb : Orbitals) if (orb.pract_infinity() > infty) infty = orb.pract_infinity();

  std::vector<double> density(infty+1, 0);

  for (int i = 0; i < density.size(); i++) {
    for (auto & orb : Orbitals) density[i] += orb.occupancy()*orb.F[i]*orb.F[i];
    density[i] *= pow(lattice->R(i), k);
  }

	Adams I(*lattice, 10);
	Result = I.Integrate(&density, 0, density.size()-1);

  return Result;
}
