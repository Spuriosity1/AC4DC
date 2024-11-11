/**
 * @file DecayRates.cpp
 * 
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
#include "DecayRates.h"
#include "Numerics.h"
#include "Wigner/wignerSymbols.h"
#include <fstream>

using namespace Constant;

int IntegrateContinuum(Grid&, Potential&, vector<RadialWF>&, RadialWF*, bool);
double A_k(int k, int L, int l_h, int l_f, int l_e, int l_c);
bool Triad(int l_a, int l_b, int l_c);
//int IntegrateContinuumOnce(Grid&, Potential&, RadialWF*);

DecayRates::DecayRates(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Inp) : lattice(Lattice), orbitals(Orbitals),
u(U), input(Inp)
{
}

vector<photo> DecayRates::Photo_Ion(double omega, ofstream & log)
{
	vector<photo> Result(0);
	if (omega <= 0) return Result;

	photo PhotoTmp;
	int j = 0, N_elec = 0;
	int infinity, L;
	double Infinity = 0, k_min = sqrt(2*omega), k_max = 0, V_tmp = 0;

	vector<RadialWF> Orbitals(orbitals.size());

	for (int i = 0; i < orbitals.size(); i++)
	{
		Orbitals[i].set_N(orbitals[i].N());
		Orbitals[i].set_L(orbitals[i].L());        // Resets occupancy
		Orbitals[i].Energy = orbitals[i].Energy;
		if (orbitals[i].occupancy() != 0)
		{
			N_elec += orbitals[i].occupancy();
			if (j < orbitals[i].pract_infinity()) { j = orbitals[i].pract_infinity(); Infinity = lattice.R(j); }
			if (omega + orbitals[i].Energy >= 0.25)
			{
				if (0.5*k_min*k_min > omega + orbitals[i].Energy) k_min = sqrt(2 * (omega + orbitals[i].Energy));
				if (0.5*k_max*k_max < omega + orbitals[i].Energy) k_max = sqrt(2 * (omega + orbitals[i].Energy));
			}
		}
		Orbitals[i].set_occupancy(orbitals[i].occupancy());
	}
	if (k_max <= 0) return Result;
	if (Infinity < 5*6.28/k_min) Infinity = 5*6.28/k_min;
//	Usually a finer grid is required.
	Grid Lattice(lattice.R(0), Infinity, 0.03/k_max);
//	Grid Lattice(50000, 0.001, 50, "linear");

//	Interpolate orbitals on the new grid
	Interpolation W(6);

	for (int i = 0; i < Orbitals.size(); i++) {
		if (input.Exited_Pot_Model() != "V_N-1" && orbitals[i].occupancy() == 0) continue;
		if (orbitals[i].F[0] != 0) W.RecalcWF(orbitals[i], lattice, Orbitals[i], Lattice);
		Infinity = lattice.R(orbitals[i].pract_infinity());
		j = 0;
		while (Lattice.R(j) < Infinity && j < Lattice.size() - 1) j++;
		Orbitals[i].set_infinity(j);
	}

//	Potential on the new grid.
	Potential U(&Lattice, u.NuclCharge());
	U.GenerateTrial(Orbitals);

	RadialWF Continuum(Lattice.size());

	// These lines do not do anything, it just initialises bogus values of n and l
	// so that valgrind does not complain about the comparison in Potential.cpp.
	// n=-1 is chosen to distinguish it from bound states.
	Continuum.set_N(-1);
	Continuum.set_L(0); // bogus line to make valgrind happy
	Continuum.set_occupancy(1);
	Continuum.set_infinity(Lattice.size() - 1);

	Adams I(Lattice, 10);
	vector<double> density;
	double ME = 0;

// main loop over all the possible transitions
	for (int i = 0; i < Orbitals.size(); i++)
	{
		if (Orbitals[i].occupancy() > 0)
		{
			infinity = Orbitals[i].pract_infinity();

			density.clear();
			density.resize(infinity + 1);
			for (int s = 0; s < density.size(); s++) density[s] = pow(Orbitals[i].F[s], 2);
			ME = I.Integrate(&density, 0, infinity);

			j = Orbitals[i].occupancy();
			PhotoTmp.hole = i;
			PhotoTmp.val = 0;
			if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[i].set_occupancy(j - 1);//excite one electron...
			Continuum.Energy = omega + Orbitals[i].Energy;//... into continuum
			if (Continuum.Energy > 0)
			{
				// Update here if Input.Hamiltonian == "LDA".
				if (input.Hamiltonian() == 1) {
					U.LDA_upd_dir(Orbitals);
				} else {
					U.HF_upd_dir(&Continuum, Orbitals); // This line seems to give 'conditional jump or move depends on uninitialised values'
				}
				if (input.Exited_Pot_Model() == "V_N-1") U.HF_V_N1(&Continuum, Orbitals, i, true, false);

				for (int l = Orbitals[i].L() - 1; l <= Orbitals[i].L() + 1; l += 2)
				{
					if (l >= 0)
					{
						Continuum.set_L(l);   //TODO make sure this is fine - S.P.
						if (IntegrateContinuum(Lattice, U, Orbitals, &Continuum, i) < 0) {
						log << "====================================================================" << endl;
							log << "Continuum didn't converge: " << endl;
							for (int i = 0; i < orbitals.size(); i++)
							{
								log << i + 1 << ") n = " << orbitals[i].N() << " l = " << orbitals[i].L()
									<< " Energy = " << orbitals[i].Energy << " Occup = " << orbitals[i].occupancy() << endl;
							}
							log.flush();
						}

						density.clear();
						density.resize(infinity + 1);
						if (input.Gauge() == "length") {
								for (int s = 0; s < density.size(); s++)	{
								density[s] = Lattice.R(s)*Continuum.F[s] * Orbitals[i].F[s];
							}
							ME = I.Integrate(&density, 0, infinity);
						}
						else {
						// Velocity gauge.
							double ang_coeff =  0.5*(Orbitals[i].L() - Continuum.L())*(Orbitals[i].L() + Continuum.L() + 1);
								for (int s = 0; s < density.size(); s++)	{
								density[s] = Continuum.F[s] *(Orbitals[i].G[s] + ang_coeff * Orbitals[i].F[s]/Lattice.R(s)) ;
							}
							ME = I.Integrate(&density, 0, infinity)/omega;
						}

						if (Orbitals[i].L() > l) { L = Orbitals[i].L();	}
						else { L = l; }
						PhotoTmp.val += 4*Pi*Pi*Alpha*omega*j*ME*ME*L / ((2*Orbitals[i].L()+1) * 3);
					}
				}
			}
			if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[i].set_occupancy(j);
			Result.push_back(PhotoTmp);
		}
	}

	return Result;
}

// vectorized version...
vector<vector<photo>> DecayRates::Photo_Ion(std::vector<double> omega_array, ofstream & log)
{
	size_t dim = omega_array.size();
	vector<vector<photo>> Result_array(dim,vector<photo>(0));


	photo PhotoTmp;
	int j = 0, N_elec = 0;
	int infinity, L;
	//vector<double> Infinity(dim,0), V_tmp(dim,0);
	double Infinity = 0;
	vector<double> V_tmp(dim,0);
	//vector<bool> found_result(dim,0); // tracks when results have been found.
	vector<RadialWF> Orbitals(orbitals.size());

		std::vector<double> k_min(dim), k_max(dim,0);
		for (size_t m = 0; m < dim; m++){
			k_min[m] = sqrt(2*omega_array[m]);
		}
		for (int i = 0; i < orbitals.size(); i++)
		{
			Orbitals[i].set_N(orbitals[i].N());
			Orbitals[i].set_L(orbitals[i].L());        // Resets occupancy
			Orbitals[i].Energy = orbitals[i].Energy;
			
			for (size_t m = 0; m < dim; m++){
				if (orbitals[i].occupancy() != 0)
				{
					N_elec += orbitals[i].occupancy();
					if (j < orbitals[i].pract_infinity()) { j = orbitals[i].pract_infinity(); Infinity = lattice.R(j); }
					if (omega_array[m] + orbitals[i].Energy >= 0.25)
					{
						if (0.5*k_min[m]*k_min[m] > omega_array[m] + orbitals[i].Energy) k_min[m] = sqrt(2 * (omega_array[m] + orbitals[i].Energy));
						if (0.5*k_max[m]*k_max[m] < omega_array[m] + orbitals[i].Energy) k_max[m] = sqrt(2 * (omega_array[m] + orbitals[i].Energy));
					}
				}
			}

			Orbitals[i].set_occupancy(orbitals[i].occupancy());
		}
		double max_k_max = 0;
		for (size_t m = 0; m < dim; m++){
			//if (k_max[m] <= 0) found_result[m] = true;
			if (Infinity < 5*6.28/k_min[m]) Infinity = 5*6.28/k_min[m];
			if (max_k_max < k_max[m]) max_k_max = k_max[m];
		}
	//	Usually a finer grid is required.
	//	Grid Lattice(50000, 0.001, 50, "linear");
	Grid Lattice(lattice.R(0), Infinity, 0.03/max_k_max);

//	Interpolate orbitals on the new grid
	Interpolation W(6);

	for (int i = 0; i < Orbitals.size(); i++) {
		if (input.Exited_Pot_Model() != "V_N-1" && orbitals[i].occupancy() == 0) continue;
		if (orbitals[i].F[0] != 0) W.RecalcWF(orbitals[i], lattice, Orbitals[i], Lattice);
		Infinity = lattice.R(orbitals[i].pract_infinity());
		j = 0;
		while (Lattice.R(j) < Infinity && j < Lattice.size() - 1) j++;
		Orbitals[i].set_infinity(j);
	}

//	Potential on the new grid.
	Potential U(&Lattice, u.NuclCharge());
	U.GenerateTrial(Orbitals);

	RadialWF Continuum(Lattice.size());

	// These lines do not do anything, it just initialises bogus values of n and l
	// so that valgrind does not complain about the comparison in Potential.cpp.
	// n=-1 is chosen to distinguish it from bound states.
	Continuum.set_N(-1);
	Continuum.set_L(0); // bogus line to make valgrind happy
	Continuum.set_occupancy(1);
	Continuum.set_infinity(Lattice.size() - 1);

	Adams I(Lattice, 10);
	vector<double> density;
	double ME = 0;

// main loop over all the possible transitions, peformed largely independently of omega.
	std::vector<photo> omega_independent_Result(0);
	std::vector<std::vector<double>> continuum_energies(dim);
	for (int i = 0; i < Orbitals.size(); i++)
	{
		if (Orbitals[i].occupancy() > 0)
		{
			infinity = Orbitals[i].pract_infinity();

			density.clear();
			density.resize(infinity + 1);
			for (int s = 0; s < density.size(); s++) density[s] = pow(Orbitals[i].F[s], 2);
			ME = I.Integrate(&density, 0, infinity);

			j = Orbitals[i].occupancy();
			PhotoTmp.hole = i;
			PhotoTmp.val = 0;
			std::vector<photo> phototmps(dim,PhotoTmp);
			if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[i].set_occupancy(j - 1);//excite one electron...
			double max_continuum_energy = 0;
			for (size_t m = 0; m < dim; m++){
				double cont_energy = omega_array[m] + Orbitals[i].Energy;
				if (cont_energy > max_continuum_energy) max_continuum_energy = cont_energy;
			}
			Continuum.Energy = max_continuum_energy;
			if (Continuum.Energy > 0)
			{
				// Update here if Input.Hamiltonian == "LDA".
				if (input.Hamiltonian() == 1) {
					U.LDA_upd_dir(Orbitals);
				} else {
					U.HF_upd_dir(&Continuum, Orbitals); // This line seems to give 'conditional jump or move depends on uninitialised values'
				}
				if (input.Exited_Pot_Model() == "V_N-1") U.HF_V_N1(&Continuum, Orbitals, i, true, false);
				
				for (int l = Orbitals[i].L() - 1; l <= Orbitals[i].L() + 1; l += 2)
				{
					if (l >= 0)
					{
						for (size_t m = 0; m < dim; m++){
							Continuum.Energy = omega_array[m] + Orbitals[i].Energy;
							if (Continuum.Energy <= 0) continue;
							Continuum.set_L(l);   //TODO make sure this is fine - S.P.
							if (IntegrateContinuum(Lattice, U, Orbitals, &Continuum, i) < 0) {
							log << "====================================================================" << endl;
								log << "Continuum didn't converge: " << endl;
								for (int i = 0; i < orbitals.size(); i++)
								{
									log << i + 1 << ") n = " << orbitals[i].N() << " l = " << orbitals[i].L()
										<< " Energy = " << orbitals[i].Energy << " Occup = " << orbitals[i].occupancy() << endl;
								}
								log.flush();
							}

							density.clear();
							density.resize(infinity + 1);
							if (input.Gauge() == "length") {
									for (int s = 0; s < density.size(); s++)	{
									density[s] = Lattice.R(s)*Continuum.F[s] * Orbitals[i].F[s];
								}
								ME = I.Integrate(&density, 0, infinity);
							}
							else {
							// Velocity gauge.
								double ang_coeff =  0.5*(Orbitals[i].L() - Continuum.L())*(Orbitals[i].L() + Continuum.L() + 1);
									for (int s = 0; s < density.size(); s++)	{
									density[s] = Continuum.F[s] *(Orbitals[i].G[s] + ang_coeff * Orbitals[i].F[s]/Lattice.R(s)) ;
								}
								ME = I.Integrate(&density, 0, infinity)/omega_array[m]; //;
							}

							if (Orbitals[i].L() > l) { L = Orbitals[i].L();	}
							else { L = l; }
							phototmps[m].val += 4*Pi*Pi*Alpha*omega_array[m]*j*ME*ME*L / ((2*Orbitals[i].L()+1) * 3);; 
						}
					}
				}
			}
			if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[i].set_occupancy(j);
			for (size_t m = 0; m < dim; m++){
				Result_array[m].push_back(phototmps[m]);
			}
		}
	}
	// Fill in for all omegas TODO This logic is wrong at present.
	for (size_t m = 0; m < dim; m++){
		// Result_array[m] = omega_independent_Result;
		// for (int i = 0; i < omega_independent_Result.size(); i++)
		// {
		// 	if (continuum_energies[m][i] < 0){
		// 		Result_array[m][i].val = 0;
		// 	}
		// 	// We multiply by ME twice and multiply by omega.
		// 	if (input.Gauge() == "length")
		// 		Result_array[m][i].val *= omega_array[m];
		// 	else
		// 		Result_array[m][i].val /= omega_array[m];
				
		// }
		// Remove the zero elements
		Result_array[m].erase(
			std::remove_if(
				Result_array[m].begin(), 
				Result_array[m].end(),
				[](photo const & e) { return e.val == 0; }
			), 
			Result_array[m].end()
		); 
	}

	return Result_array;
}


vector<fluor> DecayRates::Fluor()
{
	vector<fluor> Result;
	vector<double> density;
	int N_h = 0;
	int N_j = 0;
	int L_max = 0;
	double ME = 0;
	fluor Tmp;

	Adams I(lattice, 10);

	for (int i = 0; i < orbitals.size(); i++)
	{
		if (orbitals[i].occupancy() < 4 * orbitals[i].L() + 2)//check if there is a hole
		{
			N_h = 4 * orbitals[i].L() + 2 - orbitals[i].occupancy();
			for (int j = i + 1; j < orbitals.size(); j++)
			{
				if (orbitals[j].L() == orbitals[i].L() + 1 || orbitals[j].L() == orbitals[i].L() - 1)//selection rules
				{
					N_j = orbitals[j].occupancy();
					if (N_j > 0 && orbitals[j].Energy - orbitals[i].Energy > 0)//check if there are any electrons
					{
						Tmp.hole = i;
						Tmp.fill = j;
						if (orbitals[i].pract_infinity() > orbitals[j].pract_infinity()) { L_max = orbitals[j].pract_infinity(); }
						else { L_max = orbitals[i].pract_infinity();; }
						density.clear();
						density.resize(L_max);

						if (input.Gauge() == "length") {
							for (int s = 0; s < density.size(); s++)	{
								density[s] = lattice.R(s) * orbitals[i].F[s] * orbitals[j].F[s];
							}
							ME = I.Integrate(&density, 0, density.size()-1);// Electric dipole matrix element
						}
						else {
						// Velocity gauge.
							double ang_coeff =  0.5*(orbitals[i].L() - orbitals[j].L())*(orbitals[i].L() + orbitals[j].L() + 1);
								for (int s = 0; s < density.size(); s++)	{
								density[s] = orbitals[j].F[s] *(orbitals[i].G[s] + ang_coeff * orbitals[i].F[s]/lattice.R(s)) ;
							}
							ME = I.Integrate(&density, 0, density.size()-1)/(orbitals[j].Energy - orbitals[i].Energy);
						}

						if (orbitals[i].L() > orbitals[j].L()) { L_max = orbitals[i].L(); }
						else { L_max = orbitals[j].L(); }
						Tmp.val = 2 * N_h * N_j * pow(Alpha*(orbitals[j].Energy - orbitals[i].Energy), 3)*L_max / (3 * (2 * orbitals[i].L() + 1)* (2 * orbitals[j].L() + 1))*ME*ME;
						Result.push_back(Tmp);
					}
				}
			}
		}
	}

	return Result;
}

void DecayRates::set_Max_occ(vector<int> & Max_occ, vector<RadialWF> orbitals){
	if (Max_occ.size() != orbitals.size()) {
		Max_occ.resize(orbitals.size());
		for (size_t i = 0; i < orbitals.size(); i++) {
			Max_occ[i] = max(4*orbitals[i].L() + 2,orbitals[i].occupancy()); // The second term in the max() function in case of shell approximation. We assume a user would never want to start a simulation with electrons missing from any shell other than the valence one, and this also assumes that electrons being excited into higher shells is insignificant.
		}
	}
}

vector<auger> DecayRates::Auger(vector<int> Max_occ, ofstream & log)
{
	vector<auger> Result;
	auger Tmp;
	int N_h;
	double N_fe, E_cont, V_tmp, Infinity;
	assert(Max_occ.size() == orbitals.size());
	//Create new lattice, core and potential suitable for continuum states
	double k = 0;
	for (int i = 0; i < orbitals.size(); i++)//find k that corresponds to shortest wavelength
	{
		if (orbitals[i].occupancy() < Max_occ[i]) {
			k = sqrt(-2 * orbitals[i].Energy);
			break;
		}
	}
	int N_elec = 0, allowed = 0;

	vector<RadialWF> Orbitals(orbitals.size());

	for (int i = orbitals.size()-1; i >=0; i--)
	{
		Orbitals[i].set_N(orbitals[i].N());
		Orbitals[i].set_L(orbitals[i].L());
		Orbitals[i].Energy = orbitals[i].Energy;
		Orbitals[i].set_occupancy(orbitals[i].occupancy());
		Orbitals[i].set_infinity(0);
		if (orbitals[i].occupancy() != 0 && allowed == 0) {
			allowed = i;//find highest occupied orbital
		}
		N_elec += orbitals[i].occupancy();
	}
	//	Usually a finer grid is required.
	Grid Lattice(lattice.R(0), 50., 0.03 / k);

	//	Interpolate orbitals on the new grid
	Interpolation W(6);

	//check if there are electrons above current orbital that can fill
	//empty orbital. If above > 1 auger is possible and even hollow orbital should be interpolated
	for (int i = 0; i < orbitals.size(); i++) {
		//if (input.Exited_Pot_Model() != "V_N-1" && orbitals[i].occupancy() == 0) continue;
		if (input.Exited_Pot_Model() == "V_N-1" || i <= allowed) W.RecalcWF(orbitals[i], lattice, Orbitals[i], Lattice);
		Infinity = lattice.R(orbitals[i].pract_infinity());
		int j = 0;
		while (Lattice.R(j) < Infinity && j < Lattice.size() - 1) j++;
		Orbitals[i].set_infinity(j);
	}
	//	potential on the new grid
	Potential U(&Lattice, u.NuclCharge(), u.Type());
	//U.GenerateTrial(Orbitals);

	bool select;
	int h_occ = 0, e_occ = 0, f_occ = 0;

	for (int h = 0; h < allowed; h++)
	{
		if (Orbitals[h].occupancy() < Max_occ[h])//check if there is a hole
		{
			N_h = Max_occ[h] - Orbitals[h].occupancy();
			Tmp.hole = h;
			h_occ = Orbitals[h].occupancy();
			if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[h].set_occupancy(h_occ + 1);
			for (int f = h + 1; f <= allowed; f++)//fill the hole
			{
				Tmp.fill = f;
				f_occ = Orbitals[f].occupancy();
				if (f_occ == 0) continue;
				if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[f].set_occupancy(f_occ - 1);
				for (int e = f; e <= allowed; e++)//eject
				{
					e_occ = Orbitals[e].occupancy();
					if (e_occ == 0) continue;
					if (e == f && input.Exited_Pot_Model() == "V_N-1no" && e_occ < 1) continue;
					if (e == f && input.Exited_Pot_Model() != "V_N-1no" && e_occ < 2) continue;
					if (e == f && input.Exited_Pot_Model() != "V_N-1no") e_occ--;
					double T = 0;
					Tmp.eject = e;
					Tmp.val = 0;
					E_cont = Orbitals[f].Energy + Orbitals[e].Energy - Orbitals[h].Energy;
					//Check if the state reaches asymptotic within the lattice
					if (E_cont < 5 * (u.NuclCharge() - N_elec + 1) /Lattice.R(Lattice.size()-1)) continue;
					/*if (E_cont < 0) continue;
					if (E_cont < 5 * (u.NuclCharge() - N_elec + 1) /Lattice.R(Lattice.size()-1)) {
						E_cont = 5 * (u.NuclCharge() - N_elec + 1) /Lattice.R(Lattice.size()-1);
					}*/
					if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[e].set_occupancy(e_occ - 1);
					Tmp.energy = E_cont;

					if (e == f)	{
						N_fe = 1.* e_occ * f_occ / (4 * Orbitals[e].L() + 2) / (4 * Orbitals[e].L() + 1);
						T = 1 / sqrt(2);
					} else {
						N_fe = 1.* e_occ * f_occ / (4 * Orbitals[e].L() + 2) / (4 * Orbitals[f].L() + 2);
						T = 1;
					}

					N_fe *= Constant::Pi*N_h / (2 * Orbitals[h].L() + 1);
					int min_L_cont;
					if (abs(Orbitals[e].L() - Orbitals[f].L()) <= Orbitals[h].L() && Orbitals[e].L() + Orbitals[f].L() >= Orbitals[h].L()) { min_L_cont = 0; }
					if (abs(Orbitals[e].L() - Orbitals[f].L()) > Orbitals[h].L()) { min_L_cont = abs(Orbitals[e].L() - Orbitals[f].L()) - Orbitals[h].L(); }
					if (Orbitals[e].L() + Orbitals[f].L() < Orbitals[h].L()) { min_L_cont = Orbitals[h].L() - Orbitals[e].L() - Orbitals[f].L(); }
					RadialWF Continuum(Lattice.size());
					Continuum.Energy = E_cont;
					Continuum.set_occupancy(1);
					Continuum.set_N(-1);
					Continuum.set_infinity(Lattice.size() - 1);
					Continuum.set_L(0);
					if (input.Hamiltonian() == 1) U.LDA_upd_dir(Orbitals);
					else U.HF_upd_dir(&Continuum, Orbitals);
					if (input.Exited_Pot_Model() == "V_N-1") U.HF_V_N1(&Continuum, Orbitals, f, true, false);
					for (int l_E = min_L_cont; l_E <= Orbitals[e].L() + Orbitals[f].L() + Orbitals[h].L(); l_E++)
					{
						//sum over all posible Continuum states
						Continuum.set_L(l_E); //TODO make sure this is fine -S.P.
						if (IntegrateContinuum(Lattice, U, Orbitals, &Continuum, f) < 0) {
							log << "Continuum didn't converge: " << endl;
							for (int i = 0; i < orbitals.size(); i++)
							{
								log << i + 1 << ") n = " << orbitals[i].N() << " l = " << orbitals[i].L()
									<< " Energy = " << orbitals[i].Energy << " Occup = " << orbitals[i].occupancy() << endl;
							}
							log.flush();
						}

						double R_k_dir = 0;//direct coulomb integral
						double R_k_exc = 0;//exchange coulomb integral
						int L_min = abs(Orbitals[e].L() - Orbitals[f].L());
						int L_max = Orbitals[e].L() + Orbitals[f].L();
						//ME is an array of Matrix Elements with LS coupling scheme
						//ME[L][S]
						vector<vector<double> > M_LS(L_max - L_min + 1, vector<double>(2, 0));
						for (int K = 0; K <= 2*max(max(Orbitals[f].L(), Orbitals[e].L()), max(Orbitals[h].L(), l_E)); K++)//min(abs(Orbitals[e].L() - Orbitals[f].L()), abs(Orbitals[h].L() - l_E))
						{
							//3j symbol conditions for direct and exchange coulomb integrals
							//direct
							select = (!Triad(Orbitals[h].L(), K, Orbitals[f].L())
								  || (Orbitals[h].L() + K + Orbitals[f].L()) % 2 == 1
								  || !Triad(Orbitals[e].L(), K, Continuum.L())
								  || (Orbitals[e].L() + K + Continuum.L()) % 2 == 1);
							if (select) { R_k_dir = 0; }
							else R_k_dir = U.R_k(K, Orbitals[h], Continuum, Orbitals[f], Orbitals[e]);
							//exchange
							select = (!Triad(Orbitals[h].L(), K, Orbitals[e].L())
								  || (Orbitals[h].L() + K + Orbitals[e].L()) % 2 == 1
								  || !Triad(Orbitals[f].L(), K, Continuum.L())
								  || (Orbitals[f].L() + K + Continuum.L()) % 2 == 1);
							if (select) { R_k_exc = 0; }
							else R_k_exc = U.R_k(K, Orbitals[h], Continuum, Orbitals[e], Orbitals[f]);
							//loop over angular momenta
							if (R_k_dir == 0 && R_k_exc == 0) continue;
							for (int L = L_min; L <= L_max; L++)
							{
								double A_dir = A_k(K, L, Orbitals[h].L(), Continuum.L(), Orbitals[f].L(), Orbitals[e].L());
								double A_exc = A_k(K, L, Orbitals[h].L(), Continuum.L(), Orbitals[e].L(), Orbitals[f].L());
								if (A_dir == 0 && A_exc == 0) continue;
								for (int S = 0; S <= 1; S++)
								{
									M_LS[L-L_min][S] += T*pow(-1, L + Orbitals[f].L() + Continuum.L())*(R_k_dir*A_dir + pow(-1, L + S)*R_k_exc*A_exc);
								}
							}
						}
						for (int L = L_min; L <= L_max; L++)
						{
							for (int S = 0; S <= 1; S++) {
								//if (f == 1 && e == 2) printf("%1.6E  L = %d  S = %d  l-b = %2d\n", M_LS[L - L_min][S], L, S, Continuum.L());
								Tmp.val += N_fe*(2 * L + 1)*(2 * S + 1) * M_LS[L - L_min][S] * M_LS[L - L_min][S];
							}
						}

					}
					if (Tmp.val != 0) Result.push_back(Tmp);
					if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[e].set_occupancy(e_occ);
				}
				if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[f].set_occupancy(f_occ);
			}
			if (input.Exited_Pot_Model() == "V_N-1no") Orbitals[h].set_occupancy(h_occ);
		}
	}

	return Result;
}


vector<double> DecayRates::FT_density(double Q_min, double Q_max, int Q_size)
{
	// Fourier Transform of radially symmetric density. Essentially, an integral:
	// int_0^{infty} dr density(r) sinc(2*pi*Q*r)
	double Q = Q_min;
	double dQ = (Q_max - Q_min)/(Q_size - 1);

	Grid Inp_lattice = lattice;

	vector<double> Output(Q_size, 0);
	vector<double> Inp_density = u.make_density(orbitals);

	// find practical infty of density
	int infty = lattice.size()-1;
	double cutoff = 1e-6, tmp;
	while (Inp_density[infty] < cutoff) infty--;

	// Check if a finer mesh is required. Criteria - 
	// 2*M_PI*Q*dR(infty) = 1/20 - at least 20 points per period.
	Adams IT(Inp_lattice, 6);

	vector<double> integr(infty+1, 0);

	for (int m = 0; m < Q_size; m++) {
		if (2*M_PI*Inp_lattice.dR(infty)*Q > 0.05) {
			// Create a new lattice and do interpolation.
			double dR_max = 0.05/(2*M_PI*Q_max);
			Inp_lattice = Grid(lattice.R(0), lattice.R(infty), dR_max);

			Interpolation IN(6);
			vector<double> density = Inp_density;
			Inp_density.clear();
			Inp_density.resize(Inp_lattice.size(), 0);
			for (int i = 0; i < Inp_lattice.size(); i++) {
				Inp_density[i] = IN.get_value(density, lattice.all_R(), Inp_lattice.R(i))[0];
			}
			integr.clear();
			integr.resize(Inp_lattice.size(), 0);
		}

		if (Q == 0) {
			for (int i = 0; i < integr.size(); i++) {
				integr[i] = Inp_density[i];
			}
		} else {
			for (int i = 0; i < integr.size(); i++) {
				tmp = 2*M_PI*Inp_lattice.R(i)*Q;
				integr[i] = Inp_density[i]*sin(tmp)/tmp;
			}
		}

		Output[m] = IT.Integrate(&integr, 0, integr.size()-1);
		Q += dQ;
	}

	return Output;
}


DecayRates::~DecayRates()
{
}

// Correct continuum integrator. Check if potential accounts for tail correction.

int DecayRates::IntegrateContinuum(Grid &Lattice, Potential &U, vector<RadialWF> &Core, RadialWF* Current, int c)
{
	// This function does outwards integration to find Continuum Wavefunction. It integrates the Hartree-Fock potential U, which is unchanged
	// till point "infinity" where semiclassical momentum P = sqrt(2 * Current->Energy) is withing one percent of the
	// current momentum given by P_new = sqrt(2*Current->Energy - 2*U.V[i] - L()*(L() + 1) / R(i)^2)
	// after that semiclassical approack used. At that point the following equations are solved:
	//
	//  N*sin(P*R(infinity) + phase) = Current->F[infinity]
	//  N*cos(P*R(infinity) + phase) = Current->G[infinity] / P
	//
	// to find Norm (N) and Phase (phase) and continue solution analytically to the end of coordinate grid and normalize it.
	// The functions are Energy normalized.

	double correction, P = sqrt(2 * Current->Energy);
	double accuracy = pow(10, -2), Norm = 10, new_Norm = 5, Phase = 5, new_Phase = 10;
	int infinity = Lattice.size() - 1;// Infinity for continuum wave is defined when it reaches it's asymptotic.
	int Core_infinity = 0;

	// Try to optimize practical infinity.
	for (int i = 0; i < Core.size(); i++) {
		if (Core[i].pract_infinity() > Core_infinity) Core_infinity = Core[i].pract_infinity();
	}

	for (int i = 0; i < Lattice.size(); i++) {
		if (Current->Energy > -5 * U.V[i]) {
			infinity = i;
			break;
		}
	}

	if (Core_infinity > infinity) Current->set_infinity(Core_infinity);
	else Current->set_infinity(infinity);

	Adams NumIntgr(Lattice, 10);

	for (int i = 0; i < Lattice.size(); i++)
	{
		NumIntgr.B[i] = 1.;
		NumIntgr.C[i] = -2 * (Current->Energy - U.V[i] - 0.5*Current->L()*(Current->L() + 1) / Lattice.R(i) / Lattice.R(i));
		U.Exchange[i] = 0;
	}

	// Loop until exchange doesn't shift the phase between iterations.
	int m = 0;
	while (fabs(new_Phase / Phase - 1) > accuracy || fabs(new_Norm / Norm - 1) > accuracy || m < 2)
	{
		if (m > 20) break;
		m++;
		Phase = new_Phase;
		Norm = new_Norm;

		correction = Current->F[0];
		if (correction)	{
			Current->F[0] = correction;
			Current->G[0] = Current->F[0] * (Current->L() + 1 + U.V[0] * Lattice.R(0)*Lattice.R(0) / (Current->L() + 1)) / Lattice.R(0);
		} else {
			Current->F[0] = pow(Lattice.R(0), (Current->L() + 1));
			Current->G[0] = pow(Lattice.R(0), Current->L())*(Current->L() + 1 + U.V[0] * Lattice.R(0)*Lattice.R(0) / (Current->L() + 1));
		}

		for (int i = 0; i <= Current->pract_infinity(); i++) {
			NumIntgr.Y[i] = 2 * U.Exchange[i];
		}

		NumIntgr.StartAdams(Current, 0, true);
		NumIntgr.Integrate(Current, 0, Current->pract_infinity());
		infinity = Current->pract_infinity();
		while (Current->G[infinity - 1] * Current->G[infinity] > 0) { infinity--; }

		new_Phase = atan(P*Current->F[infinity] / Current->G[infinity]);
		new_Norm = 0.5 * fabs(Current->F[infinity])*sqrt(2*Pi*P);

		Current->scale(1 / new_Norm);
		if (input.Hamiltonian() == 0) U.HF_upd_exc(Current, Core);
		if (input.Exited_Pot_Model() == "V_N-1") U.HF_V_N1(Current, Core, c, false, true);
		if (U.Exchange[0] == 0) break;
	}
	if (m > 20) return -1;

	//Current->set_infinity(Lattice.size() - 1);
	return infinity;
}


double A_k(int k, int L, int la, int lb, int lc, int ld)
{
	//calculates <la||C_k||lc><lb||C_k||ld> [ la lb L ]
	// 									    [ ld lc k ]

	double Result;
	Result = Constant::Wigner3j(lc, k, la, 0, 0, 0);
	if (Result == 0) return Result;
	Result *= Constant::Wigner3j(ld, k, lb, 0, 0, 0);
	if (Result == 0) return Result;
	Result *= WignerSymbols::wigner6j(la, lb, L, ld, lc, k);
	if ((la + lb) % 2 == 1) Result = -Result;
	Result *= sqrt((2 * la + 1)*(2 * lb + 1)*(2 * lc + 1)*(2 * ld + 1));
	return Result;
}


bool Triad(int l_a, int l_b, int l_c)
{
	bool Result = true;
	if (l_a > l_b + l_c) { Result = false; }
	if (l_a < abs(l_b - l_c)) { Result = false; }
	if (l_b > l_a + l_c) { Result = false; }
	if (l_b < abs(l_a - l_c)) { Result = false; }
	if (l_c > l_b + l_a) { Result = false; }
	if (l_c < abs(l_b - l_a)) { Result = false; }

	return Result;
}
