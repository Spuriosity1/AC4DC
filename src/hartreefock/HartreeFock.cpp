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
#include "HartreeFock.hpp"
#include <algorithm>

using namespace std;

int SetBoundaryValues(Grid*, RadialWF*, Potential*);
int SetBoundaryValuesApprox(Grid*, RadialWF*, Potential*);

HartreeFock::HartreeFock(Grid &Lattice, std::vector<RadialWF> &Orbitals, Potential &Potential, const HFInput& Inp, std::ofstream & log) : lattice(&Lattice)
{
//==========================================================================================================
// Estimate starting energies using Slater rules.
	float s = 0, p = 1;
	int N_elec_shell = 0;// Current shell.
	int N_elec_n1 = 0, N_elec_n2 = 0;// Shell with n-1 and n-2 respectively.

	Potential.GenerateTrial(Orbitals);

	for (int i = 0; i < Orbitals.size(); i++)
	{
		N_elec_shell = Orbitals[i].occupancy() - 1;
		N_elec_n1 = 0;
		N_elec_n2 = 0;
		for (int j = 0; j < Orbitals.size(); j++)
		{
			if (j == i) { continue; }// No double counting.
			if (Orbitals[i].L() <= 1)// [s] and [s,p] groups.
			{
				if (Orbitals[j].N() == Orbitals[i].N() && Orbitals[j].L() < 2) { N_elec_shell += Orbitals[j].occupancy(); }
				if (Orbitals[j].N() == Orbitals[i].N() - 1) { N_elec_n1 += Orbitals[j].occupancy(); }
				if (Orbitals[j].N() < Orbitals[i].N() - 1) { N_elec_n2 += Orbitals[j].occupancy(); }
			}
			else if (Orbitals[i].L() > Orbitals[j].L() && Orbitals[j].N() <= Orbitals[i].N())// [d] and [f] groups.
			{
				N_elec_n2 += Orbitals[j].occupancy();
			}
		}
		s = 0.35*N_elec_shell + 0.85*N_elec_n1 + 1.0*N_elec_n2;
		Orbitals[i].Energy = -0.5*(Potential.NuclCharge() - s)*(Potential.NuclCharge() - s) / Orbitals[i].N() / Orbitals[i].N();
	}
//==========================================================================================================
// Find initial Guess for wavefunctions. Check if there is a single orbital occupied. If there is
// only one electron (hydrogenic case) solve problem.

	double Norm = 0.;
	std::vector<double> E_rel_change(Orbitals.size(), 1);
	double duration;
	int infinity = 0;

	Adams I(Lattice, 10);
	int m = 0;
	double E_max_error = 0;

	E_max_error = 1;

// Check if the system has a single occupied orbital
	int check_orb = 0, single = -1;
	for (int i = 0; i < Orbitals.size(); i++)
	{
		if (Orbitals[i].occupancy() != 0)
		{
			check_orb++;
			single = i;
		}
	}
	if (check_orb != 1) single =-1;

	Potential.GenerateTrial(Orbitals);

	if (check_orb == 1 && Orbitals[single].occupancy() == 1) Potential.Reset();

	for (int i = 0; i < Orbitals.size(); i++) {
		Master(&Lattice, &Orbitals[i], &Potential, Inp.master_tolerance, log);
	}

//==========================================================================================================
	// Solve a single shell problem. Single shell Hartree-Fock/Hartree-fock-Slater.
	std::vector<RadialWF> Orbitals_old = Orbitals;
	std::vector<double> V_old = Potential.V;
	double V_tmp = 0;

	if (check_orb == 1 && Orbitals[single].occupancy() > 1 && Inp.ham_model == HFInput::ham_mod_t::HF) {
		// A single occupied orbital. HF (with Exchange) solution is here.
		// Afterwards unoccupied orbitals are calculated.
		for (int i = 0; i < Orbitals.size(); i++) {
			if (Orbitals[i].occupancy() != 0) single = i;
		}

		while (E_rel_change[single] > Inp.HF_tolerance) {
			if (m > 20) {
				log << "Starting approximation does not converge... " << endl;
				break;
			}

			Potential.HF_upd_dir(&Orbitals[single], Orbitals);
			for (int j = 0; j < Lattice.size(); j++) {
				Potential.V[j] = p*Potential.V[j] + (1 - p)*V_old[j];
			}
			if (m == 0) p = 0.5;
			if (m == 6) p = 0.8;

			Master(&Lattice, &Orbitals[single], &Potential, Inp.master_tolerance, log);
			E_rel_change[single] = fabs(Orbitals[single].Energy / Orbitals_old[single].Energy - 1);

			Orbitals_old[single] = Orbitals[single];
			V_old = Potential.V;
			m++;
		}
		m = 0;
	}
	if (check_orb == 1 && Orbitals[single].occupancy() == 1) {
		Potential.Reset();
		for (auto& Orb: Orbitals) Master(&Lattice, &Orb, &Potential, Inp.master_tolerance, log);
	} else {
		p = 0.5;
		// There is more than one orbital.

		// Set up starting approximation for the Hartree-Fock equations.
		// Algorithm follows W. Johnson, but with Local Exchange.

		double LDA_tolerance = Inp.no_exch_tolerance;
		if (Inp.ham_model == HFInput::ham_mod_t::LDA) LDA_tolerance = Inp.HF_tolerance;
		else LDA_tolerance = Inp.no_exch_tolerance;
		bool Final_Check = false;
		Orbitals_old = Orbitals;

		//=======================================================================================
		// Hartree-Fock loop without exchange. Local exchange is evaluated here.
		while (E_max_error > LDA_tolerance || m < 2)
		{
			if (m > 20 * Orbitals.size()) {
				log << "Starting approximation does not converge... " << endl;
				log.flush();
				break;
			}

			Potential.LDA_upd_dir(Orbitals);

			// Smooth start for m = 0. Introduce latter tail correction. May help iof the original routine diverges.
			/*
			if (m == 0)	{
				for (int j = 0; j < Lattice.size(); j++) {
					Potential.V[j] = p*Potential.V[j] + (1 - p)*Potential.Trial[j];
				}
			} else {
				for (int j = 0; j < Lattice.size(); j++) {
					Potential.V[j] = p*Potential.V[j] + (1 - p) * V_old[j];
				}
			}*/

			for (int i = 0; i < Orbitals.size(); i++)
			{
				/*if (i == single) {
					E_rel_change[i] = 0;
					continue;
				}*/
				if (E_rel_change[i] < E_max_error && m != 0) continue;

				if (Master(&Lattice, &Orbitals[i], &Potential, Inp.master_tolerance, log)) {
					// Master didn't converge. This is bad. Return to old solution and
					// iterated second worst instead.
					for (auto& orb: Orbitals) log << orb.occupancy() << " ";
					log << endl;
					Orbitals[i] = Orbitals_old[i];
					E_max_error = 0;
					for (int j = 0; j < Orbitals.size(); j++) {
						if (E_rel_change[j] > E_max_error && j != i) E_max_error = E_rel_change[j];
					}
					Orbitals[i].Energy *= 1 + 0.99*E_max_error;
				}
				else {
					E_rel_change[i] = fabs(Orbitals[i].Energy / Orbitals_old[i].Energy - 1);
					MixOldNew(&Orbitals[i], &Orbitals_old[i]);
					Orbitals_old[i] = Orbitals[i];
				}
			}

			V_old = Potential.V;

			// Find new largest relative error.
			E_max_error = 0;
			for (int i = 0; i < Orbitals.size(); i++) {
				if (Orbitals[i].Energy != Orbitals_old[i].Energy) {
					E_rel_change[i] = fabs(Orbitals[i].Energy / Orbitals_old[i].Energy - 1);
				}
				if (E_max_error < E_rel_change[i]) {
					E_max_error = E_rel_change[i];
				}
				Orbitals_old[i].Energy = Orbitals[i].Energy;
			}
			m++;
			// Final check. Varry all orbitals at once.
			if ( !(E_max_error > LDA_tolerance || m < 2) && !Final_Check) {
				Final_Check = true;
				E_max_error = 2*LDA_tolerance;
				for (std::vector<double>::iterator it = E_rel_change.begin(); it != E_rel_change.end(); ++it) *it = 1;
			}
			else Final_Check = false;

		}

		//=============================================================================================
		if (Inp.ham_model == HFInput::ham_mod_t::HF) {
		// Hartree-Fock method. Add Exhange potential.

			std::vector<std::vector<double>> Exchange_old(Orbitals.size(), std::vector<double>(Lattice.size(), 0));
			std::vector<std::vector<double>> Direct_old(Orbitals.size(), std::vector<double>(Lattice.size(), 0));
			Orbitals_old = Orbitals;

			double correction_scaling = 1, change_cs = 1;
			E_max_error = 1;
			m = 0;
			double max_norm_dev = 1;

			// Hartree-Fock loop with exchange
			while (E_max_error > Inp.HF_tolerance || m < 1)
			{
				if (m > Inp.max_HF_iterations) { break; }

				for (int i = 0; i < Orbitals.size(); i++)
				{
					if (i == single) continue;
					if (E_rel_change[i] < E_max_error && m != 0) continue;

					Potential.HF_upd_dir(&Orbitals[i], Orbitals_old);
					Potential.HF_upd_exc(&Orbitals[i], Orbitals_old);

					if (E_max_error < Inp.HF_tolerance * 100) p = 0.8;
					else p = 0.5;
					if (m == 0) p = 1;

					for (int j = 0; j < Lattice.size(); j++) {
						Potential.V[j] = p*Potential.V[j] + (1 - p)*Direct_old[i][j];
						Potential.Exchange[j] = p*Potential.Exchange[j] + (1 - p)*Exchange_old[i][j];
					}

					if (m != 0) {
						Master(&Lattice, &Orbitals[i], &Potential, Inp.master_tolerance, log);
						double numerator = 0, denominator = 0;
						for (int j = 0; j < max(Orbitals[i].pract_infinity(), Orbitals_old[i].pract_infinity()); j++) {
							numerator += Potential.Exchange[j] * Orbitals[i].F[j] * Lattice.dR(j);
							denominator += Orbitals[i].F[j] * Orbitals_old[i].F[j] * Lattice.dR(j);
						}
						Orbitals[i].Energy += numerator / denominator;
					}

					int max_iter = 0;
					correction_scaling = 1;
					do
					{
						max_iter++;
						if (max_iter > 20)
						{
							int new_max = 0;
							for (int j = 0; j < E_rel_change.size(); j++)
							{
								if (j == i) continue;
								if (E_rel_change[new_max] < E_rel_change[j]) new_max = j;
							}
							Orbitals[i].Energy = Orbitals_old[i].Energy*(1+E_rel_change[i]);
							Orbitals[i].F = Orbitals_old[i].F;
							if (E_rel_change[new_max] > Inp.HF_tolerance) break;
						}
						Orbitals[i].Energy *= correction_scaling;
						GreensMethod P(&Lattice, &Orbitals[i], &Potential);
						if (Orbitals[i].check_nodes() == Orbitals[i].GetNodes()) {
							correction_scaling = 1;
							E_rel_change[i] = fabs(Orbitals[i].Energy / Orbitals_old[i].Energy - 1);
							Exchange_old[i] = Potential.Exchange;
							Direct_old[i] = Potential.V;
							if (m != 0) {
								Orbitals_old[i] = Orbitals[i];
								break;
							}
						}
						else {// New wavefunction have incorrect n. Iterate untill it is correct.
							change_cs = 0.2*(Orbitals[i].check_nodes() - Orbitals[i].GetNodes()) / Orbitals[i].N();
							if (fabs(change_cs) <= 0.2) correction_scaling = 1 + change_cs;
							else correction_scaling = 1 + 0.2*change_cs/fabs(change_cs);
						}
					} while (correction_scaling != 1);
				}
				E_max_error = 0;
				for (int i = 0; i < Orbitals.size(); i++) {
					// Find new largest relative error.
					/*if (Orbitals[i].Energy != Orbitals_old[i].Energy) {
						E_rel_change[i] = fabs(Orbitals[i].Energy / Orbitals_old[i].Energy - 1);
					}*/
					if (E_max_error < E_rel_change[i]) {
						E_max_error = E_rel_change[i];
					}
				}
				m++;
				if (E_max_error < Inp.HF_tolerance && !Final_Check)
				{
					for (int i = 0; i < Orbitals.size(); i++)
					{
						E_rel_change[i] = 1;
					}
					E_max_error = 2 * Inp.HF_tolerance;
					Final_Check = true;
				}
			}

		}
	}

	if (m >= Inp.max_HF_iterations && log.is_open()) {
		log << "====================================================================" << endl;
		log << "Too many iterations in Hartree-Fock" << endl;
		for (int i = 0; i < Orbitals.size(); i++) {
			log << i + 1 << ") n = " << Orbitals[i].N() << " l = " << Orbitals[i].L()
				<< " Energy = " << Orbitals[i].Energy << " Occup = " << Orbitals[i].occupancy() << endl;
		}
		log.flush();
	}
}

double HartreeFock::OrthogonalityTest(std::vector<RadialWF> &Orbitals)
{
	Adams I(*lattice, 10);

	double Result = 0;
	double ort = 0;
	std::vector<double> density(lattice->size(), 0);

	for (int i = 0; i < Orbitals.size(); i++)
	{
		for (int j = i + 1; j < Orbitals.size(); j++)
		{
			if (Orbitals[i].L() != Orbitals[j].L()) continue;
			int max_R = max(Orbitals[i].pract_infinity(), Orbitals[j].pract_infinity());
			for (int k = 0; k <= max_R; k++)
			{
				density[k] = Orbitals[i].F[k] * Orbitals[j].F[k];
			}
			ort = fabs(I.Integrate(&density, 0, max_R));
			if (Result < ort) Result = ort;
		}
	}

	return Result;
}

HartreeFock::~HartreeFock()
{
}

int SetBoundaryValues(Grid* Lattice, RadialWF* Psi, Potential* U)
{
	//Near the origin the behaviour is determined by the centrifugal term l(l+1)/r^2, therefore asymptotic is same for Coulomb and finite size potentials.
	int Correction = 1;
	double exchange_correction = pow(10, -10);
	if (Psi->F[0] < 0) { Correction = -1; }
	Psi->F[0] = Correction*pow(Lattice->R(0), (Psi->L() + 1));
	Psi->G[0] = Correction*pow(Lattice->R(0), Psi->L())*(Psi->L() + 1 + U->V[0] * Lattice->R(0)*Lattice->R(0) / (Psi->L() + 1));

	unsigned int Turn = 1;
	unsigned int infinity = Lattice->size() - 1;
	double S = 0, P = 0;

	if (Psi->Energy < 0.)
	{

		while (Psi->Energy > U->V[Turn] && Turn < Lattice->size() - 1)
		{
			Turn++;
		}

		//set boundary conditions for outwards integration
		infinity = Turn;
		Psi->set_turn(Turn);
		while (S < 22. || U->Exchange[infinity] > exchange_correction)
		{
			S += sqrt(U->V[infinity] - Psi->Energy)*Lattice->dR(infinity);
			infinity++;
			if (infinity >= (Lattice->size() - 1)) { break; }
		}

		//set boundary values for inwards integration
		P = sqrt(2 * (U->V[infinity] - Psi->Energy));
		Psi->F[infinity] = pow(10., -14);
		Psi->G[infinity] = -Psi->F[infinity] * P;
	}

	Psi->set_infinity(infinity);

	return infinity;
}

int SetBoundaryValuesApprox(Grid * Lattice, RadialWF * Psi, Potential* U)
{
	//Create first 10 (max_adams_size) points for integration using asymptotics.
	int Correction = 1;
	double exchange_correction = pow(10, -8);
	if (Psi->F[0] < 0) { Correction = -1; }
	Psi->F[0] = Correction*pow(Lattice->R(0), (Psi->L() + 1));
	Psi->G[0] = Correction*pow(Lattice->R(0), Psi->L())*(Psi->L() + 1 + U->V[0] * Lattice->R(0)*Lattice->R(0) / (Psi->L() + 1));

	int Turn = 1;
	int infinity = Lattice->size() - 1;
	int adams_max_order = 10;

	double sigma, lambda;
	int order = 3;
	std::vector<double> a(order, 0);
	std::vector<double> b(order, 0);
	double S = 0;

	if (Psi->Energy < 0.)
	{

		while (Psi->Energy > U->V[Turn] && Turn < Lattice->size() - 1)
		{
			Turn++;
		}

		//set boundary conditions for outwards integration
		infinity = Turn;
		Psi->set_turn(Turn);
		while (S < 17. || U->Exchange[infinity] > exchange_correction)
		{
			S += sqrt(U->V[infinity] - Psi->Energy)*Lattice->dR(infinity);
			if (infinity >= (Lattice->size() - 1)) { break; }
			infinity++;
		}

		//set boundary values for inwards integration
		lambda = sqrt(-2 * Psi->Energy);
		sigma = (U->V[infinity] * Lattice->R(infinity) - 1.) / lambda;
		a[0] = 1;
		b[0] = -lambda;
		for (int i = 1; i < a.size(); i++)
		{
			a[i] = a[i - 1] * (Psi->L()*(Psi->L() + 1) - (sigma - i)*(sigma - i + 1)) / (2 * i*lambda);
			b[i] = a[i - 1] * ((sigma + i)*(sigma - i + 1) - Psi->L()*(Psi->L() + 1)) / (2 * i);
		}

		for (int Inf = infinity; Inf > infinity - adams_max_order; Inf--)
		{
			for (int k = order - 1; k > 0; k--)
			{
				Psi->F[Inf] += a[k];
				Psi->G[Inf] += b[k];
				Psi->F[Inf] /= Lattice->R(Inf);
				Psi->G[Inf] /= Lattice->R(Inf);
			}
			Psi->F[Inf] += a[0];
			Psi->G[Inf] += b[0];
			S = pow(Lattice->R(Inf), sigma)*exp(-lambda*Lattice->R(Inf));
			Psi->F[Inf] *= S;
			Psi->G[Inf] *= S;

		}
	}

	Psi->set_infinity(infinity);

	return infinity;
}

int HartreeFock::Master(Grid* Lattice, RadialWF* Psi, Potential* U, double Epsilon, std::ofstream & log)
{
	// This routine does the same as "Master" by W. Johnson.
	// 1) Takes Psi.Energy_0 and integrates the HF equations inwards and outwards
	// 2) On each iteration the energy is corrected to the point when |Psi.Energy(i) - Psi.Energy(i-1)| < Epsilon
	// 3) Normalizes Psi in the end of the routine
	double F_left, G_left = 1., F_right, G_right = -1., E_tmp = Psi->Energy, Norm = 0.0, old_Energy;
	double E_low = -0.5*U->NuclCharge()*U->NuclCharge() / Psi->N() / Psi->N(), E_high = 0;
	int Turn = 1, infinity = 0, NumNodes = 0.;
	std::vector<double> density;

	int Alarm = 0;
	Adams NumIntgr(*Lattice, 10);

	while (Psi->Energy > U->V[Turn] && Turn < Lattice->size() - 1) Turn++;


	for (int i = 0; i < Lattice->size(); i++)
	{
		NumIntgr.B[i] = 1.;
		NumIntgr.C[i] = -2 * (Psi->Energy - U->V[i] - 0.5*Psi->L()*(Psi->L() + 1) / Lattice->R(i) / Lattice->R(i));
		Psi->F[i] = 0.;
		Psi->G[i] = 0.;
	}

	while (fabs(E_tmp / Psi->Energy) > Epsilon)
	{
		Alarm++;
		if (Alarm > 25) {
			log << "Oops, Master failed to converge." << endl
				<< "Practical infinity: " << Lattice->R(infinity) << endl
				<< "n = " << Psi->N() << " l = " << Psi->L() << " Energy = " << Psi->Energy << endl;
			return 1;
		}

		//Find the point where in and out integrated functions are to be stiched. The point is the largest maximum before classical Turning point.
		//Should be done prior to the outwards integration.

		infinity = SetBoundaryValuesApprox(Lattice, Psi, U);//Approx
		Turn = Psi->turn_pt();
		if (Psi->pract_infinity() == Lattice->size() - 1) {
			E_tmp = -2 * Psi->Energy;
			for (int i = 0; i < Lattice->size(); i++) {
				NumIntgr.C[i] += E_tmp;
			}
			Psi->Energy *= 2;
			continue;
		}

		NumIntgr.StartAdams(Psi, 0, true);

		density.clear();
		density.resize(infinity + 1);

		NumIntgr.Integrate(Psi, 0, Turn);

		Norm = 0.0;

		F_left = Psi->F[Turn];
		G_left = Psi->G[Turn] / F_left;

		NumIntgr.StartAdams(Psi, infinity, false);
		NumIntgr.Integrate(Psi, infinity, Turn);

		F_right = Psi->F[Turn];
		G_right = Psi->G[Turn] / F_right;

		for (int i = 0; i <= infinity; i++) {
			if (i >= Turn) {
				Psi->F[i] /= F_right;
				Psi->G[i] /= F_right;
			} else {
				Psi->F[i] /= F_left;
				Psi->G[i] /= F_left;
			}

			density[i] = Psi->F[i] * Psi->F[i];
		}

		Norm = NumIntgr.Integrate(&density, 0, infinity) + Psi->F[0] * Psi->F[0] * Lattice->R(0) / (2 * Psi->L() + 3);//last term accounts for WF between 0 and Lattice.R(0)
		E_tmp = 0.5*((G_right - G_left) / Norm);
		NumNodes = Psi->check_nodes();

		old_Energy = Psi->Energy;
		if (NumNodes == Psi->GetNodes()) {
			if (Psi->Energy - E_tmp > E_high) {
				if (E_low < 1.5*Psi->Energy) E_low = 1.5*Psi->Energy;
				Psi->Energy = 0.5*Psi->Energy + 0.5*E_high;
			}
			else if (Psi->Energy - E_tmp < E_low) {
				if (E_high > 0.75*Psi->Energy) E_high = 0.75*Psi->Energy;
				Psi->Energy = 0.5*Psi->Energy + 0.5*E_low;
			}
			else {
				if (E_tmp > 0) E_high = 0.5*E_high + 0.5*Psi->Energy;
				else E_low = 0.5*E_low + 0.5*Psi->Energy;
				if (fabs(E_tmp) > 0.1*fabs(Psi->Energy)) E_tmp = 0.1*fabs(Psi->Energy*E_tmp)/E_tmp;
				Psi->Energy -= E_tmp;
			}
		} else {
			if (NumNodes < Psi->GetNodes()) {
				if (Psi->Energy > E_low) E_low = Psi->Energy;
			} else {
				if (Psi->Energy < E_high) E_high = Psi->Energy;
			}
			Psi->Energy = 0.5*(E_low + E_high);
		}
		/*
		if (NumNodes == Psi->GetNodes()) {
			if (Psi->Energy - E_tmp > E_high) { Psi->Energy = 0.5*Psi->Energy + 0.5*E_high; }
			else if (Psi->Energy - E_tmp < E_low) { Psi->Energy = 0.5*Psi->Energy + 0.5*E_low; }
			else { Psi->Energy -= E_tmp; }
		} else {
			if (NumNodes < Psi->GetNodes())	{
				if (Psi->Energy > E_low) { E_low = Psi->Energy; }
			} else {
				if (Psi->Energy < E_high) { E_high = Psi->Energy; }
			}
			Psi->Energy = 0.5*(E_low + E_high);
		}
		*/
		E_tmp = -2 * (Psi->Energy - old_Energy);
		for (int i = 0; i < Lattice->size(); i++) {
			NumIntgr.C[i] += E_tmp;
		}
		E_tmp *= 0.5;

	}

	Psi->set_infinity(infinity);
	//normalize the answer
	Psi->scale(1. / sqrt(Norm));

	return 0;
}

GreensMethod::GreensMethod(Grid* Lattice, RadialWF* Psi, Potential* U) : Adams(*Lattice, 10), lattice(Lattice), psi(Psi), u(U)
{
	int infinity, track_sign, NumNodes = -1;
	double dEnergy = Psi->Energy, E_tolerance = pow(10, -8), Norm_tolerance = pow(10, -10);
	double W, Norm = 0, correct = 0;
	std::vector<double> Integrand(Lattice->size(), 0);

	RadialWF Psi_O = *Psi;//regular at the origin
	RadialWF Psi_Inf = *Psi;//regular at infinity
	std::vector<double> Green_O;
	std::vector<double> Green_Inf;

	if (Psi->F[0] < 0) { track_sign = -1; }
	else { track_sign = 1; }

	for (int i = 0; i < Lattice->size(); i++)
	{
		B[i] = 1.;
		C[i] = -2 * (Psi->Energy - U->V[i] - 0.5*Psi->L()*(Psi->L() + 1) / Lattice->R(i) / Lattice->R(i));
		Y[i] = 0.;
	}

	infinity = SetBoundaryValuesApprox(Lattice, &Psi_O, U);//Approx
	StartAdams(&Psi_O, 0, true);
	Integrate(&Psi_O, 0, infinity);

	infinity = SetBoundaryValuesApprox(Lattice, &Psi_Inf, U);
	Integrate(&Psi_Inf, infinity, 0);

	Integrand.clear();
	Integrand.resize(Lattice->size());

	for (int i = 0; i <= infinity; i++)
	{
		W = Psi_O.F[i] * Psi_Inf.G[i] - Psi_Inf.F[i] * Psi_O.G[i];
		Y[i] = -2.*U->Exchange[i] / W;
	}

	Green_O = GreenOrigin(&Psi_O);
	Green_Inf = GreenInfinity(&Psi_Inf);

	for (int i = 0; i <= infinity; i++)
	{
		Psi->F[i] = Psi_Inf.F[i] * Green_O[i] + Psi_O.F[i] * Green_Inf[i];
		Psi->G[i] = Psi_Inf.G[i] * Green_O[i] + Psi_O.G[i] * Green_Inf[i];
	}

	double hNorm = 0, lNorm = 0;
	double dEnergy_old = dEnergy, HydroEnergy = -0.5*U->NuclCharge()*U->NuclCharge() / Psi->N() / Psi->N();
	double hEnergy = 0, lEnergy = 0;
	//by this point we have unnormalized solution.
	//refining the solution following Johnson's variation of parameter
	while (fabs(dEnergy) > E_tolerance)
	{
		//by this point we have unnormalized solution.
		//refining the solution following Johnson's variation of parameter
		Psi->set_infinity(infinity);
		if (Psi->check_nodes() != Psi->GetNodes()) break;

		for (int i = 0; i <= infinity; i++)
		{
			W = Psi_O.F[i] * Psi_Inf.G[i] - Psi_Inf.F[i] * Psi_O.G[i];
			Y[i] = Psi->F[i] / W;
			Integrand[i] = Psi->F[i] * Psi->F[i];
		}
		Norm = Integrate(&Integrand, 0, infinity);

		Green_O = GreenOrigin(&Psi_O);
		Green_Inf = GreenInfinity(&Psi_Inf);

		for (int i = 0; i <= infinity; i++)
		{
			Integrand[i] = (Psi_Inf.F[i] * Green_O[i] + Psi_O.F[i] * Green_Inf[i])*Psi->F[i];
		}

		correct = Integrate(&Integrand, 0, infinity);

		dEnergy_old = dEnergy;
		dEnergy = 0.25*(Norm - 1.) / correct;
		if (fabs(dEnergy) > 0.2*fabs(Psi->Energy)) dEnergy = 0.2*dEnergy / fabs(dEnergy)*fabs(Psi->Energy);
		if (dEnergy*dEnergy_old < 0 && dEnergy_old != 1 && fabs(dEnergy) > fabs(dEnergy_old)) dEnergy = -0.5*dEnergy_old;

		if (Norm > 1 && hNorm == 0)
		{
			hEnergy = Psi->Energy;
			hNorm = Norm;
		}

		if (Norm < 1 && lNorm == 0)
		{
			lEnergy = Psi->Energy;
			lNorm = Norm;
		}

		if ( hNorm != 0 && lNorm != 0 && (fabs(dEnergy/Psi->Energy) > 0.05
			|| (Psi->Energy + dEnergy < HydroEnergy)) )
		{
			if (Norm > 1)
			{
				if (hEnergy > Psi->Energy) hEnergy = Psi->Energy;
				dEnergy = 0.5*(lEnergy - Psi->Energy);
				Psi->Energy = 0.5*(Psi->Energy + lEnergy);
			}
			else
			{
				if (lEnergy < Psi->Energy) lEnergy = Psi->Energy;
				dEnergy = 0.5*(hEnergy - Psi->Energy);
				Psi->Energy = 0.5*(Psi->Energy + hEnergy);
			}
		}
		else Psi->Energy += dEnergy;

		if (Psi->Energy < HydroEnergy) Psi->Energy = HydroEnergy;

		for (int i = 0; i < lattice->size(); i++)
		{
			if (i <= Psi->pract_infinity()) {
				Psi->F[i] -= 2 * dEnergy*(Psi_Inf.F[i] * Green_O[i] + Psi_O.F[i] * Green_Inf[i]);
				Psi->G[i] -= 2 * dEnergy*(Psi_Inf.G[i] * Green_O[i] + Psi_O.G[i] * Green_Inf[i]);
			} else {
				Psi->F[i] = 0.;
				Psi->G[i] = 0.;
			}
		}

		//use the last bit again
		if (fabs(dEnergy) > E_tolerance)
		{
			for (int i = 0; i < Lattice->size(); i++)
			{
				B[i] = 1.;
				C[i] -= 2 * dEnergy;
				Y[i] = 0.;
				Psi_O.F[i] = 0;
				Psi_O.G[i] = 0;
				Psi_Inf.G[i] = 0;
				Psi_Inf.F[i] = 0;
			}

			infinity = SetBoundaryValuesApprox(Lattice, &Psi_O, U);//Approx
			StartAdams(&Psi_O, 0, true);
			Integrate(&Psi_O, 0, infinity);

			infinity = SetBoundaryValuesApprox(Lattice, &Psi_Inf, U);
			//	I.StartAdams(&Psi_Inf, infinity, false);
			Integrate(&Psi_Inf, infinity, 0);

			Integrand.clear();
			Integrand.resize(infinity + 1);

			Psi_O.Energy = Psi->Energy;
			Psi_Inf.Energy = Psi->Energy;
		}
	}
	if (Psi->F[0] * track_sign < 0) { Psi->scale(-1); }
}

std::vector<double> GreensMethod::GreenOrigin(RadialWF * Psi_O)
{
	//this function returns the following integral
	//int_(0)^(Lattice.R(end_pt)) Psi.F*Y*Lattice.dR,
	// where Y includes the Wronskian.
	std::vector<double> Result;
	double Func_tmp;

	Result.clear();
	Result.resize(lattice->size());

	Result[0] = 0.5*lattice->dR(0) * Psi_O->F[0] * Y[0];//trapezoid rule for first interval [0...Lattice.dR(0)]

	for (int i = 1; i < Adams_N; i++)
	{
		Result[i] = Result[i - 1] + 0.5*Lattice.dR(i) * Psi_O->F[i] * Y[i] +
			0.5*Lattice.dR(i - 1) * Psi_O->F[i - 1] * Y[i - 1];
	}

	for (int i = Adams_N; i <= Psi_O->pract_infinity(); i++)
	{
		Func_tmp = Result[i - 1] + Adams_Coeff[0] * Lattice.dR(i) * Psi_O->F[i] * Y[i];

		for (int j = 1; j < Adams_N; j++)
		{
			Func_tmp += Adams_Coeff[j] * Psi_O->F[i - j] * Y[i - j] * Lattice.dR(i - j);
		}

		Result[i] = Func_tmp;
	}

	return Result;
}

std::vector<double> GreensMethod::GreenInfinity(RadialWF * Psi_Inf)
{
	//this function returns the following integral
	//int_(Lattice.R(start_pt))^(Lattice.R(end_pt)) Psi.F*Y*Lattice.dR,
	// where Y is inhomogenous term in the original ODE. 1/W multiplier is ignored, since it is constant and resultant function will be normalized. It should be constant, but just in case it is recalculated in each grid point

	std::vector<double> Result(Lattice.size(), 0.);
	double Func_tmp;
	int Infty = Psi_Inf->pract_infinity();

	Result[Infty] = 0.5*Lattice.dR(Infty) * Psi_Inf->F[Infty] * Y[Infty];

	for (int i = Infty - 1; i > (Infty - Adams_N); i--)
	{
		Result[i] = Result[i + 1] + 0.5*Lattice.dR(i) * Psi_Inf->F[i] * Y[i] +
			0.5*Lattice.dR(i + 1) * Psi_Inf->F[i + 1] * Y[i + 1];
	}

	for (int i = (Infty - Adams_N); i >= 0; i--)
	{
		Func_tmp = Result[i + 1] + Adams_Coeff[0] * Lattice.dR(i) * Psi_Inf->F[i] * Y[i];

		for (int j = 1; j < Adams_N; j++)
		{
			Func_tmp += Adams_Coeff[j] * Psi_Inf->F[i + j] * Y[i + j] * Lattice.dR(i + j);
		}

		Result[i] = Func_tmp;
	}

	return Result;
}

double HartreeFock::Conf_En(std::vector<RadialWF> &Orbitals, Potential &U)
{
	// Calculate total energy of electronic configuration.
	// Average over configurations is assumed.
	double Result = 0;
	for (auto& Orb: Orbitals) Result += Orb.occupancy()*Orb.Energy;
	// Get orbital Coulomb interaction energies.
	double wght = 0; // Complete shell occupancy.
	double eAB = 0;
	double angular = 0;
	for (int a = 0; a < Orbitals.size(); a++) {
		for (int b = a; b < Orbitals.size(); b++) {
			if (b != a) wght = 1.*Orbitals[b].occupancy();
			else wght = 0.5*(4*Orbitals[b].L() + 2)*(Orbitals[b].occupancy() - 1)/(4*Orbitals[b].L() + 1);
			if (wght == 0) continue;

			eAB = U.R_k(0, Orbitals[a], Orbitals[b], Orbitals[a], Orbitals[b]);
			for (int k = 0; k <= Orbitals[a].L() + Orbitals[b].L(); k++) {
				if ((k + Orbitals[a].L() + Orbitals[b].L()) % 2 != 0) continue;
				angular = Constant::Wigner3j(Orbitals[a].L(), k, Orbitals[b].L(), 0, 0, 0);
				angular = 0.5*angular*angular;
				eAB -= angular*U.R_k(k, Orbitals[a], Orbitals[b], Orbitals[b], Orbitals[a]);
			}
			Result -= Orbitals[a].occupancy()*wght*eAB;
		}
	}

	return Result;
}

double HartreeFock::Conf_En(std::vector<RadialWF> &Orbitals, std::vector<RadialWF> &Virtual, Potential &U)
{
	// Calculate total energy of electronic configuration.
	// Average over configurations is assumed.
	// It's WRONG, need to do through the R_k(abab), not the potential!!
	double Result = 0;
	// Test for hydrogenic syste - only one electron.
	int OneElec = 0;

	std::vector<RadialWF*> Occupied(0);
	for (int i = 0; i < Orbitals.size(); i++) {
		if (Orbitals[i].occupancy() != 0) {
			Occupied.push_back(&Orbitals[i]);
			OneElec += Orbitals[i].occupancy();
		}
	}
	for (int i = 0; i < Virtual.size(); i++) {
		if (Virtual[i].occupancy() != 0) {
			Occupied.push_back(&Virtual[i]);
			OneElec += Orbitals[i].occupancy();
		}
	}
	if (OneElec == 1) return Occupied[0]->Energy;

	for (auto& Orb: Orbitals) Result += Orb.occupancy()*Orb.Energy;
	// Get orbital Coulomb interaction energies.
	double wght = 0; // Complete shell occupancy.
	double eAB = 0;
	double angular = 0;
	for (int a = 0; a < Occupied.size(); a++) {
		for (int b = a; b < Occupied.size(); b++) {
			wght = Occupied[a]->occupancy();
			if (b != a) wght *= Occupied[b]->occupancy();
			else wght *= 0.5*(4*Occupied[b]->L() + 2)*(Occupied[b]->occupancy() - 1)/(4*Occupied[b]->L() + 1);
			if (wght == 0) continue;

			eAB = U.R_k(0, *Occupied[a], *Occupied[b], *Occupied[a], *Occupied[b]);
			for (int k = 0; k <= Occupied[a]->L() + Occupied[b]->L(); k++) {
				if ((k + Occupied[a]->L() + Occupied[b]->L())%2 != 0) continue;
				angular = Constant::Wigner3j(Occupied[a]->L(), k, Occupied[b]->L(), 0, 0, 0);
				angular = 0.5*angular*angular;
				eAB -= angular*U.R_k(k, *Occupied[a], *Occupied[b], *Occupied[b], *Occupied[a]);
			}

			Result -= wght*eAB;
		}
	}

	return Result;
}

void HartreeFock::MixOldNew(RadialWF * New_Orbital, RadialWF * Old_Orbital)
{
	// Stabilize convergence of the HF routine. Mix both orbitals.
	std::vector<double> density(lattice->size(), 0);
	Adams I(*lattice, 5);
	double Norm = 0;
	int Infty = max(New_Orbital->pract_infinity(), Old_Orbital->pract_infinity());

	for (int i = 0; i < lattice->size(); i++) {
		New_Orbital->F[i] += Old_Orbital->F[i];
		density[i] = New_Orbital->F[i] * New_Orbital->F[i];
	}

	Norm = I.Integrate(&density, 0, Infty);
	New_Orbital->scale(1./sqrt(Norm));
}
