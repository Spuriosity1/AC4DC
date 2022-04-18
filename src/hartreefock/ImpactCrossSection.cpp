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

#include "ImpactCrossSection.hpp"
#include <cmath>
#include "Constant.hpp"
#include <cassert>

double ImpactCrossSection::sigmaBEB(double T, double B, double u, int occ)
{
	if (T < B) return 0;
	// Electron impact ionization cross-section
	// in Binary Encounter Bethe approximation : Kim, Rudd, PRA 50(5), 3954 (1994)
	const double t = T/B;
	if (t < 1) return 0;
	const double lnt = log(t);
	const double S = Constant::Pi*occ/B/B;
	constexpr double Q = 1;

	// double Result = S/(t + u + 1)*(0.5*(1 - 1./t/t)*lnt + (1 - 1/t - lnt/(1+t)));        // old
	double Result = (Q-2)*(1./t -1. + lnt /(1.+t)) + 0.5*Q*(1 - 1./t/t)*lnt; 			// new
	Result *= S/(t + u + 1); // new
	// double Result =S*(-(((-2 + Q)*(-1 + t))/t) + ((-2 + Q)*log(t))/(1 + t) + (Q*(-1 + t*t)*Log(t))/(2.*t*t))/(1 + t + u); // Mathematica garbage

	// #warning Using fake sigma! Edit ImpactCrossSection.cpp to fix
	// Result = u * (T-B)*(T-B)*(T-B)/24;
	return Result;
}


double ImpactCrossSection::sigmaBEBw1(double T, double B, double u, int occ)
{
	// Calculates the integral \int_0^(p*p/2 - B) dW W d(sigma)/dW using gaussian quadrature.
	double Result = 0;
	double Wmax = 0.5*(T - B);
	// Critically important to capture the interval where d(sigma)/dW substantially != 0.
	if (Wmax <= 0) return 0;
	double W = 0, k = 0.5*Wmax;
	// -1 .. 1 -> 0 .. Wmax [ W = 0.5*Wmax*(x + 1) ]
	for (int i = 0; i < 13; i++) {
		W = k*(gaussX_13[i] + 1);
		Result += gaussW_13[i]*W*ImpactCrossSection::DsigmaBEB(T, W, B, u, occ);
	}

	Result *= 0.5*Wmax;
	return Result;
}


// Represents dsigma(W | T)/dW, with normalistation
// sigma(T) = \int_0^{(T-B)/2} dsigma(W | T)/dW dW
double ImpactCrossSection::DsigmaBEB(double T, double W, double B, double u, int occ)
{
	// T - impactor electron energy.
	// W - ejected electron energy.
	// B - binding energy.
	// u - dimensionless kinetic energy, stored as 'kin' in EIIparams.
	// occ - occupancy of target orbital (overall multiplicative constant)
	// i - orbital index from which W is ejected.
	// Q[i] are set to 1. See commented functions if need to be calculated.
	// see eq. (52) of Kim & Rudd 1994 

	if ( T < B ) return 0; // Enforce the fact that impactor energy is greater than binding energy
	if ( W <0 || W > T-B ) return 0; // Enforce energy conservation

	const double t = T/B;
	const double w = W/B;

	const double S = Constant::Pi*occ/(B*B); // 4 pi a0^2 N R^2/B^2, atomic units (R = Rydeberg energy = 0.5 Ha)
	constexpr double Q = 1; // Approximation that avoids need to calculate M_i from paper. This can be restored if needed.

	const double x = 1./(w + 1);
	const double y = 1./(t - w);
	double Result = -(2-Q)*(x + y)/(t+1) + (2-Q)*(x*x + y*y) + Q*log(t)*(x*x*x + y*y*y);
	Result *= S/(t + u + 1);
	Result /= B; // THIS LINE IS KEY
	// Converts between dsigma/dw and dsigma/dW
	// the latter is needed

	// #warning Fake Dsigma in use!
	// Result = u*(W*(W - T+B) + (T-B)*(T-B)/4);
	return Result;
}


/*
void Plasma::Get_Ni(Grid & Lattice, std::vector<RadialWF> & Orbitals, std::vector<RadialWF> &Virtuals)
{
	// NOTE: for accurate calculation of N[i] parameter for BED model
	//       for orbital "i" there has to be virtual orbitals of at least
	//       two shells above orbitals[i].N().


	// Treat every orbital as if other don't exist.
	// Unless there is a hole, than account for that oscillator.
	N.clear();
	N.resize(Orbitals.size(), 0.);

	// First thing - check if there are holes in occupied orbitals.
	// If so, add the weighted (fractional occupancy) oscillator stregths
	// to N_i.
	// Use a tuple with (double, int, int) structure.
	std::vector<fluor> VirtOscStrg;
	fluor Osc_tmp;
	double ME = 0;
	int infinity = 0;
	std::vector<double> density(Lattice.size(), 0.);

	Adams AI(Lattice, 5);

	for (int i = 0; i < Orbitals.size(); i++) {
		N[i] = Orbitals[i].occupancy();
		for (int j = i + 1; j < Orbitals.size(); j++) {
			if (abs(Orbitals[j].L() - Orbitals[i].L()) != 1) continue;
			if (Orbitals[i].occupancy() == 4*Orbitals[i].L() + 2 &&
				Orbitals[j].occupancy() == 4*Orbitals[j].L() + 2) continue;

			if (Orbitals[i].occupancy() <= 0 &&
				Orbitals[j].occupancy() <= 0) continue;

			infinity = min(Orbitals[i].pract_infinity(), Orbitals[j].pract_infinity());
			for (int s = 0; s <= infinity; s++) {
				density[s] = Orbitals[i].F[s] * Lattice.R(s) *Orbitals[j].F[s];
			}
			ME = AI.Integrate(&density, 0, infinity);
			ME = 2*(Orbitals[j].Energy - Orbitals[i].Energy)*max(Orbitals[j].L(), Orbitals[i].L())*ME*ME
			/3;
			N[i] -= Orbitals[i].occupancy() * (1. - 1.*Orbitals[j].occupancy()/(4*Orbitals[j].L() + 2)) *
					  ME/(2*Orbitals[i].L() + 1);
			N[j] += Orbitals[j].occupancy() * (1. - 1.*Orbitals[i].occupancy()/(4*Orbitals[i].L() + 2)) *
					  ME/(2*Orbitals[j].L() + 1);
		}
	}

	// At this point we have N_i containing all bound-continuum oscillator strengths.
	// as well as the bound-virtual bound oscillators. The latter is calculated and
	// approximated by asymptotic series below.
	for (int i = 0; i < Orbitals.size(); i++) {
		if (Orbitals[i].occupancy() == 0) continue;
		for (int v = 0; v < Virtuals.size(); v++) {
			if (abs(Virtuals[v].L() - Orbitals[i].L()) != 1) continue;
			infinity = Orbitals[i].pract_infinity();
			for (int j = 0; j <= infinity; j++) {
				density[j] = Orbitals[i].F[j] * Lattice.R(j) * Virtuals[v].F[j];
			}
			ME = AI.Integrate(&density, 0, infinity);

			Osc_tmp.val = 2*(Virtuals[v].Energy - Orbitals[i].Energy)*max(Virtuals[v].L(), Orbitals[i].L())*ME*ME/3/(2*Orbitals[i].L() + 1);
			Osc_tmp.fill = v;
			Osc_tmp.hole = i;
			VirtOscStrg.push_back(Osc_tmp);
		}
	}

	// Now we have a set of bound-virtual bound Oscillator Strengths.
	// Work out N_i using sum rules and hydrogenic-style asymptotic for
	// bound-bound Oscillator Strength sum beyond the calculated values.
	std::vector<fluor*> OscStrg_i;
	int last_virt_ind = 0, last_Osc_ind = 0;
	int orb_ind = 0, virt_ind = 0;
	double BoundOscStrg = 0;
	double AsmpRowSumm = Constant::RiemannZeta3;
	for (int i = 0; i < Orbitals.size(); i++) {
		for (int l = Orbitals[i].L() - 1; l <= Orbitals[i].L() + 1; l += 2) {
			OscStrg_i.clear();
			BoundOscStrg = 0;
			AsmpRowSumm = Constant::RiemannZeta3;
			for (auto& elem: VirtOscStrg) {
				orb_ind = elem.hole;
				virt_ind = elem.fill;
				if (orb_ind == i && Virtuals[virt_ind].L() == l)	OscStrg_i.push_back(&elem);
			}
			if (OscStrg_i.size() == 0) continue;
			last_virt_ind = -1;//OscStrg_i[0]->fill;
			last_Osc_ind = -1;
			for (int lst = 0; lst < OscStrg_i.size(); lst++) {
				if (OscStrg_i[lst]->fill > last_virt_ind) {
					last_virt_ind = OscStrg_i[lst]->fill;
					last_Osc_ind = lst;
				}
			}
			if (last_Osc_ind == -1) continue;
			// Start with asymptotic bound-bound contribution.
			BoundOscStrg = OscStrg_i[last_Osc_ind]->val * pow(Virtuals[last_virt_ind].N(), 3);
			// Multiply pre-factor by RiemannZeta3 and remove double counted bound-bound transitions.
			for (int corrN = 1; corrN < Virtuals[last_virt_ind].N(); corrN++) AsmpRowSumm -= 1./pow(corrN, 3);
			BoundOscStrg *= AsmpRowSumm;
			// Add the remaining bound-bound Osc. Strg.
			for (int os = 0; os < OscStrg_i.size(); os++) {
				if (os == last_Osc_ind) continue;
				BoundOscStrg += OscStrg_i[os]->val;
			}
			BoundOscStrg *= Orbitals[i].occupancy();
			N[i] -= BoundOscStrg;
		}
	}
}


void Plasma::Get_Qi(Grid & Lattice, std::vector<RadialWF> & Orbitals, std::vector<RadialWF> &Virtuals)
{
	// NOTE: for accurate calculation of Q[i] parameter for BEB model
	//       for orbital "i" there has to be virtual orbitals of at least
	//       two shells above orbitals[i].Q().
	//
	//		 f_{i-n}/(E_n - E_i) = 2/3/(2*l_i+1)|<i||r||n>|^2  for averaged
	//                                                         over projection OS
	//
	// Q_i = 2*\int_0^{infty} dw df_i/dw 1/(1 + w) / orbitals[i].occupancy();
	// Q_i = |E_i|(T1 + T2 + T3)
	// T1 = -2/3/(2*l_i+1) sum_n l_> |<i||r||n>|^2 - "n" includes virtuals except last, holes in core.
	// T2 = -2/3/(2*l_i+1) sum_n l_> |<i||r||n>|^2 - "n" from last virtual to infinity.
	//                                               Asymptotic for |<i||r||n>|^2 is used.
	// T3 = 2/3/(2*l_i+1)<i||r^2||i>
	//
	// l_> = max(l_i, l_n)


	// Treat every orbital as if other don't exist.
	// Unless there is a hole, than account for that oscillator.
	Q.clear();
	Q.resize(Orbitals.size(), 0.);
	std::vector<double> T1(Orbitals.size(), 0.);
	std::vector<double> T2(Orbitals.size(), 0.);
	std::vector<double> T3(Orbitals.size(), 0.);
	std::vector<fluor> CoreCoreDipl;
	std::vector<fluor> CoreVirtDipl;

	fluor Tmp;
	double ME = 0;
	int infinity = 0;
	std::vector<double> density(Lattice.size(), 0.);
	Adams AI(Lattice, 5);
	// T3 terms.
	for (int i = 0; i < Orbitals.size(); i++) {
		infinity = Orbitals[i].pract_infinity();
		for (int s = 0; s <= infinity; s++) {
			density[s] = Orbitals[i].F[s] * Lattice.R(s) * Lattice.R(s) * Orbitals[i].F[s];
		}
		T3[i] = (2*Orbitals[i].L()+1)*AI.Integrate(&density, 0, infinity);
	}

	for (int i = 0; i < Orbitals.size(); i++) {
		for (int j = i + 1; j < Orbitals.size(); j++) {
			if (Orbitals[i].occupancy() == 0 && Orbitals[j].occupancy() == 0) continue;
			if (abs(Orbitals[i].L() - Orbitals[j].L()) != 1) continue;
			infinity = min(Orbitals[i].pract_infinity(), Orbitals[j].pract_infinity());
			for (int s = 0; s <= infinity; s++) {
				density[s] = Orbitals[i].F[s] * Lattice.R(s) * Orbitals[j].F[s];
			}
			ME = AI.Integrate(&density, 0, infinity);
			Tmp.val = ME*ME*max(Orbitals[i].L(), Orbitals[j].L());
			Tmp.fill = j;
			Tmp.hole = i;
			CoreCoreDipl.push_back(Tmp);
		}
	}

	// Assemble core-core contributions to T1.
	for (auto& cco: CoreCoreDipl) {
		int i = cco.hole;
		int j = cco.fill;
		T1[i] -= cco.val * (1. - 1.*Orbitals[j].occupancy()/(4*Orbitals[j].L() + 2));
		T1[j] -= cco.val * (1. - 1.*Orbitals[i].occupancy()/(4*Orbitals[i].L() + 2));
	}

	// Calculate oscillator strength for core-virtual bound states and get the last virtual state for each orbital.
	std::vector<int> last_virt;// [j1, j2, ...] - maximum "N" quantum number of give L = virtuals[j].L().
	int N_max = 0;
	bool check = false;
	for (int i = 0; i < Virtuals.size(); i++) {
		N_max = i;
		for (int j = i+1; j < Virtuals.size(); j++) {
			if (Virtuals[i].L() == Virtuals[j].L() && Virtuals[i].N() < Virtuals[j].N()) N_max = j;
		}
		check = false;
		for (auto& lv: last_virt) {
			if (Virtuals[lv].L() == Virtuals[N_max].L()) check = true;
		}
		if (!check) last_virt.push_back(N_max);
	}

	for (int i = 0; i < Orbitals.size(); i++) {
		if (Orbitals[i].occupancy() == 0) continue;
		for (int j = 0; j < Virtuals.size(); j++) {
			if (abs(Orbitals[i].L() - Virtuals[j].L()) != 1) continue;
			infinity = min(Orbitals[i].pract_infinity(), Virtuals[j].pract_infinity());
			for (int s = 0; s <= infinity; s++) {
				density[s] = Orbitals[i].F[s] * Lattice.R(s) * Virtuals[j].F[s];
			}
			ME = AI.Integrate(&density, 0, infinity);
			Tmp.val = max(Virtuals[j].L(), Orbitals[i].L())*ME*ME;
			Tmp.fill = j;
			Tmp.hole = i;
			CoreVirtDipl.push_back(Tmp);
		}
	}
	// Finish assembling T1;
	for (auto& cvo: CoreVirtDipl) {
		int i = cvo.hole;
		int j = cvo.fill;
		check = false;
		for (auto& lv: last_virt) {
			if (lv == j) check = true;
		}
		if (!check) T1[i] -= cvo.val;
	}

	// Work out T2 using hydrogenic-style asymptotic for
	// core-virtual dipole ME. Sum to infinity, see Bethe book for details.
	int last_virt_ind = 0, last_Osc_ind = 0;
	int orb_ind = 0, virt_ind = 0;
	double SumOsc = 0;
	double AsmpRowSum = Constant::RiemannZeta3;
	for (int i = 0; i < Orbitals.size(); i++) {
		SumOsc = 0;
		for (int l = Orbitals[i].L() - 1; l <= Orbitals[i].L() + 1; l += 2) {
			if (l < 0) continue;

			last_virt_ind = -1;
			for (auto& vl: last_virt) {
				if (Virtuals[vl].L() == l) {
					last_virt_ind = vl;
					break;
				}
			}
			if (last_virt_ind == -1) continue;

			last_Osc_ind = -1;
			for (int b = 0; b < CoreVirtDipl.size(); b++) {
				if (CoreVirtDipl[b].fill == last_virt_ind && CoreVirtDipl[b].hole == i) {
					last_Osc_ind = b;
					break;
				}
			}
			if (last_Osc_ind == -1) continue;

			AsmpRowSum = Constant::RiemannZeta3;

			SumOsc = CoreVirtDipl[last_Osc_ind].val * pow(Virtuals[last_virt_ind].N(), 3);
			for (int corrN = 1; corrN < Virtuals[last_virt_ind].N(); corrN++) {
				AsmpRowSum -= 1./pow(corrN, 3);
			}

			T2[i] -= SumOsc * AsmpRowSum;
		}
	}
	for (int i = 0; i < Orbitals.size(); i++) {
		Q[i] = 2*(T1[i] + T2[i] + T3[i])*(-1)*Orbitals[i].Energy*2/3/(2*Orbitals[i].L() + 1);
	}
	// Comment in Q[i] orbitals[i].occupancy() cancells out according to eq. (28) in the paper.
}
*/

/*
void Plasma::setup_EII(std::vector<RadialWF> &Virtual, double k_min, double k_max)
{

	// Calculate oscillatory strengths dF_dW in a given range of p_e.
	//dF_dW.clear();
	//dF_dW.resize(orbitals.size());
	//for (auto& df: dF_dW) df.resize(p_e.size());

	if (k_min < 0.5) k_min = 0.5; // For low values of p_e calculate p_e = 1. a.u. (w = 0.5 a.u. = 13.6 eV);

	double Infinity = 0;
	for (auto& Orb: orbitals) {
		if(lattice.R(Orb.pract_infinity()) > Infinity && Orb.occupancy() != 0) Infinity = lattice.R(Orb.pract_infinity());
	}

	if (Infinity < 5*6.28/k_min) Infinity = 5*6.28/k_min;
	// Interpolation lattice.
	CntLattice = Grid(Infinity, 0.03/k_max, u.NuclCharge());

	Interpolation W(6);

	CntOrbitals.clear();
	CntOrbitals.resize(orbitals.size());

	for (int i = 0; i < orbitals.size(); i++) {
		int j = 0;
		if (input.Exited_Pot_Model() != "V_N-1" && orbitals[i].occupancy() == 0) continue;
		CntOrbitals[i].Energy = orbitals[i].Energy;
		CntOrbitals[i].set_N(orbitals[i].N());
		CntOrbitals[i].set_L(orbitals[i].L());
		CntOrbitals[i].set_occupancy(orbitals[i].occupancy());
		if (orbitals[i].F[0] != 0) W.RecalcWF(orbitals[i], lattice, CntOrbitals[i], CntLattice);
		Infinity = lattice.R(orbitals[i].pract_infinity());
		j = 0;
		while (CntLattice.R(j) < Infinity && j < CntLattice.size() - 1) j++;
		CntOrbitals[i].set_infinity(j);
	}

	CntU = Potential(&CntLattice, u.NuclCharge());
	CntU.LDA_upd_dir(CntOrbitals);

	RadialWF Cont(CntLattice.size());
	RadialWF dCont(CntLattice.size());
	Cont.set_N(-1);
	Cont.set_occupancy(1);
	Cont.set_infinity(CntLattice.size() - 1);
	Cont.set_L(0);

	// Update here if Input.Hamiltonian == "LDA".
	if (input.Hamiltonian() == 0) U.HF_upd_dir(&Cont, CntOrbitals);
	else U.LDA_upd_dir(CntOrbitals);

	double ME = 0, dME = 0;
	Adams I1(CntLattice, 5);

	for (auto& l: L) {
		for (int j = 0; j < p_e.size(); j++) {
			Cont.set_L(l);
			Cont.Energy = 0.5*p_e[j]*p_e[j];
			dCont = Cont;
			dCont.Energy *= 1.01;

			if (IntegrateContinuum(CntLattice, U, CntOrbitals, &Cont) < 0
				|| IntegrateContinuum(CntLattice, U, CntOrbitals, &dCont) < 0) {
				log << "====================================================================" << endl;
				log << "setupe_EII: Continuum didn't converge: " << endl;
				for (int i = 0; i < orbitals.size(); i++)
				{
					log << i + 1 << ") n = " << orbitals[i].N() << " l = " << orbitals[i].L()
						<< " Energy = " << orbitals[i].Energy << " Occup = " << orbitals[i].occupancy() << endl;
				}
				log.flush();
			}

			for (int i = 0; i < orbitals.size(); i++) {
				if (orbitals[i].occupancy() <= 0) continue;
				if (abs(l - orbitals[i].L()) != 1) continue;

				infinity = CntOrbitals[i].pract_infinity();
				density.clear();
				density.resize(infinity + 1);

				if (input.Gauge() == "length") {
					for (int s = 0; s < density.size(); s++)	{
						density[s] = CntLattice.R(s) * Cont.F[s] * CntOrbitals[i].F[s];
					}
					ME = I1.Integrate(&density, 0, infinity);
					for (int s = 0; s < density.size(); s++)	{
						density[s] = CntLattice.R(s) * dCont.F[s] * CntOrbitals[i].F[s];
					}
					dME = I1.Integrate(&density, 0, infinity);
				}
				else { // Velocity gauge.
					double ang_coeff =  0.5*(CntOrbitals[i].L() - Cont.L())*(CntOrbitals[i].L() + Cont.L() + 1);
					for (int s = 0; s < density.size(); s++)	{
						density[s] = Cont.F[s] *(CntOrbitals[i].G[s] + ang_coeff * CntOrbitals[i].F[s]/CntLattice.R(s)) ;
					}
					ME = I1.Integrate(&density, 0, infinity)/(Cont.Energy - CntOrbitals[i].Energy);
					for (int s = 0; s < density.size(); s++)	{
						density[s] = dCont.F[s] *(CntOrbitals[i].G[s] + ang_coeff *CntOrbitals[i].F[s]/CntLattice.R(s)) ;
					}
					dME = I1.Integrate(&density, 0, infinity)/(dCont.Energy -CntOrbitals[i].Energy);
				}

				double A = 2*(Cont.Energy -CntOrbitals[i].Energy)/3/(2*CntOrbitals[i].L() + 1) * ME * ME;
				A -= 2*(dCont.Energy -CntOrbitals[i].Energy)/3/(2*CntOrbitals[i].L() + 1) * dME * dME;
				A /= (Cont.Energy - dCont.Energy);
				if (Cont.L() >CntOrbitals[i].L()) A *= Cont.L();
				else A *=CntOrbitals[i].L();

				dF_dW[i][j] += A;
			}
		}
	}

}
*/
