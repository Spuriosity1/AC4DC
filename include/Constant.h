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

/* This header contains some Physical constants,
functions for calculations of Wigner 3j, 6j symbols, and Clebsh-Gordan coefficients,
and some data containers used througout the code. */
#include <vector>>
#include <string>
#include <cassert>


namespace Constant
{
	const double Pi = 3.1415926535897932384626433832795;
	const double Alpha = 0.0072973525698;
	const double Alpha2 = Constant::Alpha * Constant::Alpha;
	const double Fm = 1.8897261246 / 100000;

	// Conversions
	// Read these as "One au is 0.024 fs"
	const double fs_per_au = 0.024188843265857;//femtosecond in au
	const double eV_per_Ha = 27.211385;//electron volts in atomic units
	const double J_per_eV = 1.60217662e-19;
	// const double Intensity_in_au = 6.434;// x10^15 W/cm^2 //3.50944758;//x10^2 W/nm^2
	// const double Jcm2_per_Haa02 = 1./6.4230434293;//J/cm^2 
	const double au2_in_barn = 5.2917721067*5.2917721067*1000000;//atomic units to Barns.
	const double au2_in_Mbarn = 5.2917721067*5.2917721067;//atomic units to Mega Barns.
	const double RiemannZeta3 = 1.202056903159594;
	const double Angs_per_au = 0.52917721067; // Bohr radius = 1 atomic unit in Angstrom.
	const double cm_per_au = Angs_per_au*1e-8;
	const double kb_eV = 8.617333262145e-5; // Boltzmann constant, electronvolt per Kelvin
	const double kb_Ha = 8.617333262145e-5/eV_per_Ha; // Boltzmann constant, Ha per Kelvin

	const double Jcm2_per_Haa02 = J_per_eV/cm_per_au/cm_per_au * eV_per_Ha;

	double Wigner3j(double, double, double, double, double, double);
}

namespace CustomDataType
{
	struct photo //photoionisation rate
	{
		double val;//value of the rate in a.u.
		int hole;//orbital with hole
		double energy; // Ha
	};

	struct fluor//fluorescence rate
	{
		double val;//value of the rate in a.u.
		int fill;//filled orbital
		int hole;//orbital with hole
	};

	struct auger//auger decay rate
	{
		double val;
		int hole;
		int fill;
		int eject;
		double energy = 0;
	};

	struct ffactor//form factor for Q_mesh values in FormFactor class
	{
		int index;
		std::vector<double> val;
	};

	struct polarize
	{
		std::vector<int> reference;// Contains only occupancies of Orbitals, no Virtual included.
		double refEnergy;// Energy of the reference configuration.
		/* E1-selected excited configurations.
		excited[][0] - orbital in reference from which to excite.
		excited[][1] - orbital in reference to which to excite.
		*/
		std::vector<std::vector<int>> excited;// Includes both Orbitals and Virtual.
		std::vector<double> extEnergy;
		std::vector<double> Dipoles; // Reduced transition dipole matrix elements form 'reference' to 'excited'.
		int index = 0;
	};
}

typedef std::vector<double> bound_t;


namespace RateData {

	struct EIIdata
	{
		int init; // initial state
		std::vector<int> fin; // final states
		std::vector<int> occ; // occupancy of state
		std::vector<float> ionB; // ion binding energy
		std::vector<float> kin; // u for atom in this state (see Kim and Rudd BEB for details)

		void resize(size_t n)
		{
			fin.resize(n);
			occ.resize(n);
			ionB.resize(n);
			kin.resize(n);
		}

		void push_back(int f, int o, float B, float U) {
			fin.push_back(f);
			occ.push_back(o);
			ionB.push_back(B);
			kin.push_back(U);
		}

		size_t size() {
			#ifdef DEBUG
			assert(fin.size() == occ.size());
			assert(fin.size() == ionB.size());
			assert(fin.size() == kin.size());
			#endif
			return fin.size();
		}
	};

	// Though this structure is identical (but for the names) to EIIdata, it is made deliberately incompatible
	// to prevent confusion.
	typedef EIIdata InverseEIIdata;

	// Reorganises a EIIData tree by final index rather than initial
    // Used for Q_TBR
	std::vector<InverseEIIdata> inverse(const std::vector<EIIdata>& eiiVec);

	struct Rate
	{
		double val = 0;
		long int from = 0;
		long int to = 0;
		double energy = 0;
	};

	struct Atom
	{
		std::vector<std::string> index_names = std::vector<std::string>(0);
		std::string name = "";
		double nAtoms = 1.;// atomic number density
		// double R = 189.; // 100nm focal spot radius.
		unsigned int num_conf = 1;
		std::vector<RateData::Rate> Photo = std::vector<RateData::Rate>(0);
		std::vector<RateData::Rate> Fluor = std::vector<RateData::Rate>(0);
		std::vector<RateData::Rate> Auger = std::vector<RateData::Rate>(0);
		std::vector<RateData::EIIdata> EIIparams = std::vector<RateData::EIIdata>(0);
	};

	bool ReadRates(const std::string & input, std::vector<RateData::Rate> & PutHere);
	bool ReadEIIParams(const std::string & input, std::vector<RateData::EIIdata> & PutHere);
	void WriteRates(const std::string& fname, const std::vector<RateData::Rate>& rateVector);
	void WriteEIIParams(const std::string& fname, const std::vector<RateData::EIIdata>& eiiVector);
}


static const double Moulton_5[5] = { 251. / 720., 646. / 720., -264. / 720., 106. / 720., -19. / 720. }; //Adams-Moulton method
static const double Bashforth_5[5] = { 1901. / 720., -1378. / 360., 109. / 30., -637. / 360., 251. / 720. }; //Adams-Bashforth method
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
