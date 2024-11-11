/**
 * @file Constant.h
 * @brief This header contains some Physical constants,
 * @details note
functions for calculations of Wigner 3j, 6j symbols, and Clebsh-Gordan coefficients,
and some data containers used throughout the code. 
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

#include <vector>
#include <string>
#include <cassert>
using namespace std;

namespace Constant
{
	const double Pi = 3.1415926535897932384626433832795;
	const double Alpha = 0.0072973525698;
	const double Alpha2 = Constant::Alpha * Constant::Alpha;
	const double Fm = 1.8897261246 / 100000;

	//**// Conversions
	// Read these as "One au is 0.024 fs".
	// Note unit conversion direction is consistent, i.e. always divide to convert to atomic units.
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
	//**//
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
		vector<double> val;
	};

	struct polarize
	{
		vector<int> reference;// Contains only occupancies of Orbitals, no Virtual included.
		double refEnergy;// Energy of the reference configuration.
		/* E1-selected excited configurations.
		excited[][0] - orbital in reference from which to excite.
		excited[][1] - orbital in reference to which to excite.
		*/
		vector<vector<int>> excited;// Includes both Orbitals and Virtual.
		vector<double> extEnergy;
		vector<double> Dipoles; // Reduced transition dipole matrix elements form 'reference' to 'excited'.
		int index = 0;
	};

	struct bound_transport //
	{
		int from_heavy; // Index of configuration of heavy atom.
		int from_light; // Index of configuration of light atom.
		int to_heavy; // Index of configuration for heavy atom corresponding to allowed configuration with lowest energy for the complex
		int to_light; // Index of configuration for light atom corresponding to allowed configuration with lowest energy for the complex
	};	
	// Idea: If donate electron, go to donator_index. If receive electron, go to receiver index. 
	struct energy_config
	{
		int index;
		double valence_energy; // Energy of valence electron - for our bound transport correction the electron automatically chooses to hop in to a valence orbital if its binding energy is stronger.
		int receiver_index;
		int donator_index;
	};	
}

typedef std::vector<double> bound_t; // TODO I'm debating removing this since there are lots of std::vector<double> declarations that this makes confusing -S.P.


namespace RateData {

	struct EIIdata
	{
		//EIIdata() : init(0), fin(vector<int>(0)), occ(vector<int>(0)),ionB(vector<float>(0)),kin(vector<float>(0)) {}

		vector<int> fin; // final states
		int init; // initial state
		vector<int> occ; // occupancy of state
		vector<float> ionB; // ion binding energy
		vector<float> kin; // u for atom in this state (see Kim and Rudd BEB for details)

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
	vector<InverseEIIdata> inverse(const vector<EIIdata>& eiiVec);

	struct Rate
	{
		double val = 0;
		long int from = 0;
		long int to = 0;
		double energy = 0;
	};

	struct Atom
	{
		bool bound_free_excluded = false; // Whether to skip calculation of EII and TBR for this species.
		vector<string> index_names = vector<string>(0);
		std::string name = "";
		double nAtoms = 1.;// atomic number density
		// double R = 189.; // 100nm focal spot radius.
		unsigned int num_conf = 1;
		vector<RateData::Rate> Photo = vector<RateData::Rate>(0);
		vector<RateData::Rate> Fluor = vector<RateData::Rate>(0);
		vector<RateData::Rate> Auger = vector<RateData::Rate>(0);
		vector<RateData::EIIdata> EIIparams = vector<RateData::EIIdata>(0);
		// Tacked on energy_config here.
		vector<CustomDataType::energy_config> EnergyConfig = vector<CustomDataType::energy_config>(0);
	};

	bool ReadRates(const string & input, vector<RateData::Rate> & PutHere);
	bool ReadEIIParams(const string & input, vector<RateData::EIIdata> & PutHere);
	bool ReadDecayRates(const string & rate_location, const string & rate_file_type, vector<RateData::Rate> & PutHere,int num_allowed_configs);
	bool InterpolateRates(const string & rate_location, const string & rate_file_type, vector<RateData::Rate> & PutHere, double photon_energy,double allowed_interp = 1000);
	void WriteRates(const string& fname, const vector<RateData::Rate>& rateVector);
	void WriteEIIParams(const string& fname, const vector<RateData::EIIdata>& eiiVector);
}


static const double Moulton_5[5] = { 251. / 720., 646. / 720., -264. / 720., 106. / 720., -19. / 720. }; //Adams-Moulton method
static const double Bashforth_5[5] = { 1901. / 720., -1378. / 360., 109. / 30., -637. / 360., 251. / 720. }; //Adams-Bashforth method

// Guassian Quadrature, you may also wish to check SplineIntegral.h
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


// Sourced from https://pomax.github.io/bezierinfo/legendre-gauss.html - S.P. Better results than order 10 of same dec. places, usage is just slower. 
static const double gaussW_64[64] = {0.048690957,0.048690957,0.048575467,0.048575467,0.048344762,0.048344762,0.047999389,0.047999389,0.047540166,
0.047540166,0.046968183,0.046968183,0.046284797,0.046284797,0.045491628,0.045491628,0.044590558,0.044590558,0.043583725,0.043583725,0.042473515,
0.042473515,0.041262563,0.041262563,0.039953741,0.039953741,0.038550153,0.038550153,0.037055129,0.037055129,0.035472213,0.035472213,0.033805162,
0.033805162,0.032057928,0.032057928,0.030234657,0.030234657,0.028339673,0.028339673,0.02637747,0.02637747,0.024352703,0.024352703,0.022270174,0.022270174,
0.020134823,0.020134823,0.017951716,0.017951716,0.01572603,0.01572603,0.013463048,0.013463048,0.011168139,0.011168139,0.00884676,0.00884676,0.006504458,
0.006504458,0.004147033,0.004147033,0.001783281,0.001783281};

static const double gaussX_64[64] = {-0.024350293,0.024350293,-0.072993122,0.072993122,-0.121462819,0.121462819,-0.16964442,0.16964442,-0.217423644,0.217423644,
-0.264687162,0.264687162,-0.311322872,0.311322872,-0.357220158,0.357220158,-0.402270158,0.402270158,-0.446366017,0.446366017,-0.489403146,0.489403146,-0.531279464,
0.531279464,-0.571895646,0.571895646,-0.611155355,0.611155355,-0.648965471,0.648965471,-0.685236313,0.685236313,-0.71988185,0.71988185,-0.752819907,0.752819907,
-0.783972359,0.783972359,-0.813265315,0.813265315,-0.840629296,0.840629296,-0.865999398,0.865999398,-0.889315446,0.889315446,-0.910522137,0.910522137,-0.929569172,
0.929569172,-0.946411375,0.946411375,-0.9610088,0.9610088,-0.973326828,0.973326828,-0.983336254,0.983336254,-0.991013371,0.991013371,-0.996340117,0.996340117,
-0.999305042,0.999305042};

/// Ensures we have consistent removal of decimal places which is important for loading sims.
/// Necessary due to high precision of Constant::fs_per_au
static const double loading_t_precision = 9;