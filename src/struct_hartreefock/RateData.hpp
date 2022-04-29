#pragma once

#include <iostream>
#include <filesystem>
#include <vector>
// The containers storing rates as used by rate equation solver


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
		// flags for specifying what to save
		// static const char PHOTO = 0x01;
		// static const char FLUOR = 0x02;
		// static const char AUGER = 0x04;
		// static const char EII   = 0x08;
		std::vector<std::string> index_names = std::vector<std::string>(0);
		std::string name = "";
		double nAtoms = 1.;// atomic number density
		// double R = 189.; // 100nm focal spot radius.
		unsigned int num_conf = 1;

		// the true rate data
		std::vector<RateData::Rate> Photo = std::vector<RateData::Rate>(0);
		std::vector<RateData::Rate> Fluor = std::vector<RateData::Rate>(0);
		std::vector<RateData::Rate> Auger = std::vector<RateData::Rate>(0);
		std::vector<RateData::EIIdata> EIIparams = std::vector<RateData::EIIdata>(0);
	};

	// Takes the dtailed balance dual
	std::vector<InverseEIIdata> inverse(const std::vector<EIIdata>& eiiVec);

	// Functional-style save functions modeled on printf/scanf
	void save_csv(const std::filesystem::path& root, const std::vector<Rate>& rates);
	void save_csv(const std::filesystem::path& root, const std::vector<EIIdata>& rates);
	void save_csv(const std::filesystem::path& root, const Atom& a);

	void load_csv(const std::filesystem::path& root, std::vector<Rate>& rates);
	void load_csv(const std::filesystem::path& root, std::vector<EIIdata>& rates);
	void load_csv(const std::filesystem::path& root, Atom& a);
}