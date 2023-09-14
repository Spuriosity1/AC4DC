#ifndef RATEDATA_INCLUDE_H
#define RATEDATA_INCLUDE_H

#include <string>
#include <vector>
#include <nlohmann/json.hpp>

using json= nlohmann::json;

namespace RateData
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

	//form factor for Q_mesh values in FormFactor class
	struct ffactor
	{	
		inline static ffactor make_ffactor(int idx, const std::vector<double>& val){
			ffactor r;
			r.index = idx;
			r.val = val;
			return r;
		}
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

	struct bound_transport //
	{
		int from_heavy; // Index of configuration of heavy atom.
		int from_light; // Index of configuration of light atom.
		int to_heavy; // Index of configuration for heavy atom corresponding to allowed configuration with lowest energy for the complex
		int to_light; // Index of configuration for light atom corresponding to allowed configuration with lowest energy for the complex
	};


	/**
	 * @brief Represents a valence electron process for bound transport
	 * 
	 *  Idea: If donate electron, go to donator_index. If receive electron, go to receiver index. 
	 */
	struct energy_config
	{
		int index;
		double valence_energy; // Energy of valence electron - for our bound transport correction the electron automatically chooses to hop in to a valence orbital if its binding energy is stronger.
		int receiver_index;
		int donator_index;
	};	

	struct EIIdata
	{
		//EIIdata() : init(0), fin(std::vector<int>(0)), occ(std::vector<int>(0)),ionB(std::vector<float>(0)),kin(std::vector<float>(0)) {}

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
		double atomic_density = 1.;// atomic number density
		// double R = 189.; // 100nm focal spot radius.
		unsigned int num_conf = 1;
		std::vector<Rate> Photo = std::vector<Rate>(0);
		std::vector<Rate> Fluor = std::vector<Rate>(0);
		std::vector<Rate> Auger = std::vector<Rate>(0);
		std::vector<EIIdata> EIIparams = std::vector<EIIdata>(0);
		// Tacked on energy_config here.
		std::vector<energy_config> EnergyConfig = std::vector<energy_config>(0);
	};

	void from_json(const json& j, energy_config& ec);
	void from_json(const json& j, Rate& rate);
	void from_json(const json& j, EIIdata& rate);
	void from_json(const json& j, Atom& atom);

	void to_json(json& dest, const energy_config& ec);
	void to_json(json& dest, const Rate& r);
	void to_json(json& dest, const EIIdata& eii);
	void to_json(json& dest, const Atom& a);

	
	// TODO
	// void save_atom_hdf5(const std::string& ofile, const Atom& );
	// void load_atom_hdf5(const std::string& ofile, const Atom& );

	// deprecated (maybe just altogether broken?)
	bool ReadRates(const std::string & input, std::vector<Rate> & PutHere);
	bool ReadEIIParams(const std::string & input, std::vector<EIIdata> & PutHere);
	void WriteRates(const std::string& fname, const std::vector<Rate>& ratevector);
	void WriteEIIParams(const std::string& fname, const std::vector<EIIdata>& eiiVector);
}

#endif