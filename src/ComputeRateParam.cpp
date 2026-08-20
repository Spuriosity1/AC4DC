/**
 * @file ComputeRateParam.cpp
 * @brief @copybrief ComputeRateParam.h
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
#include "ComputeRateParam.h"

inline bool exists_test(const std::string&);
vector<double> generate_dT(int);
vector<double> generate_T(vector<double>&);
vector<double> generate_I(vector<double>&, double, double);
void SmoothOrigin(vector<double>&, vector<double>&);


inline bool exists_test(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}


using namespace CustomDataType;

// Called when atomic input is relevant.
int ComputeRateParam::SolveFrozen(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	// Calculates cross sections needed to solve rate equations.
	// Final_occ defines the lowest possible occupancies for the initial orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	assert(Max_occ.size() == orbitals.size());
	//DecayRates::set_Max_occ(Max_occ,orbitals);  // used when max_occ isn't valid. 

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return 1;

	// Make an output folder, if it doesn't exist yet
	if (!exists_test("./output")) {
		mkdir("output", ACCESSPERMS);
	}


	string RateLocation = "./output/" + input.Name() + "/Xsections/";

	// Make a data folder inside output
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	// Make a rate-containing folder inside the data folder
	if (!exists_test(RateLocation)) {
		mkdir(RateLocation.c_str(), ACCESSPERMS);
	}

	cout << "Check if there are pre-calculated rates..." << endl;
	bool have_Pht = RateData::ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool have_Flr = RateData::ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool have_Aug = RateData::ReadRates(RateLocation + "Auger.txt", Store.Auger);
	bool saveFF = true; // TODO: Make MolInp control this

	if (have_Pht) printf("Photoionization rates found. Reading...\n");
	if (have_Flr) printf("Fluorescence rates found. Reading...\n");
	if (have_Aug) printf("Auger rates found. Reading...\n");

	string PolarFileName = "./output/Polar_" + input.Name() + ".txt";

	if ( true || !have_Pht || !have_Flr || !have_Aug )
	{
		cout <<"======================================================="<<endl;
		cout << "Total number of configurations: " << dimension << endl;
		cout <<" Beginning Hartree-Fock Frozen calculations... "<<endl;
		cout <<"======================================================="<<endl;
		RateData::Rate Tmp;
		vector<RateData::Rate> LocalPhoto(0);
		vector<RateData::Rate> LocalFluor(0);
		vector<RateData::Rate> LocalAuger(0);
		vector<ffactor> LocalFF(0);

		#pragma omp parallel default(none) num_threads(input.Num_Threads()) \
		shared(cout, runlog, have_Aug, have_Flr, have_Pht, saveFF) private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor, LocalFF)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0;i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				cout << "[HF Frozen] configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (size_t j = 0;j < Orbitals.size(); j++)
				{
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				Potential U(&Lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(Lattice, Orbitals, U, input, runlog);

				DecayRates Transit(Lattice, Orbitals, u, input);

				// ======= Experimental =========
				if (saveFF) LocalFF.push_back({i, Transit.FT_density()});
				Tmp.from = i;

				if (!have_Pht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega(), runlog);
					for (size_t k = 0;k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega() - Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (!have_Flr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (size_t k = 0;k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (!have_Aug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (size_t k = 0;k < Auger.size(); k++)
						{
							if (Auger[k].val <= 0) continue;
							Tmp.val = Auger[k].val;
							Tmp.to = i - hole_posit[Auger[k].hole] + hole_posit[Auger[k].fill] + hole_posit[Auger[k].eject];
							Tmp.energy = Auger[k].energy;
							LocalAuger.push_back(Tmp);
						}
					}
				}
			}

			#pragma omp critical
			{
				Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
				FF.insert(FF.end(), LocalFF.begin(), LocalFF.end());
			}
		}

		sort(Store.Photo.begin(), Store.Photo.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(FF.begin(), FF.end(), [](ffactor A, ffactor B) { return (A.index < B.index); });
		GenerateRateKeys(Store.Auger);

		if (!have_Pht) {
			string dummy = RateLocation + "Photo.txt";
			RateData::WriteRates(dummy, Store.Photo);
		}
		if (!have_Flr) {
			string dummy = RateLocation + "Fluor.txt";
			RateData::WriteRates(dummy, Store.Fluor);
		}
		if (!have_Pht) {
			string dummy = RateLocation + "Auger.txt";
			RateData::WriteRates(dummy, Store.Auger);
		}
		if (saveFF) {
			string dummy = RateLocation + "Form_Factor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& ff : FF) {
				for (size_t i = 0;i < ff.val.size(); i++) fprintf(fl, "%3.5f ", ff.val[i]);
				fprintf(fl, "\n");
			}
			fclose(fl);
		}
	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	ofstream config_out(IndexTrslt);
	for (size_t i = 0;i < Index.size(); i++) {
		for (size_t j = 0;j < Max_occ.size(); j++) {
			config_out << Max_occ[j] - Index[i][j] << " ";
		}
		config_out << endl;
	}

 	return dimension;
}

// Called for molecular inputs.
// Computes molecular collision parameters.
RateData::Atom ComputeRateParam::SolvePlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, vector<bool> shell_check, ofstream & runlog)
{
	// Uses BEB model to compute fundamental
	// EII, Auger, Photoionisation and Fluorescence rates
	// Final_occ defines the lowest possible occupancies for the initial orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return Store;

	Store.num_conf = dimension;

	// Rates are always computed in-memory here; file IO is owned by atomic_rate_data
	// (the Stage-1 binary), which serializes the returned Store to HDF5. The legacy
	// have_* caching flags are retained as always-false so the compute block runs.
	const bool have_Aug = false, have_EII = false, have_Pht = false, have_Flr = false;
	const bool saveFF = true;

	if (!have_Aug || !have_EII || !have_Pht || !have_Flr || saveFF)
	{
		cout <<"======================================================="<<endl;
		cout << "Total number of configurations: " << dimension << endl;
		cout <<"Beginning Hartree-Fock BEB calculations for missing parameters " <<endl;
		cout <<"======================================================="<<endl;
		RateData::Rate Tmp;
		vector<RateData::Rate> LocalPhoto(0);
		vector<RateData::Rate> LocalFluor(0);
		vector<RateData::Rate> LocalAuger(0);
		vector<ffactor> LocalFF(0);
		vector<energy_config> LocalEnergyConfig(0); 
		// Electron impact ionization orbital enerrgy storage.
		RateData::EIIdata tmpEIIparams;
		int MaxBindInd = 0;
		// Slippery assumption - electron impact cannot ionize more than the XFEL photon.
		// Guard the scan: if every orbital is frozen (nothing is photoionizable, e.g. the
		// photon energy is below even the valence binding energy) this would otherwise run
		// off the end of Final_occ/orbitals.
		while (MaxBindInd < (int)orbitals.size()
		       && Final_occ[MaxBindInd] == orbitals[MaxBindInd].occupancy()) MaxBindInd++;
		tmpEIIparams.kin.clear();
		tmpEIIparams.kin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.ionB.clear();
		tmpEIIparams.ionB.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.fin.clear();
		tmpEIIparams.fin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.occ.clear();
		tmpEIIparams.occ.resize(orbitals.size() - MaxBindInd, 0);
		vector<RateData::EIIdata> LocalEIIparams(0);

		density.clear();

	  	#pragma omp parallel default(none) num_threads(input.Num_Threads())\
		shared(cout, runlog, shell_check, MaxBindInd, have_Aug, have_Flr, have_Pht, saveFF) \
		private(Tmp, LocalPhoto, LocalAuger, LocalFluor, LocalEnergyConfig, LocalEIIparams, tmpEIIparams, LocalFF) \
		firstprivate(Max_occ)  // I believe this should be shared, but it was in private() before (which I believe is a mistake, since this meant the value of max_occ was lost) so I'm putting it here to be safe. -S.P.
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0;i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				cout << "[HF BEB] configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (size_t j = 0;j < Orbitals.size(); j++) {
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
					// Store shell flag if orbital corresponds to shell
					if(shell_check[j]) Orbitals[j].flag_shell(true);
				}
				// Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				// Change Lattice to lattice for electron density evaluation.
				Potential U(&lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(lattice, Orbitals, U, input, runlog);

				// EII parameters to store for Later BEB model calculation.
				tmpEIIparams.init = i;
				int size = 0;
				for (int n = MaxBindInd; n < static_cast<int>(Orbitals.size()); n++) if (Orbitals[n].occupancy() != 0) size++;
				tmpEIIparams.kin = U.Get_Kinetic(Orbitals, MaxBindInd);
				tmpEIIparams.ionB = vector<float>(size, 0);
				tmpEIIparams.fin = vector<int>(size, 0);
				tmpEIIparams.occ = vector<int>(size, 0);
				size = 0;
				//tmpEIIparams.inds.resize(tmpEIIparams.vec2.size(), 0);
				for (int j = MaxBindInd; j < static_cast<int>(Orbitals.size()); j++) {
					if (Orbitals[j].occupancy() == 0) continue;
					int old_occ = Orbitals[j].occupancy();
					Orbitals[j].set_occupancy(old_occ - 1);
					tmpEIIparams.fin[size] = mapOccInd(Orbitals);
					tmpEIIparams.occ[size] = old_occ;
					Orbitals[j].set_occupancy(old_occ);
					tmpEIIparams.ionB[size] = float(-1*Orbitals[j].Energy);
					tmpEIIparams.kin[size] /= tmpEIIparams.ionB[size];
					size++;
				}
				LocalEIIparams.push_back(tmpEIIparams);

				bool calc_bound_transport = true;
				if (calc_bound_transport){
					assert(Max_occ.size() == Orbitals.size());
					size = 0;
					double valence_energy = 0;
					for (int j = MaxBindInd; j < static_cast<int>(Orbitals.size()); j++) {
						if (Orbitals[j].occupancy() == 0) continue;
						valence_energy = -tmpEIIparams.ionB[size];
						size++;
					}
					int receiver_orbital = -1;
					int donator_orbital = -1;
					int receiver_occupancy= -1;
					int donator_occupancy = -1;
					for (int j = MaxBindInd; j < static_cast<int>(Orbitals.size()); j++) {
						if (Orbitals[j].occupancy() == 0) continue;
						donator_orbital = j;
						receiver_orbital = j;
					}				
					if (Max_occ[receiver_orbital] == Orbitals[receiver_orbital].occupancy())
						receiver_orbital += 1;
					int donator_idx;
					int receiver_idx;
					// Find index of config after electron is transported to this config
					if (receiver_orbital >= static_cast<int>(Orbitals.size()))
						// maximum orbital is filled - Transport is disallowed.
						receiver_idx = -1;
					else{
						int old_occ = Orbitals[receiver_orbital].occupancy();
						Orbitals[receiver_orbital].set_occupancy(old_occ +1);
						receiver_idx = mapOccInd(Orbitals);
						Orbitals[receiver_orbital].set_occupancy(old_occ);
					}
					// Find index of config after electron is transported away.
					if (donator_orbital == -1){
						// Empty shells - Transport is disallowed (I don't think such a config would be passed in this loop anyway).
						donator_idx = -1;
					}
					else{
						int old_occ = Orbitals[donator_orbital].occupancy();
						Orbitals[donator_orbital].set_occupancy(old_occ - 1);
						donator_idx = mapOccInd(Orbitals);
						Orbitals[donator_orbital].set_occupancy(old_occ);
					}
					assert(valence_energy < 0);
					LocalEnergyConfig.push_back(energy_config{(int)i,valence_energy,receiver_idx,donator_idx});  // if valence energy of atom 1 at receiver_idx + valence energy of atom 2 at donator_idx  << current valence energies, transport occurs.
				}
				DecayRates Transit(lattice, Orbitals, u, input);

				if (saveFF) LocalFF.push_back({i, Transit.FT_density()});
				Tmp.from = i;

				if (!have_Pht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega(), runlog);
					for (size_t k = 0;k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega() + Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (!have_Flr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (size_t k = 0;k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (!have_Aug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (size_t k = 0;k < Auger.size(); k++)
						{
							if (Auger[k].val <= 0) continue;
							Tmp.val = Auger[k].val;
							Tmp.to = i - hole_posit[Auger[k].hole] + hole_posit[Auger[k].fill] + hole_posit[Auger[k].eject];
							Tmp.energy = Auger[k].energy;
							LocalAuger.push_back(Tmp);
						}
					}
				}
			}

			#pragma omp critical
			{
				Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
				Store.EIIparams.insert(Store.EIIparams.end(), LocalEIIparams.begin(), LocalEIIparams.end());
				FF.insert(FF.end(), LocalFF.begin(), LocalFF.end());
				Store.EnergyConfig.insert(Store.EnergyConfig.end(), LocalEnergyConfig.begin(), LocalEnergyConfig.end());
			}
		}

		sort(Store.Photo.begin(), Store.Photo.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(FF.begin(), FF.end(), [](ffactor A, ffactor B) { return (A.index < B.index); });
		sort(Store.EIIparams.begin(), Store.EIIparams.end(), [](RateData::EIIdata A, RateData::EIIdata B) {return (A.init < B.init);});
		GenerateRateKeys(Store.Auger);
	}

	// Populate human-readable configuration names (consumed downstream by ac4dc).
	Store.index_names.clear();
	for (size_t i = 0; i < Index.size(); i++) {
		Store.index_names.push_back(InterpretIndex(i));
	}

 	return Store;
}


int ComputeRateParam::Symbolic(const string & input, const string & output)
{
	if (Store.Photo.size() == 0)
	{
		if (exists_test(input)) {
			ifstream Rates_in(input);
			ofstream Rates_out(output);

			RateData::Rate Tmp;
			char type;

			while (!Rates_in.eof())
			{
				string line;
				getline(Rates_in, line);

				stringstream stream(line);
				stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

				Rates_out << Tmp.val << " " << InterpretIndex(Tmp.from) << " " << InterpretIndex(Tmp.to) << " " << Tmp.energy << endl;
			}

			Rates_in.close();
			Rates_out.close();

			return 0;
		}
		else return 1;
	} else {
		ofstream Rates_out(output);

		for (auto& v : Store.Photo) {
			Rates_out << v.val << " " << InterpretIndex(v.from) << " " << InterpretIndex(v.to) << " " << v.energy << endl;
		}

		Rates_out.close();

		return 0;
	}

}
string ComputeRateParam::InterpretIndex(int i)
{
	// Outputs electronic configuration referenced in the i^th entry of
	// EIIparams, Auger, Photo, Fluor
	// LaTeX output format for direct insertion into table.
	string Result;
	if (!Index.empty()) {
		if (Index[i].size() == orbitals.size())	{
			ostringstream tmpstr;
			tmpstr << '$';
			for (size_t j = 0;j < orbitals.size(); j++) {
				tmpstr << orbitals[j].N();
				switch (orbitals[j].L()) {
				case 0:
					tmpstr << 's';
					break;
				case 1:
					tmpstr << 'p';
					break;
				case 2:
					tmpstr << 'd';
					break;
				case 3:
					tmpstr << 'f';
					break;
				default:
					tmpstr << 'x';
					break;
				}
				tmpstr << "^{" << orbitals[j].occupancy() - Index[i][j] <<'}';
			}
			tmpstr << '$';
			Result = tmpstr.str();
		}
	}
	return Result;
}

int ComputeRateParam::Charge(int Iconf)
{
	int Result = 0;
	for (size_t j = 0;j < Index[Iconf].size(); j++) Result += Index[Iconf][j] - Index[0][j];
	return Result;
}

bool ComputeRateParam::SetupIndex(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	if (orbitals.size() != Final_occ.size())
	{
		runlog << "Final occupancies should be provided for all orbitals." << endl;
		return false;
	}
	// Work out the number of allowed configurations
	dimension = 1;
	int orbitals_size = orbitals.size();
	int l_hole = 0;//lowest energy orbital with allowed holes
	hole_posit.clear();
	hole_posit.resize(orbitals.size());

	vector<int> max_holes(orbitals.size(), 0);
	for (int i = orbitals_size - 1; i >= 0; i--)
	{
		max_holes[i] = Max_occ[i] - Final_occ[i] + 1;
		if (Max_occ[i] == Final_occ[i]) l_hole++;
		hole_posit[i] = dimension;
		dimension *= max_holes[i];
	}

	//configuration information. Index[i][j] denotes number of holes in oorbital [j], that together form configuration [i].
	Index.resize(dimension);
	for (auto& v : Index) {
		v.resize(orbitals.size());
	}

	for (int i = 0;i < dimension; i++)
	{
		int tmp = i;
		for (int j = 0;j < orbitals_size; j++)
		{
			Index[i][j] = tmp / hole_posit[j];
			if (Index[i][j] > Max_occ[j]) { Index[i][j] = Max_occ[j]; }
			tmp -= Index[i][j] * hole_posit[j];
		}
	}
	return true;
}

ComputeRateParam::~ComputeRateParam()
{
}

vector<double> ComputeRateParam::generate_dT(int num_elem)//default time interval
{
	vector<double> Result(num_elem, 0);
	double tmp = 1;
	for (int i = 0;i < num_elem; i++) {
	tmp = fabs(1.*i / (num_elem-1) - 0.5) + 0.01;
		Result[i] = tmp;
	}
	return Result;
}

vector<double> ComputeRateParam::generate_T(vector<double>& dT)//default time
{
	vector<double> Result(dT.size(), 0);
	vector<double> Bashforth_4{ 55. / 24., -59. / 24., 37. / 24., -9. / 24. }; //Adams�Bashforth method
	for (int i = 1; i < static_cast<int>(Bashforth_4.size()); i++)//initial few points
	{
		Result[i] = Result[i - 1] + dT[i - 1];
	}
	for (int i = Bashforth_4.size(); i < static_cast<int>(Result.size()); i++)//subsequent points
	{
		Result[i] = Result[i - 1];
		for (size_t j = 0;j < Bashforth_4.size(); j++)
		{
			Result[i] += Bashforth_4[j] * dT[i - j - 1];
		}
	}

	return Result;
}

vector<double> ComputeRateParam::generate_I(vector<double>& Time, double Fluence, double Sigma)//intensity of the Gaussian X-ray pulse
{

	vector<double> Result(Time.size(), 0);

	double midpoint = 0.5*(T.back() + T[0]);
	double denom = 2*Sigma*Sigma;
	double norm = 1./sqrt(denom*Constant::Pi);
	//include the window function to make it exactly 0 at the beginning and smoothly increase toward Gaussian


  int smooth = T.size()/10;
  double tmp = 0;
	for (int i = 0;i < static_cast<int>(Time.size()); i++)
	{
    Result[i] = Fluence * norm * exp(-(Time[i] - midpoint)*(Time[i] - midpoint) / denom);
    if (i < smooth) {
      tmp = fabs(T[i] - T[0]) / fabs(T[smooth] - T[0]);
      Result[i] *= tmp*tmp*(3 - 2 * tmp);
    }
  }

  /*
  for (size_t i = 0;i < Time.size(); i++)
	{
    Result[i] = Fluence / T.back();
	}
  */



	return Result;
}

int ComputeRateParam::extend_I(vector<double>& Intensity, double new_max_T, double step_T)
{
  // Uniform mesh is added.
  double last_T = T.back(), dT_last = dT.back();
  //double MidT = 0.5*(T.back() - T[0]);
  //double denom = 2*Sigma*Sigma;
  //double norm = 1./sqrt(denom*Constant::Pi);

  while (last_T < new_max_T) {
    dT.push_back(dT_last);
    T.push_back(last_T + dT_last);
    last_T = T.back();
    Intensity.push_back(0);
  }
  return 0;
}

void SmoothOrigin(vector<double> & T, vector<double> & F)
{
	int smooth = T.size() / 10;
	for (int i = 0;i < smooth; i++)
	{
		F[i] *= (T[i] / T[smooth])*(T[i] / T[smooth])*(3 - 2 * (T[i] / T[smooth]));
	}
}

/**
 * @brief Intensity profile normalized to 1.
 * Time is assumbed to be in FEM.
 * 
 * @return vector<double> 
 */
vector<double> ComputeRateParam::generate_G()
{
	double Sigma = input.Width()/(2*sqrt(2*log(2.)));

	return generate_I(T, 1, Sigma);
}

void ComputeRateParam::GenerateRateKeys(vector<RateData::Rate> & ToSort)
{
	int CurrentFrom = 0;
	int start = 0;
	RatesFromKeys.push_back(0);
	for (int i = 1; i < static_cast<int>(ToSort.size()); i++) {
		if (ToSort[i].from != CurrentFrom) {
			CurrentFrom = ToSort[i].from;
			start = RatesFromKeys.back();
			RatesFromKeys.push_back(i);
			sort(ToSort.begin() + start, ToSort.begin() + RatesFromKeys.back(), sortRatesTo);
		}
	}
}

int ComputeRateParam::mapOccInd(vector<RadialWF> & Orbitals)
{
	int Result = 0;
	for (size_t j = 0;j < hole_posit.size(); j++)	{
		Result += (orbitals[j].occupancy() - Orbitals[j].occupancy())*hole_posit[j];
	}

	return Result;
}

double ComputeRateParam::T_avg_RMS(vector<pair<double, int>> conf_RMS)
{
  // Calculate pulse-averaged root mean square radius of an atom.
  double tmp = 0;

  vector<double> intensity = generate_G();
  if (P.size()-1 != density.size()) return -1;
  for (int m = 0; m < static_cast<int>(T.size()); m++) {
    tmp = 0;
    for (size_t i = 0;i < conf_RMS.size(); i++) tmp += P[i][m]*conf_RMS[i].first;
    intensity[m] *= tmp;
  }

  Grid Time(T, dT);
  Adams I(Time, 10);

  return I.Integrate(&intensity, 0, T.size()-1);
}


double ComputeRateParam::T_avg_Charge()
{
  // Calculate pulse-averaged charge of atom.
  double tmp = 0;

  vector<double> intensity = generate_G();
  for (int m = 0; m < static_cast<int>(T.size()); m++) {
    tmp = 0;
    for (size_t i = 0;i < charge.size(); i++) tmp += (input.Nuclear_Z() - i)*charge[i][m];
    intensity[m] *= tmp;
  }

  Grid Time(T, dT);
  Adams I(Time, 10);

  return I.Integrate(&intensity, 0, T.size()-1);
}


// Unfinished or archived:
//
// int ComputeRateParam::SetupAndSolve(ofstream & runlog)
// {
// 	// initial conditions for rate equation
// 	// first index represents the configuration
// 	// second index represents the time. The routine will expand it itself unless
// 	// you provide first adams_n points yourself.
// 	double fluence = input.Fluence();
// 	double Sigma = input.Width()/ (2*sqrt(2*log(2.)));
// 	int T_size = input.TimePts();
// 	vector<double> InitCond(dimension, 0);
// 	InitCond[0] = 1;
// 	P.clear();
//
// 	vector<double> Intensity;
// 	double scaling_T = 1;
//
// 	int converged = 1;
// 	Plasma Mxwll(T.size());
// 	while (converged != 0)
// 	{
// 		dT.clear();
// 		dT = generate_dT(T_size);
// 		T.clear();
// 		T = generate_T(dT);
// 		scaling_T = 4 * input.Width() / T.back();
// 		for (size_t i = 0;i < T.size(); i++) {
// 			T[i]  *= scaling_T;
// 			dT[i] *= scaling_T;
// 		}
// 		Intensity = generate_I(T, fluence, Sigma);
// 		// SmoothOrigin(T, Intensity);
// 		IntegrateRateEquation Calc(dT, T, Store, InitCond, Intensity);
// 		converged = Calc.Solve(0, 1, input.Out_T_size());
// 		if (converged == 0) {
// 			cout << "[Rates] Final number of time steps: " << T_size << endl;
// 			P = Calc.GetP();
// 			T.clear();
// 			T = Calc.GetT();
// 			dT.clear();
// 			dT = vector<double>(T.size(), 0);
// 		    for (int m = 1; m < T.size(); m++) dT[m-1] = T[m] - T[m-1];
// 		    dT[T.size()-1] = dT[T.size()-2];
// 		} else {
// 			cout << "[Rates] Diverged at step: " << converged << " of " << T_size << endl;
// 			cout << "Halving timestep..." << endl;
// 			T_size *= 2;
// 		}
// 	}
//
//
// 	Intensity.clear();
// 	Intensity = generate_I(T, fluence, Sigma);
// 	SmoothOrigin(T, Intensity);
//
// 	int tmp = 0;
// 	for (size_t i = 0;i < Index.back().size(); i++) tmp += Index.back()[i] - Index.begin()->at(i);
// 	charge.clear();
// 	charge.resize(tmp + 1);
// 	for (auto& ch: charge) ch = vector<double>(T.size(), 0);
// 	for (size_t i = 0;i < P.size(); i++)
// 	{
// 		tmp = Charge(i);
// 		for (int m = 0; m < P[i].size(); m++) {
// 			charge[tmp][m] += P[i][m];
// 		}
// 	}
//
// 	double t_tmp = 0;
// 		for (int m = 0; m < T.size(); m++) {
// 		T[m] = (T[m]-0.5*T.back())*Constant::fs_per_au;
// 		dT[m] *= Constant::fs_per_au;
// 	}
//
// 	if (input.Write_Charges()) {
// 		cout << "Writing charges..."<<endl;
// 		string ChargeName = "./output/Charge_" + input.Name() + ".txt";
// 		ofstream charge_out(ChargeName);
// 		double chrg_tmp = 0;
// 		for (int m = 0; m < T.size(); m++) {
// 			for (size_t i = 0;i < charge.size(); i++) {
// 				chrg_tmp = charge[i][m];
// 				if (chrg_tmp <= 0.00000001) charge_out << 0 << " ";
// 				else charge_out << chrg_tmp << " ";
// 			}
// 			charge_out << T[m];
// 			if (m != T.size() - 1) charge_out << endl;
// 		}
//
// 		charge_out.close();
// 	}
//
// 	if (input.Write_Intensity()) {
// 		cout << "Writing intensity..."<<endl;
// 		string IntensityName = "./output/Intensity_" + input.Name() + ".txt";
// 		ofstream intensity_out(IntensityName);
//
// 		double I_max = *max_element(begin(Intensity), end(Intensity));
// 		for (int m = 0; m < T.size(); m++) {
// 			intensity_out << Intensity[m] / I_max << " " << T[m];
// 			if (m != T.size() - 1) intensity_out << endl;
// 		}
//
// 		intensity_out.close();
// 	}
//
//
//
// 	return 0;
// }
//
//
// int ComputeRateParam::SetupAndSolve(MolInp & Input, ofstream & runlog)
// {
// 	// initial conditions for rate equation
// 	// first index represents the configuration
// 	// second index represents the time. The routine will expand it itself unless
// 	// you provide first adams_n points yourself.
// 	double fluence = Input.Fluence();
// 	double Sigma = Input.Width()/ (2*sqrt(2*log(2.)));
// 	int T_size = Input.ini_T_size();
//
// 	P.clear();
//
// 	vector<double> Intensity;
// 	double scaling_T = 1;
//
// 	int converged = 1;
// 	Plasma Mxwll(T.size());
// 	while (converged != 0)
// 	{
// 		dT.clear();
// 		dT = generate_dT(T_size);
// 		T.clear();
// 		T = generate_T(dT);
// 		scaling_T = 4*input.Width() / T.back();
// 		for (size_t i = 0;i < T.size(); i++) {
//       //T[i] = T[i]-0.5*T.back();
// 			T[i] *= scaling_T;
// 			dT[i] *= scaling_T;
// 		}
// 		Intensity = generate_I(T, fluence, Sigma);
//
// 		//SmoothOrigin(T, Intensity);
//  		Mxwll.resize(T.size());
// 		IntegrateRateEquation Calc(dT, T, Input.Store, Mxwll, Intensity);
//
//     	T_size = T.size();
//
// 		converged = Calc.Solve(Mxwll, Input.Store, Input.Out_T_size());
//
// 		if (converged == 0)
// 		{
// 			cout << "Final number of time steps: " << T_size << endl;
// 			P = Calc.GetP();
// 			T.clear();
// 			T = Calc.GetT();
// 			dT.clear();
// 			dT = generate_dT(T.size());
// 		}
// 		else
// 		{
// 			cout << "Diverged at step: " << converged << " of " << T_size << endl;
// 			T_size *= 2;
// 		}
// 	}
//
// 	Intensity.clear();
// 	Intensity = generate_I(T, fluence, Sigma);
//
// 	double t_tmp = 0;
// 	double I_max = *max_element(begin(Intensity), end(Intensity));
//
// 	for (int m = 0; m < T.size(); m++) {
// 		T[m] = (T[m]-0.5*T.back())*Constant::fs_per_au;
// 		dT[m] *= Constant::fs_per_au;
// 	}
//
// 	int shift = 0;
// 	vector<int> P_to_charge(0);
//
//   // Aggregate and output charges, plasma parameters, and other parameters into an output.
//   // Charges.
//   vector<vector<double>> AllAtomCharge(Input.Atomic.size(), vector<double>(T.size(), 0));
//
// 	for (int a = 0; a < Input.Atomic.size(); a++) {
//
// 		// Occupancies associated with the atom "a".
// 		vector<vector<double*>> map_p(Input.Store[a].num_conf);
// 		for (size_t i = 0;i < Input.Store[a].num_conf; i++) {
// 			map_p[a].push_back(P[i + shift].data());
// 		}
//
// 		// Work out charges of all configurations.
// 		double chrg_tmp = 0;
// 		int tmp = Input.Index[a][0].size();
// 		P_to_charge.clear();
// 		P_to_charge.resize(Input.Store[a].num_conf, 0);
// 		for (size_t i = 0;i < Input.Index[a].size(); i++) {
// 			for (size_t j = 0;j < tmp; j++) P_to_charge[i] += Input.Index[a][i][j];
// 			chrg_tmp = 0;
// 		}
//
// 		charge.clear();
// 		charge.resize(*max_element(begin(P_to_charge), end(P_to_charge)) + 1);
// 		for (auto& ch: charge) ch = vector<double>(T.size(), 0);
// 		for (size_t i = 0;i < map_p[a].size(); i++)
// 		{
// 			tmp = P_to_charge[i];
// 			for (int m = 0; m < P[i].size(); m++) charge[tmp][m] += *(map_p[a][i] + m);
// 		}
//
//     for (int m = 0; m < T.size(); m++) {
// 			for (int i = 1; i < charge.size(); i++) {
//         AllAtomCharge[a][m] += i*charge[i][m];
//       }
//     }
//
// 		shift += Input.Store[a].num_conf;
//
// 		if (Input.Write_Charges()) {
// 			string ChargeName = "./output/Charge_" + Input.Store[a].name + ".txt";
// 			ofstream charge_out(ChargeName);
// 			double chrg_tmp = 0;
// 			for (int m = 0; m < T.size(); m++) {
// 				for (size_t i = 0;i < charge.size(); i++) {
// 					chrg_tmp = charge[i][m];
// 					if (chrg_tmp <= 0.00000001) charge_out << 0 << " ";
// 					else charge_out << chrg_tmp << " ";
// 				}
// 				charge_out << T[m];
// 				if (m != T.size() - 1) charge_out << endl;
// 			}
//
// 			charge_out.close();
// 		}
//
// 	}
//
// 	if (Input.Write_MD_data()) {
// 		string MD_Data = "./output/MD_Data.txt";
// 		ofstream OutFile(MD_Data);
// 		OutFile << T.size() << endl;
// 		OutFile << "Time ";
// 		for (int a = 0; a < Input.Atomic.size(); a++) OutFile << Input.Atomic[a].Nuclear_Z() << " ";
// 		OutFile << "N(elec) E(elec)";
// 		for (int m = 0; m < T.size(); m++) {
// 			OutFile << endl << T[m] << " ";
// 			for (int a = 0; a < AllAtomCharge.size(); a++) {
// 				OutFile << AllAtomCharge[a][m] << " ";
// 			}
// 			if (m != 0) OutFile << Mxwll.state[m].N << " " << Mxwll.state[m].E;
// 			else OutFile << 0 << " " << 0;
// 		}
//
// 		OutFile.close();
// 	}
//
// 	if (Input.Write_Intensity()) {
// 		string IntensityName = "./output/Intensity_" + Input.name + ".txt";
// 		ofstream intensity_out(IntensityName);
//
// 		double I_max = *max_element(begin(Intensity), end(Intensity));
// 		for (int m = 0; m < T.size(); m++) {
// 			intensity_out << Intensity[m] / I_max << " " << T[m];
// 			if (m != T.size() - 1) intensity_out << endl;
// 		}
//
// 		intensity_out.close();
// 	}
//
// 	return 0;
// }
