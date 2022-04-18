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
#include "ComputeRateParam.hpp"

inline bool exists_test(const std::string&);
// std::vector<double> generate_dT(int);
// std::vector<double> generate_T(std::vector<double>&);
// std::vector<double> generate_I(std::vector<double>&, double, double);
// void SmoothOrigin(std::vector<double>&, std::vector<double>&);


inline bool exists_test(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}



using namespace std;

// Called when atomic input is relevant.
int ComputeRateParam::SolveFrozen(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & runlog)
{
	// Calculates cross sections needed to solve rate equations.
	// Final_occ defines the lowest possible occupancies for the initiall orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return 1;

	// Make an output folder, if it doesn't exist yet
	if (!exists_test("./output")) {
		mkdir("output", ACCESSPERMS);
	}


	std::string RateLocation = "./output/" + input.Name() + "/Xsections/";

	// Make a data folder inside output
	if (!exists_test("./output/" + input.Name())) {
		std::string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	// Make a rate-containing folder inside the data folder
	if (!exists_test(RateLocation)) {
		mkdir(RateLocation.c_str(), ACCESSPERMS);
	}

	std::cout << "Check if there are pre-calculated rates..." << endl;
	bool have_Pht = RateData::ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool have_Flr = RateData::ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool have_Aug = RateData::ReadRates(RateLocation + "Auger.txt", Store.Auger);
	bool saveFF = true; // TODO: Make MolInp control this

	if (have_Pht) printf("Photoionization rates found. Reading...\n");
	if (have_Flr) printf("Fluorescence rates found. Reading...\n");
	if (have_Aug) printf("Auger rates found. Reading...\n");

	std::string PolarFileName = "./output/Polar_" + input.Name() + ".txt";

	if ( true || !have_Pht || !have_Flr || !have_Aug )
	{
		std::cout <<"======================================================="<<"\n";
		std::cout << "Total number of configurations: " << dimension << endl;
		std::cout <<" Beginning Hartree-Fock Frozen calculations... "<<"\n";
		std::cout <<"======================================================="<<"\n";
		RateData::Rate Tmp;
		std::vector<RateData::Rate> LocalPhoto(0);
		std::vector<RateData::Rate> LocalFluor(0);
		std::vector<RateData::Rate> LocalAuger(0);
		std::vector<ffactor> LocalFF(0);

		#pragma omp parallel default(none) num_threads(input.Num_Threads()) \
		shared(std::cout, runlog, have_Aug, have_Flr, have_Pht, saveFF) private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor, LocalFF)
		{
			#pragma omp for schedule(dynamic) nowait
			for (size_t i = 0;i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				std::vector<RadialWF> Orbitals = orbitals;
				std::cout << "[HF Frozen] configuration " << i << " thread " << omp_get_thread_num() << endl;
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
					std::vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega(), runlog);
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
						std::vector<fluor> Fluor = Transit.Fluor();
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
						std::vector<auger> Auger = Transit.Auger(Max_occ, runlog);
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
			std::string dummy = RateLocation + "Photo.txt";
			RateData::WriteRates(dummy, Store.Photo);
		}
		if (!have_Flr) {
			std::string dummy = RateLocation + "Fluor.txt";
			RateData::WriteRates(dummy, Store.Fluor);
		}
		if (!have_Pht) {
			std::string dummy = RateLocation + "Auger.txt";
			RateData::WriteRates(dummy, Store.Auger);
		}
		if (saveFF) {
			std::string dummy = RateLocation + "Form_Factor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& ff : FF) {
				for (size_t i = 0;i < ff.val.size(); i++) fprintf(fl, "%3.5f ", ff.val[i]);
				fprintf(fl, "\n");
			}
			fclose(fl);
		}
	}

	std::string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	std::ofstream config_out(IndexTrslt);
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
RateData::Atom ComputeRateParam::SolvePlasmaBEB(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & runlog)
{
	// Uses BEB model to compute fundamental
	// EII, Auger, Photoionisation and Fluorescence rates
	// Final_occ defines the lowest possible occupancies for the initial orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return Store;

	Store.num_conf = dimension;


	std::string RateLocation = "./output/" + input.Name() + "/Xsections/";
	if (!exists_test("./output/" + input.Name())) {
		std::string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}
	if (!exists_test(RateLocation)) {
		std::string dirstring = "output/" + input.Name() + "/Xsections";
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	bool have_Aug, have_EII, have_Pht, have_Flr;

	bool saveFF = true; //exists_test(RateLocation + "Form_Factor.txt");

	if (recalculate) {
		have_Aug=false;
		have_EII=false;
		have_Pht=false;
		have_Flr=false;
	} else {
		// Check if there are pre-calculated rates
		have_Pht = RateData::ReadRates(RateLocation + "Photo.txt", Store.Photo);
		have_Flr = RateData::ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
		have_Aug = RateData::ReadRates(RateLocation + "Auger.txt", Store.Auger);
		have_EII = RateData::ReadEIIParams(RateLocation + "EII.json", Store.EIIparams);
		std::cout <<"======================================================="<<"\n";
		std::cout <<"Seeking rates for atom "<< input.Name() <<"\n";
		if (have_Pht) std::cout<<"Photoionization rates found. Reading..." <<"\n";
		if (have_Flr) std::cout<<"Fluorescence rates found. Reading..."<<"\n";
		if (have_Aug) std::cout<<"Auger rates found. Reading..."<<"\n";
		if (have_EII) std::cout<<"EII Parameters found. Reading..."<<"\n";
	}

	if (!have_Aug || !have_EII || !have_Pht || !have_Flr || saveFF)
	{
		std::cout <<"======================================================="<<"\n";
		std::cout << "Total number of configurations: " << dimension << endl;
		std::cout <<"Beginning Hartree-Fock BEB calculations for missing parameters " <<"\n";
		std::cout <<"======================================================="<<"\n";
		RateData::Rate Tmp;
		std::vector<RateData::Rate> LocalPhoto(0);
		std::vector<RateData::Rate> LocalFluor(0);
		std::vector<RateData::Rate> LocalAuger(0);
		std::vector<ffactor> LocalFF(0);
		// Electron impact ionization orbital enerrgy storage.
		RateData::EIIdata tmpEIIparams;
		int MaxBindInd = 0;
		// Slippery assumption - electron impact cannot ionize more than the XFEL photon.
		while(Final_occ[MaxBindInd] == orbitals[MaxBindInd].occupancy()) MaxBindInd++;
		tmpEIIparams.kin.clear();
		tmpEIIparams.kin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.ionB.clear();
		tmpEIIparams.ionB.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.fin.clear();
		tmpEIIparams.fin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.occ.clear();
		tmpEIIparams.occ.resize(orbitals.size() - MaxBindInd, 0);
		std::vector<RateData::EIIdata> LocalEIIparams(0);

		density.clear();

	  	#pragma omp parallel default(none) num_threads(input.Num_Threads())\
		shared(std::cout, runlog, MaxBindInd, have_Aug, have_Flr, have_Pht, saveFF) \
		private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor, LocalEIIparams, tmpEIIparams, LocalFF)
		{
			#pragma omp for schedule(dynamic) nowait
			for (size_t i = 0;i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				std::vector<RadialWF> Orbitals = orbitals;
				std::cout << "[HF BEB] configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (size_t j = 0;j < Orbitals.size(); j++) {
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				// Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				// Change Lattice to lattice for electron density evaluation.
				Potential U(&lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(lattice, Orbitals, U, input, runlog);

				// EII parameters to store for Later BEB model calculation.
				tmpEIIparams.init = i;
				int size = 0;
				for (int n = MaxBindInd; n < Orbitals.size(); n++) if (Orbitals[n].occupancy() != 0) size++;
				tmpEIIparams.kin = U.Get_Kinetic(Orbitals, MaxBindInd);
				tmpEIIparams.ionB = std::vector<float>(size, 0);
				tmpEIIparams.fin = std::vector<int>(size, 0);
				tmpEIIparams.occ = std::vector<int>(size, 0);
				size = 0;
				//tmpEIIparams.inds.resize(tmpEIIparams.vec2.size(), 0);
				for (int j = MaxBindInd; j < Orbitals.size(); j++) {
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

				DecayRates Transit(lattice, Orbitals, u, input);

				if (saveFF) LocalFF.push_back({i, Transit.FT_density()});
				Tmp.from = i;

				if (!have_Pht) {
					std::vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega(), runlog);
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
						std::vector<fluor> Fluor = Transit.Fluor();
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
						std::vector<auger> Auger = Transit.Auger(Max_occ, runlog);
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
			}
		}

		sort(Store.Photo.begin(), Store.Photo.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](RateData::Rate A, RateData::Rate B) { return (A.from < B.from); });
		sort(FF.begin(), FF.end(), [](ffactor A, ffactor B) { return (A.index < B.index); });
		sort(Store.EIIparams.begin(), Store.EIIparams.end(), [](RateData::EIIdata A, RateData::EIIdata B) {return (A.init < B.init);});
		GenerateRateKeys(Store.Auger);

		// Write rates to file
		

		if (!have_Pht) {
			std::string dummy = RateLocation + "Photo.txt";
			std::cout<<"Saving photoionisation rates to "<<dummy<<"..."<<"\n";
			RateData::WriteRates(dummy, Store.Photo);
		}
		if (!have_Flr) {
			std::string dummy = RateLocation + "Fluor.txt";
			std::cout<<"Saving fluorescence rates to "<<dummy<<"..."<<"\n";
			RateData::WriteRates(dummy, Store.Fluor);
		}
		if (!have_Aug) {
			std::string dummy = RateLocation + "Auger.txt";
			std::cout<<"Saving Auger rates to "<<dummy<<"..."<<"\n";
			RateData::WriteRates(dummy, Store.Auger);
		}
		if (!have_EII) {
			std::string dummy = RateLocation + "EII.json";
			std::cout<<"Saving EII data to "<<dummy<<"..."<<"\n";
			RateData::WriteEIIParams(dummy, Store.EIIparams);
		}

		if (saveFF) {
			std::string dummy = RateLocation + "Form_Factor.txt";
			std::cout<<"Saving form factor data to "<<dummy<<"..."<<"\n";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& ff : FF) {
				for (size_t i = 0;i < ff.val.size(); i++) fprintf(fl, "%3.5f ", ff.val[i]);
				fprintf(fl, "\n");
			}
			fclose(fl);
		}
	}

	std::string IndexTrslt = "./output/" + input.Name() + "/index.txt";

	std::ofstream config_out(IndexTrslt);
	config_out<<"# idx | configuration"<<"\n";
	for (size_t i = 0; i < Index.size(); i++) {
		Store.index_names.push_back(InterpretIndex(i));
		config_out << i << " | " << Store.index_names.back();
		config_out<<"\n";
	}
	config_out.close();

 	return Store;
}


int ComputeRateParam::Symbolic(const std::string & input, const std::string & output)
{
	if (Store.Photo.size() == 0)
	{
		if (exists_test(input)) {
			std::ifstream Rates_in(input);
			std::ofstream Rates_out(output);

			RateData::Rate Tmp;
			char type;

			while (!Rates_in.eof())
			{
				std::string line;
				getline(Rates_in, line);

				std::stringstream stream(line);
				stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

				Rates_out << Tmp.val << " " << InterpretIndex(Tmp.from) << " " << InterpretIndex(Tmp.to) << " " << Tmp.energy << endl;
			}

			Rates_in.close();
			Rates_out.close();

			return 0;
		}
		else return 1;
	} else {
		std::ofstream Rates_out(output);

		for (auto& v : Store.Photo) {
			Rates_out << v.val << " " << InterpretIndex(v.from) << " " << InterpretIndex(v.to) << " " << v.energy << endl;
		}

		Rates_out.close();

		return 0;
	}

}


std::string ComputeRateParam::InterpretIndex(int i)
{
	// Outputs electronic configuration referenced in the i^th entry of
	// EIIparams, Auger, Photo, Fluor
	// LaTeX output format for direct insertion into table.
	std::string Result;
	if (!Index.empty()) {
		if (Index[i].size() == orbitals.size())	{
			std::ostringstream tmpstr;
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

bool ComputeRateParam::SetupIndex(std::vector<int> Max_occ, std::vector<int> Final_occ, std::ofstream & runlog)
{
	if (orbitals.size() != Final_occ.size())
	{
		runlog << "Final occupancies should be provided for all orbitals." << endl;
		return false;
	}
	if (Max_occ.size() != orbitals.size()) {
		Max_occ.resize(orbitals.size());
		for (size_t i = 0;i < orbitals.size(); i++) {
			Max_occ[i] = 4*orbitals[i].L() + 2;
		}
	}
	// Work out the number of allowed configurations
	dimension = 1;
	int orbitals_size = orbitals.size();
	int l_hole = 0;//lowest energy orbital with allowed holes
	hole_posit.clear();
	hole_posit.resize(orbitals.size());

	std::vector<int> max_holes(orbitals.size(), 0);
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

	for (size_t i = 0;i < dimension; i++)
	{
		int tmp = i;
		for (size_t j = 0;j < orbitals_size; j++)
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

// std::vector<double> ComputeRateParam::generate_dT(int num_elem)//default time interval
// {
// 	std::vector<double> Result(num_elem, 0);
// 	double tmp = 1;
// 	for (size_t i = 0;i < num_elem; i++) {
// 	tmp = fabs(1.*i / (num_elem-1) - 0.5) + 0.01;
// 		Result[i] = tmp;
// 	}
// 	return Result;
// }

// std::vector<double> ComputeRateParam::generate_T(std::vector<double>& dT)//default time
// {
// 	std::vector<double> Result(dT.size(), 0);
// 	std::vector<double> Bashforth_4{ 55. / 24., -59. / 24., 37. / 24., -9. / 24. }; //Adams�Bashforth method
// 	for (int i = 1; i < Bashforth_4.size(); i++)//initial few points
// 	{
// 		Result[i] = Result[i - 1] + dT[i - 1];
// 	}
// 	for (int i = Bashforth_4.size(); i < Result.size(); i++)//subsequent points
// 	{
// 		Result[i] = Result[i - 1];
// 		for (size_t j = 0;j < Bashforth_4.size(); j++)
// 		{
// 			Result[i] += Bashforth_4[j] * dT[i - j - 1];
// 		}
// 	}

// 	return Result;
// }

// std::vector<double> ComputeRateParam::generate_I(std::vector<double>& Time, double Fluence, double Sigma)//intensity of the Gaussian X-ray pulse
// {

// 	std::vector<double> Result(Time.size(), 0);

// 	double midpoint = 0.5*(T.back() + T[0]);
// 	double denom = 2*Sigma*Sigma;
// 	double norm = 1./sqrt(denom*Constant::Pi);
// 	//include the window function to make it exactly 0 at the beginning and smoothly increase toward Gaussian


//   int smooth = T.size()/10;
//   double tmp = 0;
// 	for (size_t i = 0;i < Time.size(); i++)
// 	{
//     Result[i] = Fluence * norm * exp(-(Time[i] - midpoint)*(Time[i] - midpoint) / denom);
//     if (i < smooth) {
//       tmp = fabs(T[i] - T[0]) / fabs(T[smooth] - T[0]);
//       Result[i] *= tmp*tmp*(3 - 2 * tmp);
//     }
//   }

//   /*
//   for (size_t i = 0;i < Time.size(); i++)
// 	{
//     Result[i] = Fluence / T.back();
// 	}
//   */



// 	return Result;
// }

// int ComputeRateParam::extend_I(std::vector<double>& Intensity, double new_max_T, double step_T)
// {
//   // Uniform mesh is added.
//   double last_T = T.back(), dT_last = dT.back();
//   //double MidT = 0.5*(T.back() - T[0]);
//   //double denom = 2*Sigma*Sigma;
//   //double norm = 1./sqrt(denom*Constant::Pi);

//   while (last_T < new_max_T) {
//     dT.push_back(dT_last);
//     T.push_back(last_T + dT_last);
//     last_T = T.back();
//     Intensity.push_back(0);
//   }
//   return 0;
// }

// void SmoothOrigin(std::vector<double> & T, std::vector<double> & F)
// {
// 	int smooth = T.size() / 10;
// 	for (size_t i = 0;i < smooth; i++)
// 	{
// 		F[i] *= (T[i] / T[smooth])*(T[i] / T[smooth])*(3 - 2 * (T[i] / T[smooth]));
// 	}
// }

// std::vector<double> ComputeRateParam::generate_G()
// {
// 	// Intensity profile normalized to 1.
// 	// Time is assumbed to be in FEM
// 	double Sigma = input.Width()/(2*sqrt(2*log(2.)));

// 	return generate_I(T, 1, Sigma);
// }

void ComputeRateParam::GenerateRateKeys(std::vector<RateData::Rate> & ToSort)
{
	int CurrentFrom = 0;
	int start = 0;
	RatesFromKeys.push_back(0);
	for (int i = 1; i < ToSort.size(); i++) {
		if (ToSort[i].from != CurrentFrom) {
			CurrentFrom = ToSort[i].from;
			start = RatesFromKeys.back();
			RatesFromKeys.push_back(i);
			sort(ToSort.begin() + start, ToSort.begin() + RatesFromKeys.back(), sortRatesTo);
		}
	}
}

int ComputeRateParam::mapOccInd(std::vector<RadialWF> & Orbitals)
{
	int Result = 0;
	for (size_t j = 0;j < hole_posit.size(); j++)	{
		Result += (orbitals[j].occupancy() - Orbitals[j].occupancy())*hole_posit[j];
	}

	return Result;
}

// double ComputeRateParam::T_avg_RMS(std::vector<pair<double, int>> conf_RMS)
// {
//   // Calculate pulse-averaged root mean square radius of an atom.
//   double tmp = 0;

//   std::vector<double> intensity = generate_G();
//   if (P.size()-1 != density.size()) return -1;
//   for (int m = 0; m < T.size(); m++) {
//     tmp = 0;
//     for (size_t i = 0;i < conf_RMS.size(); i++) tmp += P[i][m]*conf_RMS[i].first;
//     intensity[m] *= tmp;
//   }

//   Grid Time(T, dT);
//   Adams I(Time, 10);

//   return I.Integrate(&intensity, 0, T.size()-1);
// }


// double ComputeRateParam::T_avg_Charge()
// {
//   // Calculate pulse-averaged charge of atom.
//   double tmp = 0;

//   std::vector<double> intensity = generate_G();
//   for (int m = 0; m < T.size(); m++) {
//     tmp = 0;
//     for (size_t i = 0;i < charge.size(); i++) tmp += (input.Nuclear_Z() - i)*charge[i][m];
//     intensity[m] *= tmp;
//   }

//   Grid Time(T, dT);
//   Adams I(Time, 10);

//   return I.Integrate(&intensity, 0, T.size()-1);
// }
