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
#include "Constant.h"

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


// Called for molecular inputs.
// Computes molecular collision parameters.
RateData::Atom ComputeRateParam::SolveAtomicRatesAndPlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, vector<bool> shell_check, bool calculate_secondary_ionisation, ofstream & runlog)
{
	// Uses BEB model to compute fundamental
	// EII, Auger, Photoionisation and Fluorescence rates
	// Final_occ defines the lowest possible occupancies for the initial orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	const string PHOTO = "Photo.txt";
	const string AUGER = "Auger.txt";
	const string FLUOR = "Fluor.txt";
	const string EII =  "EII.json";

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return Store;

	Store.num_conf = dimension;


	string RateLocation = "./output/" + input.Name() + "/Xsections/";
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}
	if (!exists_test(RateLocation)) {
		string dirstring = "output/" + input.Name() + "/Xsections";
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	bool have_Aug, have_EII, have_Pht, have_Flr;

	if (recalculate) { // Hartree Fock is calculated once, at molinp photon energy 
		have_Aug=false;
		have_EII = (calculate_secondary_ionisation == false);
		have_Pht=false;
		have_Flr=false;
	} else { // First time: Save photoionization data for multiple photon energies. Second time: Interpolate from data.
		// Check if there are pre-calculated rates
		have_Pht = RateData::InterpolateRates(RateLocation, PHOTO, Store.Photo, input.Omega()); // Omega dependent
		have_Flr = RateData::ReadDecayRates(RateLocation, FLUOR, Store.Fluor,dimension);  // Dependent on ionizable shells
		have_Aug = RateData::ReadDecayRates(RateLocation, AUGER, Store.Auger,dimension);  // Dependent on ionizable shells
		have_EII = (calculate_secondary_ionisation == false);
		// Not sure if the below line will work properly, it would need to ensure that the energies of the knots are as expected. Not sure it does at present.
		//have_EII = RateData::ReadEIIParams(RateLocation + EII, Store.EIIparams) || (calculate_secondary_ionisation == false); // Dependent on the spline basis for electron distribution 
		
		cout <<"======================================================="<<endl;
		cout <<"Seeking rates for atom "<< input.Name() <<endl;
		if (have_Pht) cout<<"Photoionization rates found. Reading..." <<endl;
		if (have_Flr) cout<<"Fluorescence rates found. Reading..."<<endl;
		if (have_Aug) cout<<"Auger rates found. Reading..."<<endl;
		if (have_EII && calculate_secondary_ionisation) cout<<"EII Parameters found. Reading..."<<endl;
	}
	if (!have_Aug || !have_EII || !have_Pht || !have_Flr)
	{
		vector<vector<RateData::Rate>> PhotoArray(0);
		vector<vector<RateData::Rate>> LocalPhotoArray(0);
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
		while(Final_occ[MaxBindInd] == orbitals[MaxBindInd].occupancy()) MaxBindInd++;
		tmpEIIparams.kin.clear();
		tmpEIIparams.kin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.ionB.clear();
		tmpEIIparams.ionB.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.fin.clear();
		tmpEIIparams.fin.resize(orbitals.size() - MaxBindInd, 0);
		tmpEIIparams.occ.clear();
		tmpEIIparams.occ.resize(orbitals.size() - MaxBindInd, 0);
		vector<RateData::EIIdata> LocalEIIparams(0);

		// Omegas for which photoion. rates are calculated for.
		// Note this is in eV here!
		vector<double> photoion_omegas_to_save(0); 
		
		// TODO need to only calculate for energies with same allowed configs. 
		// // If LDA, it's much cheaper to calculate a full grid all at once. So do that.
		// if (input.Hamiltonian() == 1){
		// 	if (!recalculate){
		// 		double spacing = 500; 
		// 		double min_energy = 3000;
		// 		double max_energy = 18000;
		// 		for(double k = min_energy; k<=max_energy; k+=spacing){
		// 			photoion_omegas_to_save.push_back(k);
		// 		}
		// 	}
		// }
		//#endif
		photoion_omegas_to_save.push_back(input.Omega()*Constant::eV_per_Ha); // Also calculate for this run's omega.
		PhotoArray.resize(photoion_omegas_to_save.size());


		density.clear();

		// Convert to correct units of the sim
		for (size_t j = 0; j < photoion_omegas_to_save.size(); j++)
			photoion_omegas_to_save[j]/=Constant::eV_per_Ha;

	  	#pragma omp parallel default(none) num_threads(input.Num_Threads())\
		shared(cout, runlog, shell_check, MaxBindInd, have_Aug, have_Flr, have_Pht, have_EII, PhotoArray,photoion_omegas_to_save) \
		private(Tmp, LocalPhoto, LocalPhotoArray, LocalAuger, LocalFluor, LocalEnergyConfig, LocalEIIparams, tmpEIIparams, LocalFF) \
		firstprivate(Max_occ)  // I believe this should be shared, but it was in private() before (which I believe is a mistake, since this meant the value of Max_occ was lost) so I'm putting it here to be safe. -S.P. // Edit pretty sure this was unnecessary, Max_occ isn't set here.
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0;i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				LocalPhotoArray.resize(photoion_omegas_to_save.size());
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
				for (size_t n = MaxBindInd; n < Orbitals.size(); n++) if (Orbitals[n].occupancy() != 0) size++;
				tmpEIIparams.kin = U.Get_Kinetic(Orbitals, MaxBindInd);
				tmpEIIparams.ionB = vector<float>(size, 0);
				tmpEIIparams.fin = vector<int>(size, 0);
				tmpEIIparams.occ = vector<int>(size, 0);
				size = 0;
				//tmpEIIparams.inds.resize(tmpEIIparams.vec2.size(), 0);
				for (size_t j = MaxBindInd; j < Orbitals.size(); j++) {
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
					for (size_t j = MaxBindInd; j < Orbitals.size(); j++) {
						if (Orbitals[j].occupancy() == 0) continue;
						valence_energy = -tmpEIIparams.ionB[size];
						size++;
					}
					int receiver_orbital = -1;
					int donator_orbital = -1;
					for (size_t j = MaxBindInd; j < Orbitals.size(); j++) {
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

				LocalFF.push_back({i, Transit.FT_density()});
				Tmp.from = i;

				if (!have_Pht) {
					vector<vector<photo>> PhotoIon_array = Transit.Photo_Ion(photoion_omegas_to_save, runlog);
					for (size_t j = 0; j < photoion_omegas_to_save.size(); j++){
						vector<photo> PhotoIon = PhotoIon_array[j];
						for (size_t k = 0;k < PhotoIon.size(); k++)
						{
									if (PhotoIon[k].val <= 0) continue;
							Tmp.val = PhotoIon[k].val;
							Tmp.to = i + hole_posit[PhotoIon[k].hole];
							Tmp.energy = photoion_omegas_to_save[j] + Orbitals[PhotoIon[k].hole].Energy;
							LocalPhotoArray[j].push_back(Tmp);
						}
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
				for (size_t j = 0; j < LocalPhotoArray.size(); j++){
					PhotoArray[j].insert(PhotoArray[j].end(),LocalPhotoArray[j].begin(),LocalPhotoArray[j].end());
				}
				// if(!have_Pht){
				// 	Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				// }
				if(!have_Flr){
					Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				}
				if(!have_Aug){
					Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
				}
				if(!have_EII){
					Store.EIIparams.insert(Store.EIIparams.end(), LocalEIIparams.begin(), LocalEIIparams.end());
				}
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

		// Convert omegas back to eV
		for (size_t j = 0; j < photoion_omegas_to_save.size(); j++)
			photoion_omegas_to_save[j]*=Constant::eV_per_Ha;

		// Write rates to file
		
		// UNDERSCORES ARE USED FOR PARSING FILE NAMES.
		if (!have_Pht) {
			for (size_t j = 0; j < photoion_omegas_to_save.size(); j++){
				string dummy = RateLocation + std::to_string(photoion_omegas_to_save[j])+"_"+PHOTO;
				cout<<"Saving photoionisation rates to "<<dummy<<"..."<<endl;
				RateData::WriteRates(dummy, PhotoArray[j]);
			}
			assert(RateData::InterpolateRates(RateLocation, PHOTO, Store.Photo, input.Omega()));

		}
		if (!have_Flr) {
			string dummy = RateLocation + std::to_string(dimension)+"_"+FLUOR ;
			cout<<"Saving fluorescence rates to "<<dummy<<"..."<<endl;
			RateData::WriteRates(dummy, Store.Fluor);
		}
		if (!have_Aug) {
			string dummy = RateLocation +  std::to_string(dimension)+"_"+AUGER;
			cout<<"Saving Auger rates to "<<dummy<<"..."<<endl;
			RateData::WriteRates(dummy, Store.Auger);
		}
		if (!have_EII) {
			string dummy = RateLocation +EII;
			cout<<"Saving EII data to "<<dummy<<"..."<<endl;
			RateData::WriteEIIParams(dummy, Store.EIIparams);
		}


		string dummy = RateLocation + std::to_string(input.Omega()*Constant::eV_per_Ha) + "_FormFactor.txt";
		cout<<"Saving form factor data to "<<dummy<<"..."<<endl;
		FILE * fl = fopen(dummy.c_str(), "w");
		for (auto& ff : FF) {
			for (size_t i = 0;i < ff.val.size(); i++) fprintf(fl, "%3.5f ", ff.val[i]);
			fprintf(fl, "\n");
		}
		fclose(fl);
	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";

	ofstream config_out(IndexTrslt);
	config_out<<"# idx | configuration"<<endl;
	for (size_t i = 0; i < Index.size(); i++) {
		Store.index_names.push_back(InterpretIndex(i));
		config_out << i << " | " << Store.index_names.back();
		config_out<<endl;
	}
	config_out.close();

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
			//char type;

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
	for (size_t i = 1; i < Bashforth_4.size(); i++)//initial few points
	{
		Result[i] = Result[i - 1] + dT[i - 1];
	}
	for (size_t i = Bashforth_4.size(); i < Result.size(); i++)//subsequent points
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
	for (size_t i = 0;i < Time.size(); i++)
	{
    Result[i] = Fluence * norm * exp(-(Time[i] - midpoint)*(Time[i] - midpoint) / denom);
    if (static_cast<int>(i) < smooth) {
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
	for (int i = 0; i < smooth; i++)
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
	for (size_t i = 1; i < ToSort.size(); i++) {
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
  for (size_t m = 0; m < T.size(); m++) {
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
  for (size_t m = 0; m < T.size(); m++) {
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
