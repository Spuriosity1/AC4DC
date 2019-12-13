#include "RateEquationSolver.h"
#include "HartreeFock.h"
#include "DecayRates.h"
#include "Numerics.h"
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <dirent.h>
#include <sstream>
#include <string>
#include <omp.h>
#include <algorithm>
#include "EigenSolver.h"
#include "Plasma.h"

inline bool exists_test(const std::string&);
vector<double> generate_dT(int);
vector<double> generate_T(vector<double>&);
vector<double> generate_I(vector<double>&, double, double);
void SmoothOrigin(vector<double>&, vector<double>&);

using namespace CustomDataType;

RateEquationSolver::RateEquationSolver(Grid &Lattice, vector<RadialWF> &Orbitals, Potential &U, Input & Inp) : 
lattice(Lattice), orbitals(Orbitals), u(U), input(Inp)
{
//Orbitals are HF wavefunctions. This configuration is an initial state.
//Assuming there are no unoccupied states in initial configuration!!!
}

int RateEquationSolver::SolveFrozen(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	// Solves system of rate equations exactly.
	// Final_occ defines the lowest possible occupancies for the initiall orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return 1;

	cout << "Check if there are pre-calculated rates..." << endl;
	string RateLocation = "./output/" + input.Name() + "/Rates/";
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
		dirstring += "/Rates";
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}
	if (!exists_test("./output/" + input.Name() + "/Rates")) {
		string dirstring = "output/" + input.Name() + "/Rates";
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	bool existPht = !ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool existFlr = !ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool existAug = !ReadRates(RateLocation + "Auger.txt", Store.Auger);

	if (existPht) printf("Photoionization rates found. Reading...\n");
	if (existFlr) printf("Fluorescence rates found. Reading...\n");
	if (existAug) printf("Auger rates found. Reading...\n");

	// Electron density evaluation.
	density.clear();
	density.resize(dimension - 1);
	for (auto& dens: density) dens.resize(lattice.size(), 0.);

	if (existPht || existFlr || existAug)
	{
		cout << "No rates found. Calculating..." << endl;
		cout << "Total number of configurations: " << dimension << endl;
		Rate Tmp;
		vector<Rate> LocalPhoto(0);
		vector<Rate> LocalFluor(0);
		vector<Rate> LocalAuger(0);
		//ffactor ffTmp;//formfactor tmp
		//vector<ffactor> LocalFF(0);

		omp_set_num_threads(7);
		#pragma omp parallel default(none) \
		shared(cout, runlog, existAug, existFlr, existPht) private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0; i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				cout << "configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (int j = 0; j < Orbitals.size(); j++)
				{
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				//Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				// Change Lattice to lattice for electron density evaluation. 
				Potential U(&lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(lattice, Orbitals, U, input.Hamiltonian(), runlog);
				
				DecayRates Transit(lattice, Orbitals, u, input);

				Tmp.from = i;

				if (existPht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega()/Constant::eV_in_au, runlog);
					for (int k = 0; k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega()/Constant::eV_in_au - Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (existFlr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (int k = 0; k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (existAug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (int k = 0; k < Auger.size(); k++)
						{
							if (Auger[k].val <= 0) continue;
							Tmp.val = Auger[k].val;
							Tmp.to = i - hole_posit[Auger[k].hole] + hole_posit[Auger[k].fill] + hole_posit[Auger[k].eject];
							Tmp.energy = Auger[k].energy;
							LocalAuger.push_back(Tmp);
						}
					}
				}
				//FormFactor setup
				/*int infinity = 0;
				for (auto& it: Orbitals) {
					if (it.occupancy() == 0) continue;
					if (it.pract_infinity() > infinity) infinity = it.pract_infinity();
				}

				FormFactor W(&Lattice, Orbitals, infinity);
				ffTmp.index = i;
				ffTmp.val =  W.getAllFF();
				LocalFF.push_back(ffTmp);
        */

			}

			#pragma omp critical
			{
				Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
				//FF.insert(FF.end(), LocalFF.begin(), LocalFF.end());
			}
		}

		sort(Store.Photo.begin(), Store.Photo.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		GenerateRateKeys(Store.Auger);
		
		if (existPht) {
			string dummy = RateLocation + "Photo.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Photo) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFlr) {
			string dummy = RateLocation + "Fluor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Fluor) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existPht) {
			string dummy = RateLocation + "Auger.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Auger) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}

	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	ofstream config_out(IndexTrslt);
	for (int i = 0; i < Index.size(); i++) {
		config_out << i << " | ";
		for (int j = 0; j < Max_occ.size(); j++) {
			config_out << Max_occ[j] - Index[i][j] << " ";
		}
		config_out << endl;
	}

	if (0 != input.TimePts()) {
		SetupAndSolve(runlog);
	}
	else runlog << "Numbur of time points = 0. Skipping rate equation." << endl;

 	return dimension;
}

int RateEquationSolver::SolveFrozen(vector<RadialWF> & virt, vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	// Solves system of rate equations exactly.
	// Final_occ defines the lowest possible occupancies for the initiall orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return 1;

	cout << "Check if there are pre-calculated rates..." << endl;
	string RateLocation = "./output/" + input.Name() + "/Rates/";
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	bool existPht = !ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool existFlr = !ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool existAug = !ReadRates(RateLocation + "Auger.txt", Store.Auger);

	if (existPht) printf("Photoionization rates found. Reading...\n");
	if (existFlr) printf("Fluorescence rates found. Reading...\n");
	if (existAug) printf("Auger rates found. Reading...\n");

	string PolarFileName = "./output/Polar_" + input.Name() + ".txt";

	{
		cout << "No rates found. Calculating..." << endl;
		cout << "Total number of configurations: " << dimension << endl;
		Rate Tmp;
		vector<Rate> LocalPhoto(0);
		vector<Rate> LocalFluor(0);
		vector<Rate> LocalAuger(0);
		vector<polarize> LocalMixMe(0);

		omp_set_num_threads(7);
		#pragma omp parallel default(none) \
		shared(cout, virt, runlog, existAug, existFlr, existPht) private(Tmp, Max_occ, LocalPhoto, LocalAuger, LocalFluor,  LocalMixMe)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0; i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				vector<RadialWF> Virtual = virt;
				cout << "configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (int j = 0; j < Orbitals.size(); j++)
				{
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				Potential U(&Lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(Lattice, Orbitals, U, input.Hamiltonian(), runlog);
				HF.Get_Virtual(Virtual, Orbitals, U, runlog);
				
				LocalMixMe.push_back(HF.Hybrid(Orbitals,Virtual,U));
				LocalMixMe.back().index = i;

				DecayRates Transit(Lattice, Orbitals, u, input);

				Tmp.from = i;

				if (existPht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega()/Constant::eV_in_au, runlog);
					for (int k = 0; k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega()/Constant::eV_in_au - Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (existFlr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (int k = 0; k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (existAug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (int k = 0; k < Auger.size(); k++)
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
				//FF.insert(FF.end(), LocalFF.begin(), LocalFF.end());
			}
		}

		sort(Store.Photo.begin(), Store.Photo.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		GenerateRateKeys(Store.Auger);
		
		if (existPht) {
			string dummy = RateLocation + "Photo.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Photo) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFlr) {
			string dummy = RateLocation + "Fluor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Fluor) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existPht) {
			string dummy = RateLocation + "Auger.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Auger) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	ofstream config_out(IndexTrslt);
	for (int i = 0; i < Index.size(); i++) {
		for (int j = 0; j < Max_occ.size(); j++) {
			config_out << Max_occ[j] - Index[i][j] << " ";
		}
		config_out << endl;
	}

	if (0 != input.TimePts()) {
		SetupAndSolve(runlog);

		vector<double> Displace = PerturbMe(virt, 4, 0.00001);

		string DipoleOut = "./output/" + input.Name() + "/Ravg.txt";
		ofstream Dist(DipoleOut);
		for (int m = 0; m < T.size(); m++) {
			Dist << Displace[m] << " " << T[m] << endl;
		}
		Dist.close();
	}
	else runlog << "Numbur of time points = 0. Skipping rate equation." << endl;

 	return dimension;
}

AtomRateData RateEquationSolver::SolvePlasmaBEB(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog, vector<bool> args)
{
	// Solves system of Atomic rate equations + Temperature and electron concentration in Plasma exactly.
	// Final_occ defines the lowest possible occupancies for the initiall orbital.
	// Intermediate orbitals are recalculated to obtain the corresponding rates.

  // Example args:    vector<bool> args = {Init.Calc_R_ion(), Init.Calc_Pol_ion(), Init.Calc_FF_ion()};

	if (!SetupIndex(Max_occ, Final_occ, runlog)) return Store;

	Store.num_conf = dimension;

	cout << "Check if there are pre-calculated rates..." << endl;
	string RateLocation = "./output/" + input.Name() + "/Rates/";
	if (!exists_test("./output/" + input.Name())) {
		string dirstring = "output/" + input.Name();
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}
	if (!exists_test(RateLocation)) {
		string dirstring = "output/" + input.Name() + "/Rates";
		mkdir(dirstring.c_str(), ACCESSPERMS);
	}

	bool existPht = !ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool existFlr = !ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool existAug = !ReadRates(RateLocation + "Auger.txt", Store.Auger);

	if (existPht) printf("Photoionization rates found. Reading...\n");
	if (existFlr) printf("Fluorescence rates found. Reading...\n");
	if (existAug) printf("Auger rates found. Reading...\n");

	if (true)// EII parameters are not currently stored.
	{
		cout << "No rates found. Calculating..." << endl;
		cout << "Total number of configurations: " << dimension << endl;
		Rate Tmp;
		vector<Rate> LocalPhoto(0);
		vector<Rate> LocalFluor(0);
		vector<Rate> LocalAuger(0);
		// Electron impact ionization orbital enerrgy storage.
		EIIdata tmpEIIparams;
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
		vector<EIIdata> LocalEIIparams(0);

		density.clear();
		vector<vector<double>> Loc_dens;

		omp_set_num_threads(7);
	  #pragma omp parallel default(none) \
		shared(cout, runlog, MaxBindInd, existAug, existFlr, existPht, args) \
		private(Tmp, Max_occ,  LocalPhoto, LocalAuger, LocalFluor, LocalEIIparams, Loc_dens, tmpEIIparams)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int i = 0; i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
			{
				vector<RadialWF> Orbitals = orbitals;
				cout << "configuration " << i << " thread " << omp_get_thread_num() << endl;
				int N_elec = 0;
				for (int j = 0; j < Orbitals.size(); j++) {
					Orbitals[j].set_occupancy(orbitals[j].occupancy() - Index[i][j]);
					N_elec += Orbitals[j].occupancy();
				}
				//Grid Lattice(lattice.size(), lattice.R(0), lattice.R(lattice.size() - 1) / (0.3*(u.NuclCharge() - N_elec) + 1), 4);
				// Change Lattice to lattice for electron density evaluation. 
				Potential U(&lattice, u.NuclCharge(), u.Type());
				HartreeFock HF(lattice, Orbitals, U, input.Hamiltonian(), runlog);

				if (args[0] || args[1]) {
					Loc_dens.push_back(U.make_density(Orbitals) );
					Loc_dens.back().back() = i; // Last element is used for sorting, set to 0 afterwards.
				}

				// EII parameters to store for Later BEB model calculation.
				tmpEIIparams.init = i;
				int size = 0;
				for (int n = MaxBindInd; n < Orbitals.size(); n++) if (Orbitals[n].occupancy() != 0) size++;
				tmpEIIparams.kin = U.Get_Kinetic(Orbitals, MaxBindInd);
				tmpEIIparams.ionB = vector<float>(size, 0);
				tmpEIIparams.fin = vector<int>(size, 0);
				tmpEIIparams.occ = vector<int>(size, 0);
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

				Tmp.from = i;

				if (existPht) {
					vector<photo> PhotoIon = Transit.Photo_Ion(input.Omega()/Constant::eV_in_au, runlog);
					for (int k = 0; k < PhotoIon.size(); k++)
					{
						if (PhotoIon[k].val <= 0) continue;
						Tmp.val = PhotoIon[k].val;
						Tmp.to = i + hole_posit[PhotoIon[k].hole];
						Tmp.energy = input.Omega()/Constant::eV_in_au + Orbitals[PhotoIon[k].hole].Energy;
						LocalPhoto.push_back(Tmp);
					}
				}

				if (i != 0)
				{
					if (existFlr) {
						vector<fluor> Fluor = Transit.Fluor();
						for (int k = 0; k < Fluor.size(); k++)
						{
							if (Fluor[k].val <= 0) continue;
							Tmp.val = Fluor[k].val;
							Tmp.to = i - hole_posit[Fluor[k].hole] + hole_posit[Fluor[k].fill];
							Tmp.energy = Orbitals[Fluor[k].fill].Energy - Orbitals[Fluor[k].hole].Energy;
							LocalFluor.push_back(Tmp);
						}
					}

					if (existAug) {
						vector<auger> Auger = Transit.Auger(Max_occ, runlog);
						for (int k = 0; k < Auger.size(); k++)
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
				if (Loc_dens.size() != 0) density.insert(density.end(), Loc_dens.begin(), Loc_dens.end());
			}
		}

		sort(Store.Photo.begin(), Store.Photo.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Auger.begin(), Store.Auger.end(), [](Rate A, Rate B) { return (A.from < B.from); });
		sort(Store.Fluor.begin(), Store.Fluor.end(), [](Rate A, Rate B) { return (A.from < B.from); });
    sort(density.begin(), density.end(), [](vector<double> A, vector<double> B) { return (A.back() < B.back()); });
		sort(Store.EIIparams.begin(), Store.EIIparams.end(), [](EIIdata A, EIIdata B) {return (A.init < B.init);});
		GenerateRateKeys(Store.Auger);
				
		if (existPht) {
			string dummy = RateLocation + "Photo.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Photo) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existFlr) {
			string dummy = RateLocation + "Fluor.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Fluor) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
		if (existPht) {
			string dummy = RateLocation + "Auger.txt";
			FILE * fl = fopen(dummy.c_str(), "w");
			for (auto& R : Store.Auger) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
			fclose(fl);
		}
	}

	string IndexTrslt = "./output/" + input.Name() + "/index.txt";
	ofstream config_out(IndexTrslt);
	for (int i = 0; i < Index.size(); i++) {
		config_out << i << " | ";
		for (int j = 0; j < Max_occ.size(); j++) {
			config_out << Max_occ[j] - Index[i][j] << " ";
		}
		config_out << endl;
	}

  // Ionic radii - calculate and store.
  if (args[0]) {
    int fin = lattice.size() - 2; // Last value (integer) is used for sorting in the original densities.
    vector<double> integrand(lattice.size(), 0);
    Adams NI(lattice, 6);
    AtomAuxStore.r_ion.clear();
    AtomAuxStore.r_ion.resize(density.size());

    for (int j = 0; j < density.size(); j++) {
      fin = lattice.size() - 2;
      while (density[j][fin] == 0) fin--;
      for (int i = 0; i < lattice.size(); i++) integrand[i] = lattice.R(i) * density[j][i];

      AtomAuxStore.r_ion[j] = NI.Integrate(&integrand, 0, fin)/NI.Integrate(&density[j], 0, fin);
    }
  }



  // Form-factors - calculate and store.
  if (args[2]) {
    int fin = lattice.size() - 2; // Last value (integer) is used for sorting in the original densities.
    vector<double> integrand(lattice.size(), 0);

		ffactor FF_tmp;
		vector<ffactor> LocalFF(0);
    vector<ffactor> GlobalFF(0);

    #pragma omp parallel default(none) \
		shared(GlobalFF) \
		private(fin, LocalFF, FF_tmp)
		{
			#pragma omp for schedule(dynamic) nowait
			for (int n = 0; n < density.size(); n++)// Calculate form-factor for atomic densities in parallel.
			{
        fin = lattice.size() - 2;
        while (density[n][fin] == 0) fin--;

        FormFactor W(lattice, density[n], fin);
        FF_tmp.index = n;
        FF_tmp.val = W.getAllFF();
        LocalFF.push_back(FF_tmp);
      }

      #pragma omp critical
			{
				GlobalFF.insert(GlobalFF.end(), LocalFF.begin(), LocalFF.end());
			}
    }

    sort(GlobalFF.begin(), GlobalFF.end(), [](ffactor A, ffactor B) { return (A.index < B.index); });
    
    AtomAuxStore.form_fact_ion.clear();
    AtomAuxStore.form_fact_ion.resize(density.size());

    for (int n = 0; n < density.size(); n++) AtomAuxStore.form_fact_ion[n] = GlobalFF[n].val;
  }


 	return Store;
}



bool RateEquationSolver::ReadRates(const string & input, vector<Rate> & PutHere)
{
	if (exists_test(input))
	{
		Rate Tmp;
		ifstream infile;
		infile.open(input);

		char type;
		while (!infile.eof())
		{
			string line;
			getline(infile, line);

			stringstream stream(line);
			stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

			PutHere.push_back(Tmp);
		}

		infile.close();
		return true;
	}
	else return false;
}

int RateEquationSolver::Symbolic(const string & input, const string & output)
{
	if (Store.Photo.size() == 0)
	{
		if (exists_test(input)) {
			ifstream Rates_in(input);
			ofstream Rates_out(output);

			Rate Tmp;
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

vector<double> RateEquationSolver::PerturbMe(vector<RadialWF> & Virtual, double Dist, double Einit = 0)
{
	// Two identical atoms are separated by Dist (in a.u.).
	// As they ionize, they create time dependent electric field
	// that distorts their electronic configurations.
	// This functions calculates the field at the nucleus of atom 1
	// created by the atom 2. The field is assumed to be homogenous, and a 
	// small version of CI is calculated. 
	// The function returns a position of charge centroid relative to the nucleus as a function of time.
	vector<double> Result(T.size(), 0);
	vector<double> Efield(T.size(), 0);
	vector<vector<double>> PrevTVec(0);// Reference eigenvector to keep track of ground state configuration.
	vector<RadialWF*> OrbList(0);
	for (auto& Orb: orbitals) OrbList.push_back(&Orb);
	for (auto& Vrt: Virtual) OrbList.push_back(&Vrt);

	// Calculate time dependent Efield at nucleus position Dist. For now it is fixed thorugh the pulse.
	double tmp = 0;
	for (int m = 0; m < T.size(); m++) {
		tmp = 0;
		for (int Z = 0; Z < charge.size(); Z++) {
			tmp -= Z*charge[Z][m];
		}
		Efield[m] = tmp/Dist/Dist - Einit;
	}
	// For each configuration i run the configuration mixing for given Efield[m].
	// Multiply by P[i][m] and store.
	double wght = 0; // Weight factor accounts for fractional occupancies.
	int g = 0;// Ground state orbital to be mixed.
	int e = 0;// Excited state orbital to be mixed.
	for (int m = 0; m < T.size(); m++) {
		for (int N = 0; N < MixMe.size(); N++) {
			int L = MixMe[N].excited.size()+1;
			vector<vector<double>> Fock(L, vector<double>(L, 0));
			vector<vector<double>> Overlap(L, vector<double>(L, 0));

			Fock[0][0] = MixMe[N].refEnergy;
			Overlap[0][0] = 1;
			for (int i = 1; i < L; i++) {
				g = MixMe[N].excited[i-1][0];
				e = MixMe[N].excited[i-1][1];
				/*wght = sqrt(0.5*MixMe[N].reference[g]*(4*OrbList[e]->L() + 2 - OrbList[e]->occupancy())
					   /(2*OrbList[e]->L() + 1));//(2*OrbList[g]->L() + 1)*/
				// wght is adjusted to match polarizability calculation in ground state.
				wght = 2*OrbList[g]->occupancy()*(1 - 1.*OrbList[e]->occupancy()/(4*OrbList[e]->L() + 2))/3/(2*OrbList[g]->L() + 1);
				wght = sqrt(wght);
				Overlap[i][i] = 1;
				Fock[0][i] = -Efield[m]*wght*MixMe[N].Dipoles[i-1];
				Fock[i][0] = Fock[0][i];
				Fock[i][i] = MixMe[N].extEnergy[i-1];
			}
			if (m == 0 && MixMe[N].index == 0) {
				double Pol = 0;
				for (int i = 1; i < L; i++) {
					Pol -= Fock[0][i]*Fock[0][i]/(Fock[0][0] - Fock[i][i]);
				}
				cout << "Polarizability estimate: " << Pol/Efield[0]/Efield[0] << endl;
			}

			EigenSolver W;
			W.SolveGenEig(Fock, Overlap);

			vector<double> Vals(W.EigenVals());
			vector<vector<double>> Vecs(W.EigenVecs());

			double Zavg = 0;//<Z> - expectation value of the charge diplacement
			int g = 0;
			if (m == 0) {
				// Find the state with a highest fraction of the original reference state.
				tmp = Vecs[0][0]*Vecs[0][0];
				for (int i = 1; i < Vals.size(); i++) {
					if (Vecs[i][0]*Vecs[i][0] > tmp) {
						tmp = Vecs[i][0]*Vecs[i][0];
						g = i;
					}
				}
				// Store this vector for the reference.
				PrevTVec.push_back(Vecs[g]);
			}
			else {
				// Compare the eigenvectors at T(m) to the reference one at T(m-1),
				// find the one with the largest scalar product and assign "g" to it.
				// If T(m) - T(m-1) is small the correlation should be about 1.
				double MaxCorr = 0;
				for (int i = 0; i < Vals.size(); i++) {
					MaxCorr += Vecs[0][i]*PrevTVec[N][i];
				}
				for (int i = 1; i < Vals.size(); i++) {
					tmp = 0;
					for (int j = 0; j < Vals.size(); j++) {
						tmp += Vecs[i][j]*PrevTVec[N][j];
					}
					if (fabs(tmp) > fabs(MaxCorr)) {
						MaxCorr = tmp;
						g = i;
					}
				}
				if (MaxCorr < 0) {
					for (int j = 0; j < Vals.size(); j++) Vecs[g][j] = -Vecs[g][j];
				}
				// Update PrevTVec.
				PrevTVec[N] = Vecs[g];
			}

			int k = 0;
			// Calculate Zavg for that state
			for (int i = 1; i < Vals.size(); i++) {
				k = MixMe[N].excited[i-1][0];
				e = MixMe[N].excited[i-1][1];
				/*wght = 0.25*MixMe[N].reference[k]*(4*OrbList[e]->L() + 2 - OrbList[e]->occupancy())
					   /(2*OrbList[k]->L() + 1)/(2*OrbList[e]->L() + 1);*/
				wght = 2*OrbList[g]->occupancy()*(1 - 1.*OrbList[e]->occupancy()/(4*OrbList[e]->L() + 2))/3/(2*OrbList[g]->L() + 1);
				wght = sqrt(wght);
				Zavg += 2*Vecs[g][i]*MixMe[N].Dipoles[i-1]*wght;
			}
			Result[m] += P[MixMe[N].index][m]*Zavg;
		} 
	}

	return Result;
}

vector<double> RateEquationSolver::Secular(vector<RadialWF> & Virtual, double Dist, double Einit = 0)
{
	// Alternative perturbation theory. Designed for Ilme project.
	// Two identical atoms are separated by Dist (in a.u.).
	// As they ionize, they create time dependent electric field
	// that distorts their electronic configurations.
	// This functions calculates the field at the nucleus of atom 1
	// created by the atom 2. The field is assumed to be homogenous, and a 
	// small version of CI is calculated. 
	// The function returns a position of charge centroid relative to the nucleus as a function of time.
	vector<double> Result(T.size(), 0);
	vector<double> Efield(T.size(), 0);
	vector<vector<double>> PrevTVec(0);// Reference eigenvector to keep track of ground state configuration.
	vector<RadialWF*> OrbList(0);
	for (auto& Orb: orbitals) OrbList.push_back(&Orb);
	for (auto& Vrt: Virtual) OrbList.push_back(&Vrt);

	// Calculate time dependent Efield at nucleus position Dist. For now it is fixed thorugh the pulse.
	double tmp = 0;
	for (int m = 0; m < T.size(); m++) {
		tmp = 0;
		for (int Z = 0; Z < charge.size(); Z++) {
			tmp -= Z*charge[Z][m];
		}
		Efield[m] = tmp/Dist/Dist - Einit;
	}
	// For each configuration i run the configuration mixing for given Efield[m].
	// Multiply by P[i][m] and store.
	double wght = 0; // Weight factor accounts for fractional occupancies.
	int g = 0;// Ground state orbital to be mixed.
	int e = 0;// Excited state orbital to be mixed.
	for (int m = 0; m < T.size(); m++) {
		for (int N = 0; N < MixMe.size(); N++) {
			// Locate semi-degenerate states.
			int L = MixMe[N].excited.size()+1;

			vector<vector<double>> Fock(L, vector<double>(L, 0));
			vector<vector<double>> Overlap(L, vector<double>(L, 0));

			Fock[0][0] = MixMe[N].refEnergy;
			Overlap[0][0] = 1;
			for (int i = 1; i < L; i++) {
				g = MixMe[N].excited[i-1][0];
				e = MixMe[N].excited[i-1][1];
				wght = sqrt(0.5*MixMe[N].reference[g]*(4*OrbList[e]->L() + 2 - OrbList[e]->occupancy())
					   /(2*OrbList[e]->L() + 1));//(2*OrbList[g]->L() + 1)
				wght = sqrt(wght);
				Overlap[i][i] = 1;
				Fock[0][i] = -Efield[m]*wght*MixMe[N].Dipoles[i-1];
				Fock[i][0] = Fock[0][i];
				Fock[i][i] = MixMe[N].extEnergy[i-1];
			}
			EigenSolver W;
			W.SolveGenEig(Fock, Overlap);

			vector<double> Vals(W.EigenVals());
			vector<vector<double>> Vecs(W.EigenVecs());

			double Zavg = 0;//<Z> - expectation value of the charge diplacement
			int g = 0;
			if (m == 0) {
				// Find the state with a highest fraction of the original reference state.
				tmp = Vecs[0][0]*Vecs[0][0];
				for (int i = 1; i < Vals.size(); i++) {
					if (Vecs[i][0]*Vecs[i][0] > tmp) {
						tmp = Vecs[i][0]*Vecs[i][0];
						g = i;
					}
				}
				// Store this vector for the reference.
				PrevTVec.push_back(Vecs[g]);
			}
			else {
				// Compare the eigenvectors at T(m) to the reference one at T(m-1),
				// find the one with the largest scalar product and assign "g" to it.
				// If T(m) - T(m-1) is small the correlation should be about 1.
				double MaxCorr = 0;
				for (int i = 0; i < Vals.size(); i++) {
					MaxCorr += Vecs[0][i]*PrevTVec[N][i];
				}
				for (int i = 1; i < Vals.size(); i++) {
					tmp = 0;
					for (int j = 0; j < Vals.size(); j++) {
						tmp += Vecs[i][j]*PrevTVec[N][j];
					}
					if (fabs(tmp) > fabs(MaxCorr)) {
						MaxCorr = tmp;
						g = i;
					}
				}
				if (MaxCorr < 0) {
					for (int j = 0; j < Vals.size(); j++) Vecs[g][j] = -Vecs[g][j];
				}
				// Update PrevTVec.
				PrevTVec[N] = Vecs[g];
			}

			int k = 0;
			// Calculate Zavg for that state
			for (int i = 1; i < Vals.size(); i++) {
				k = MixMe[N].excited[i-1][0];
				e = MixMe[N].excited[i-1][1];
				wght = sqrt(0.25*MixMe[N].reference[k]*(4*OrbList[e]->L() + 2 - OrbList[e]->occupancy())
					   /(2*OrbList[k]->L() + 1)/(2*OrbList[e]->L() + 1));
				Zavg += 2*Vecs[g][i]*MixMe[N].Dipoles[i-1]*wght;
			}
			Result[m] += P[MixMe[N].index][m]*Zavg;
		}
	}

	return Result;
}


int RateEquationSolver::SetupAndSolve(ofstream & runlog, int out_T_size)
{
	// initial conditions for rate equation
	// first index represents the configuration
	// second index represents the time. The routine will expand it itself unless
	// you provide first adams_n points yourself.
	double fluence = input.Fluence() / Constant::Fluence_in_au;
	double Sigma = input.Width()/Constant::fs_in_au/ (2*sqrt(2*log(2.)));
	int T_size = input.TimePts();
	vector<double> InitCond(dimension, 0);
	InitCond[0] = 1;
	P.clear();

	Store.nAtoms = 0.00744;

	vector<double> Intensity;
	double scaling_T = 1;

	int converged = 1;
	Plasma Mxwll(T.size());
	while (converged != 0)
	{
		dT.clear();
		dT = generate_dT(T_size);
		T.clear();
		T = generate_T(dT);
		scaling_T = 4 * input.Width()/Constant::fs_in_au / T.back();
		for (int i = 0; i < T.size(); i++) {
			T[i] *= scaling_T;
			dT[i] *= scaling_T;
		}
		Intensity = generate_I(T, fluence, Sigma);
		//SmoothOrigin(T, Intensity);
		//Mxwll.resize(T.size());
		IntegrateRateEquation Calc(dT, T, Store, InitCond, Intensity);
		converged = Calc.Solve( 0, 1, out_T_size);//(Mxwll,
		if (converged == 0)
		{
			cout << "Final number of time steps: " << T_size << endl;
			P = Calc.GetP();
			T.clear();
			T = Calc.GetT();
			dT.clear();
			dT = generate_dT(T.size());
		}
		else
		{
			cout << "Diverged at step: " << converged << " of " << T_size << endl;
			T_size *= 2;
		}
	}


	Intensity.clear();
	Intensity = generate_I(T, fluence, Sigma);
	SmoothOrigin(T, Intensity);

	int tmp = 0;
	for (int i = 0; i < Index.back().size(); i++) tmp += Index.back()[i] - Index.begin()->at(i);
	charge.clear();
	charge.resize(tmp + 1);
	for (auto& ch: charge) ch = vector<double>(T.size(), 0);
	for (int i = 0; i < P.size(); i++)
	{
		tmp = Charge(i);
		for (int m = 0; m < P[i].size(); m++) {
			charge[tmp][m] += P[i][m];
		}
	}
	
	double t_tmp = 0;
	string ChargeName = "./output/Charge_" + input.Name() + ".txt";
	string IntensityName = "./output/Intensity_" + input.Name() + ".txt";
	ofstream charge_out(ChargeName);
	ofstream intensity_out(IntensityName);

	double chrg_tmp = 0;
	double I_max = *max_element(begin(Intensity), end(Intensity));
	for (int m = 0; m < T.size(); m++) {
		T[m] = (T[m]-0.5*T.back())*Constant::fs_in_au;
		dT[m] *= Constant::fs_in_au;
		intensity_out << Intensity[m] / I_max << " " << T[m] << endl;

		for (int i = 0; i < charge.size(); i++) {
			chrg_tmp = charge[i][m];
			if (chrg_tmp <= 0.00000001) charge_out << 0 << " ";
			else charge_out << chrg_tmp << " ";			
		}
		
		charge_out << T[m] << endl;
	}

	intensity_out.close();
	charge_out.close();
	
  /*
	FILE * fl;
	fl = fopen("./output/Plasma.txt", "w");
	int cnt = Mxwll.N.size()/T.size();
	for (int m = 0; m < T.size(); m++) {
		fprintf(fl, "%3.5e %3.5e %3.5f %3.5f\n", Mxwll.E[m*cnt], Mxwll.N[m*cnt], Mxwll.Np[m*cnt], T[m]);
	}
	fclose(fl);
  */
	return 0;
}


int RateEquationSolver::SetupAndSolve(MolInp & Input, ofstream & runlog, int out_T_size)
{
	// initial conditions for rate equation
	// first index represents the configuration
	// second index represents the time. The routine will expand it itself unless
	// you provide first adams_n points yourself.
	double fluence = Input.Fluence()/ Constant::Fluence_in_au;
	double Sigma = Input.Width()/Constant::fs_in_au / (2*sqrt(2*log(2.)));
	int T_size = Input.ini_T_size();

	P.clear();

	vector<double> Intensity;
	double scaling_T = 1;

	int converged = 1;
	Plasma Mxwll(T.size());
	while (converged != 0)
	{
		dT.clear();
		dT = generate_dT(T_size);
		T.clear();
		T = generate_T(dT);
		scaling_T = 2.666*input.Width()/Constant::fs_in_au / T.back();//4
		for (int i = 0; i < T.size(); i++) {
      //T[i] = T[i]-0.5*T.back();
			T[i] *= scaling_T;
			dT[i] *= scaling_T;
		}
		Intensity = generate_I(T, fluence, Sigma);
    
    // Extend mesh to 100 fs (extra 20fs added as pulse starts at -20fs).
    extend_I(Intensity,120/Constant::fs_in_au, dT.back());
    
		//SmoothOrigin(T, Intensity);
 		Mxwll.resize(T.size());
		IntegrateRateEquation Calc(dT, T, Input.Store, Mxwll, Intensity);

    T_size = T.size();

		converged = Calc.Solve(Mxwll, Input.Store, T_size);

		if (converged == 0)
		{
			cout << "Final number of time steps: " << T_size << endl;
			P = Calc.GetP();
		//	T.clear();
		//	T = Calc.GetT();
		//	dT.clear();
		//	dT = generate_dT(T.size());
		}
		else
		{
			cout << "Diverged at step: " << converged << " of " << T_size << endl;
			T_size *= 2;
		}
	}
	
	// TODO: adjust charge calculations as they are adapted for a single atom!!!

//	Intensity.clear();
//	Intensity = generate_I(T, fluence, Sigma);

	double t_tmp = 0;

	string IntensityName = "./output/Intensity_" + Input.name + ".txt";

  // Output intensity profile.
	//ofstream intensity_out(IntensityName);
	double I_max = *max_element(begin(Intensity), end(Intensity));
  int MidT = 0;
  while (Intensity[MidT] != I_max && MidT < T.size()) MidT++;
  double MidTval = T[MidT];
	for (int m = 0; m < T.size(); m++) {
		T[m] = (T[m]-MidTval)*Constant::fs_in_au;
		//dT[m] *= Constant::fs_in_au;
		//intensity_out << Intensity[m] / I_max << " " << T[m] << endl;
	}
	//intensity_out.close();

	int shift = 0;
	vector<int> P_to_charge(0);

  // Aggregate and output charges, plasma parameters, and other parameters into an output.
  
  // Charges.
  vector<vector<double>> AllAtomCharge(Input.Atomic.size(), vector<double>(T.size(), 0));
  // Ionic radius.
  vector<vector<double>> R_ion(Input.Atomic.size(), vector<double>(T.size(), 0) );
  // Form-factor
  vector<vector<ffactor>> FF(Input.Atomic.size(), vector<ffactor>(T.size()));
  
	for (int a = 0; a < Input.Atomic.size(); a++) {
		
		// Occupancies associated with the atom "a".
		vector<vector<double*>> map_p(Input.Store[a].num_conf);
		for (int i = 0; i < Input.Store[a].num_conf; i++) {
			map_p[a].push_back(P[i + shift].data());
		}

    for (int m = 0; m < T.size(); m++) {
      for (int i = 0; i < Input.Store[a].num_conf; i++) {
        R_ion[a][m] += Input.AuxStore[a].r_ion[i] * *(map_p[a][i] + m);
      }
    }

    for (int m = 0; m < T.size(); m++) {
      FF[a][m].val.clear();
      FF[a][m].val.resize(Input.AuxStore[a].form_fact_ion[0].size(), 0);
      FF[a][m].index = a;

      for (int i = 0; i < Input.Store[a].num_conf-1; i++) {
        for (int q = 0; q < Input.AuxStore[a].form_fact_ion[0].size(); q++) {
          FF[a][m].val[q] += Input.AuxStore[a].form_fact_ion[i][q] * *(map_p[a][i] + m);
        }
      }
    }

		//string ChargeName = "./output/Charge_" + Input.Store[a].name + ".txt";
		
		//ofstream charge_out(ChargeName);

		// Work out charges of all configurations.
		double chrg_tmp = 0;
		int tmp = Input.Index[a][0].size();
		P_to_charge.clear();
		P_to_charge.resize(Input.Store[a].num_conf, 0);
		for (int i = 0; i < Input.Index[a].size(); i++) {
			for (int j = 0; j < tmp; j++) P_to_charge[i] += Input.Index[a][i][j];
			chrg_tmp = 0;
		}

		charge.clear();
		charge.resize(*max_element(begin(P_to_charge), end(P_to_charge)) + 1);
		for (auto& ch: charge) ch = vector<double>(T.size(), 0);
		for (int i = 0; i < map_p[a].size(); i++)
		{
			tmp = P_to_charge[i];
			for (int m = 0; m < P[i].size(); m++) charge[tmp][m] += *(map_p[a][i] + m);
		}		

    // Data for Harry.
    for (int m = 0; m < T.size(); m++) {
			for (int i = 1; i < charge.size(); i++) {
        AllAtomCharge[a][m] += i*charge[i][m];
      }
    }
    /*
		for (int m = 0; m < T.size(); m++) {
			for (int i = 0; i < charge.size(); i++) {
				chrg_tmp = charge[i][m];
				if (chrg_tmp <= 0.00000001) charge_out << 0 << " ";
				else charge_out << chrg_tmp << " ";			
			}
			
			charge_out << T[m] << endl;
		}
    */

		//charge_out.close();
		shift += Input.Store[a].num_conf;
	}

  /*
	FILE * fl;
	fl = fopen("./output/Plasma.txt", "w");
  int size = 500;
	int cnt = Mxwll.N.size()/size;//T.size();
	for (int m = 0; m < T.size(); m++) {
		fprintf(fl, "%3.5e %3.5e %3.5e %3.5f\n", Mxwll.E[m], Mxwll.N[m], Mxwll.Np[m], T[m]);//Mxwll.Np[m*cnt]
	}

	fclose(fl);*/
  {
    string MD_Data = "./output/MD_Data.txt";
    ofstream OutFile(MD_Data);
    OutFile << T.size() << endl;
    OutFile << "Time ";
    for (int a = 0; a < Input.Atomic.size(); a++) OutFile << Input.Atomic[a].Nuclear_Z() << " ";
    OutFile << "N(elec) E(elec)";
    for (int m = 0; m < T.size(); m++) {
      OutFile << endl << T[m] << " ";
      for (int a = 0; a < AllAtomCharge.size(); a++) {
        OutFile << AllAtomCharge[a][m] << " ";
      }
      if (m != 0) OutFile << Mxwll.N[m] << " " << Mxwll.E[m];
      else OutFile << 0 << " " << 0;
    }
  }

  {
    string FF_Data = "./output/FF_Data.txt";
    ofstream OutFile(FF_Data);
    OutFile << "Time-dependent form factors.\n";
    OutFile << "Atomic FF are output by blocks: row index is 'time', column index is 'q' as in 'Q_mesh'.\n";
    OutFile << "line 1 : atomic proton numbers. Matches number of blocks.\n";
    OutFile << "line 2 : Q_mesh. Matches umber of columns.\n";

    for (int a = 0; a < Input.Atomic.size(); a++) OutFile << Input.Atomic[a].Nuclear_Z() << " ";
    OutFile << "\n";
    FormFactor W(lattice, charge[0], 0); // Dummy initialization to retrive Q_mesh.
    vector<double> Q_mesh = W.getAllQ();
    for (auto & q: Q_mesh) OutFile << q << " ";

    for (int a = 0; a < Input.AuxStore.size(); a++) {
      OutFile << "\n";
      for (int m = 0; m < T.size(); m++) {
        double * pVal = FF[a][m].val.data();
        OutFile << "\n"; 
        for (int q = 0; q < Q_mesh.size(); q++) {
          OutFile << *(pVal + q) << " ";
        }  
      }
    }
  }


	return 0;
}


string RateEquationSolver::InterpretIndex(int i)
{
	// LaTeX output format for direct insertion into table.
	string Result;
	if (!Index.empty()) {
		if (Index[i].size() == orbitals.size())	{
			ostringstream tmpstr;
			tmpstr << '$';
			for (int j = 0; j < orbitals.size(); j++) {
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
				tmpstr << '^' << orbitals[j].occupancy() - Index[i][j];
			}
			tmpstr << '$';
			Result = tmpstr.str();
		}
	}
	return Result;
}

int RateEquationSolver::Charge(int Iconf)
{
	int Result = 0;
	for (int j = 0; j < Index[Iconf].size(); j++) Result += Index[Iconf][j] - Index[0][j];
	return Result;
}

bool RateEquationSolver::SetupIndex(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
{
	if (orbitals.size() != Final_occ.size())
	{
		runlog << "Final occupancies should be provided for all orbitals." << endl;
		return false;
	}
	if (Max_occ.size() != orbitals.size()) {
		Max_occ.resize(orbitals.size());
		for (int i = 0; i < orbitals.size(); i++) {
			Max_occ[i] = 4*orbitals[i].L() + 2;
		}
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

	for (int i = 0; i < dimension; i++)
	{
		int tmp = i;
		for (int j = 0; j < orbitals_size; j++)
		{
			Index[i][j] = tmp / hole_posit[j];
			if (Index[i][j] > Max_occ[j]) { Index[i][j] = Max_occ[j]; }
			tmp -= Index[i][j] * hole_posit[j];
		}
	}
	return true;
}

RateEquationSolver::~RateEquationSolver()
{
}

inline bool exists_test(const std::string& name)
{
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}

vector<double> RateEquationSolver::generate_dT(int num_elem)//default time interval
{
	vector<double> Result(num_elem, 0);
	double tmp = 1;
	for (int i = 0; i < num_elem; i++) {
	tmp = fabs(1.*i / (num_elem-1) - 0.5) + 0.01;
		Result[i] = tmp;
	}
	return Result;
}

vector<double> RateEquationSolver::generate_T(vector<double>& dT)//default time
{
	vector<double> Result(dT.size(), 0);
	vector<double> Bashforth_4{ 55. / 24., -59. / 24., 37. / 24., -9. / 24. }; //Adams�Bashforth method
	for (int i = 1; i < Bashforth_4.size(); i++)//initial few points
	{
		Result[i] = Result[i - 1] + dT[i - 1];
	}
	for (int i = Bashforth_4.size(); i < Result.size(); i++)//subsequent points
	{
		Result[i] = Result[i - 1];
		for (int j = 0; j < Bashforth_4.size(); j++)
		{
			Result[i] += Bashforth_4[j] * dT[i - j - 1];
		}
	}

	return Result;
}

vector<double> RateEquationSolver::generate_I(vector<double>& Time, double Fluence, double Sigma)//intensity of the Gaussian X-ray pulse
{

	vector<double> Result(Time.size(), 0);
	
	double midpoint = (Time.back() + T[0]) / 2;
	double denom = 2*Sigma*Sigma;
	double norm = 1./sqrt(denom*Constant::Pi);
	//include the window function to make it exactly 0 at the beginning and smoothly increase toward Gaussian
  
  
  int smooth = T.size()/10;
	for (int i = 0; i < Time.size(); i++)
	{
    Result[i] = Fluence * norm * exp(-(Time[i] - midpoint)*(Time[i] - midpoint) / denom);
    if (i < smooth) Result[i] *= (T[i] / T[smooth])*(T[i] / T[smooth])*(3 - 2 * (T[i] / T[smooth]));
  }
  
  /*
  for (int i = 0; i < Time.size(); i++)
	{
    Result[i] = Fluence / T.back();
	}
  */
  
  

	return Result;
}

int RateEquationSolver::extend_I(vector<double>& Intensity, double new_max_T, double step_T)
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

vector<double> RateEquationSolver::generate_G()
{
	// Intensity profile normalized to 1.
	// Time is assumbed to be in FEM
	double Sigma = input.Width()/Constant::fs_in_au/ (2*sqrt(2*log(2.)));

	return generate_I(T, 1, Sigma);
}

void RateEquationSolver::GenerateRateKeys(vector<Rate> & ToSort)
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

int RateEquationSolver::mapOccInd(vector<RadialWF> & Orbitals)
{
	int Result = 0;
	for (int j = 0; j < hole_posit.size(); j++)	{
		Result += (orbitals[j].occupancy() - Orbitals[j].occupancy())*hole_posit[j];
	}

	return Result;
}