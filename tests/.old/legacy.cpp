// Contains random commented-out snippets of code for future reference, should they be needed.

/*
int ComputeRateParam::SolveFrozen(vector<int> Max_occ, vector<int> Final_occ, ofstream & runlog)
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

	bool existPht = !RateData::ReadRates(RateLocation + "Photo.txt", Store.Photo);
	bool existFlr = !RateData::ReadRates(RateLocation + "Fluor.txt", Store.Fluor);
	bool existAug = !RateData::ReadRates(RateLocation + "Auger.txt", Store.Auger);

	if (existPht) printf("Photoionization rates found. Reading...\n");
	if (existFlr) printf("Fluorescence rates found. Reading...\n");
	if (existAug) printf("Auger rates found. Reading...\n");

	// Electron density evaluation.
	density.clear();
	density.resize(dimension - 1);
	for (auto& dens: density) dens.resize(lattice.size(), 0.);
  vector<pair<double, int>> conf_RMS(0);

	if (existPht || existFlr || existAug || true)
	{
		cout << "No rates found. Calculating..." << endl;
		cout << "Total number of configurations: " << dimension << endl;
		Rate Tmp;
		vector<RateData::Rate> LocalPhoto(0);
		vector<RateData::Rate> LocalFluor(0);
		vector<RateData::Rate> LocalAuger(0);

		omp_set_num_threads(input.Num_Threads());
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
				HartreeFock HF(lattice, Orbitals, U, input, runlog);

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

			}

			#pragma omp critical
			{
				Store.Photo.insert(Store.Photo.end(), LocalPhoto.begin(), LocalPhoto.end());
				Store.Fluor.insert(Store.Fluor.end(), LocalFluor.begin(), LocalFluor.end());
				Store.Auger.insert(Store.Auger.end(), LocalAuger.begin(), LocalAuger.end());
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
  config_out.close();

 	return dimension;
}
*/
