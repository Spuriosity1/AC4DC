#include "RateData.h"
#include <fstream>
#include <sstream>
#include <iostream>

using namespace std;

namespace RateData{
	vector<InverseEIIdata> inverse(const vector<EIIdata>& eiiVec)
	{
		vector<InverseEIIdata> tbrVec(eiiVec.size()+1);
		for (size_t i = 0; i < tbrVec.size(); i++) {
			tbrVec[i].init = i;
		}
		for (auto& eii : eiiVec) {
			for (size_t j=0; j<eii.fin.size(); j++) {
				InverseEIIdata& inv = tbrVec.at(eii.fin[j]);
				inv.fin.push_back(eii.init);
				inv.occ.push_back(eii.occ[j]);
				inv.ionB.push_back(eii.ionB[j]);
				inv.kin.push_back(eii.kin[j]);
			}
		}
		return tbrVec;
	}

	FILE* safe_fopen(const char *filename, const char *mode)
	{
		FILE* fp = fopen(filename, mode);
		if (fp == NULL) {
			std::cerr << "Could not open file '" << filename << "'" << std::endl;
		}
		return fp;
	}

	template<typename T>
	void read_vector(const string& s, vector<T>&v) {
		auto iss = istringstream(s);

		string str;
		T tmp;
		while (iss >> str) {
			auto ss = stringstream(str);
			ss >> tmp;
			v.push_back(tmp);
		}
	}

	string find_bracket_contents(string &src, char left, char right) {
		size_t first_idx = src.find(left);
		size_t last_idx = src.rfind(right);
		return src.substr(first_idx+1, last_idx-first_idx-1);
	}

	// Reads a ratefile input and stores the data in PutHere
	// Returns true on successful opening
	bool ReadRates(const string & input, vector<Rate>& PutHere) {
		PutHere.clear();

		Rate Tmp;
		std::ifstream infile;
		infile.open(input);
		if (!infile.is_open()) {
			return false;
		}

		while (!infile.eof())
		{
			string line;
			getline(infile, line);
			if (line[0] == '#') continue; // skip comments
			if (line[0] == '\0') continue; // skip null bytes

			stringstream stream(line);
			stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

			PutHere.push_back(Tmp);
		}

		infile.close();
		return true; // Returns true if all went well

	}

	// reads a "JSON-style" file and stores the data in an EIIdata structure.
	// Returns true on success
	// Possibly the worst parser ever written.
	bool ReadEIIParams(const string & input, vector<EIIdata> & PutHere) {
		PutHere.clear();
		std::ifstream infile;
		infile.open(input);
		if (!infile.is_open()) {
			return false;
		}

		string line;

		EIIdata tmp;
		string tmpstr;
		bool in_top_brace = false;
		bool in_eii_record = false;

		int i=0; // configuration counter
		while (!infile.eof())
		{

			getline(infile, line);
			#ifdef DEBUG_VERBOSE
			cout<<"[ DEBUG ] [ ReadEII ] "<<line<<endl;
			#endif
			if (line[0] == '#') continue; // skip comments
			if (in_top_brace) {
				if (in_eii_record) {
					if (line[0] == '}') {
						in_eii_record = false;
						PutHere.push_back(tmp);
						continue;
					}
					string pref = "  \"fin\":";
					if(line.compare(0, pref.size(), pref)==0)
					{
						read_vector<int>(find_bracket_contents(line, '[', ']'), tmp.fin);
						continue;
					}
					pref = "  \"occ\":";
					if(line.compare(0, pref.size(), pref)==0)
					{
						read_vector<int>(find_bracket_contents(line, '[', ']'), tmp.occ);
						continue;
					}
					pref = "  \"ionB\":";
					if(line.compare(0, pref.size(), pref)==0)
					{
						read_vector<float>(find_bracket_contents(line, '[', ']'), tmp.ionB);
						continue;
					}
					pref="  \"kin\":";
					if(line.compare(0, pref.size(), pref)==0)
					{
						read_vector<float>(find_bracket_contents(line, '[', ']'), tmp.kin);
						continue;
					}
				} else {
					string prefix="\"configuration ";
					if(line.compare(0, prefix.size(), prefix) == 0) {
						// New entry in the array
						in_eii_record = true;
						tmpstr = line.substr(prefix.size());
						int j = stoi(tmpstr.substr(0,tmpstr.find('"')));
						tmp.init = j;
						tmp.resize(0);
						if(i != j) cerr<<"[ ReadEII ] Unexpected index: got "<<j<<" expected "<<i<<endl;
						i++;
						continue;
					}
				}
			} else {
				if (line[0] == '{') {
					in_top_brace = true;
				}	
				continue;
			} 

			


		}

		infile.close();
		return true; // Returns true if all went to plan

	};

	void WriteRates(const string& fname, const vector<Rate>& rates) {
		FILE * fl = safe_fopen(fname.c_str(), "w");
		fprintf(fl, "# val from to energy(Ha)\n");
		for (auto& R : rates) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
		fclose(fl);
	}

	void WriteEIIParams(const string& fname, const vector<EIIdata>& rates) {
		FILE * fl = safe_fopen(fname.c_str(), "w");
		fprintf(fl, "{\n");
		bool first_entry = true;
		for (auto&R : rates) {
			if (!first_entry) {
				fprintf(fl, ",\n");
			}
			first_entry = false;
			fprintf(fl, "\"configuration %d\": {\n", R.init);
			// This will iterate in order of the inits anyway, so no need to explicitly store them
			// Use a JSON style to deal with the multidimensional data
			fprintf(fl, "  \"fin\": [");
			bool first_iter = true;
			for (auto& f : R.fin) {
				fprintf(fl, first_iter ? "%d" : ", %d", f);
				first_iter = false;
			}
			fprintf(fl, "],\n");
			fprintf(fl, "  \"occ\": [");
			first_iter=true;
			for (auto& f : R.occ) {
				fprintf(fl, first_iter ? "%d" : ", %d", f);
				first_iter = false;
			}
			fprintf(fl, "],\n");
			fprintf(fl, "  \"ionB\": [");
			first_iter=true;
			for (auto& f : R.ionB) {
				fprintf(fl, first_iter ? "%f" : ", %f", f);
				first_iter = false;
			}
			fprintf(fl, "],\n");
			fprintf(fl, "  \"kin\": [");
			first_iter=true;
			for (auto& f : R.kin) {
				fprintf(fl, first_iter ? "%f" : ", %f", f);
				first_iter = false;
			}
			fprintf(fl, "]\n");
			fprintf(fl, "}");
		}
		fprintf(fl, "\n}\n");
		fclose(fl);
	}

}; // end of namespace: RateData
