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
#include <vector>>
#include <algorithm>
#include "Constant.h"
#include <cmath>
#include <fstream>
#include <sstream>
#include <iostream>

double LogFactorialFraction(double Num, double Denom);

namespace Constant
{
	double Wigner3j(double j1, double j2, double j3, double m1, double m2, double m3)
	{
		double Result = 0.;
		double fracpart, intpart, J = j1 + j2 + j3;

		//some particular cases can be calculated in simpler and faster way.
		if (m1 + m2 + m3 != 0) return 0.;
		if (j1 < 0. || j2 < 0. || j3 < 0.) return 0.;
		fracpart = std::modf(J, &intpart);// J should be integer
		if (fracpart != 0) return 0;
		else if (m1 == 0 && m2 == 0)
		{
			// /j1 j2 j3\ - Edmonds, "Angular momentum in qantum mechanics", p 50
			// \0  0  0 /
			fracpart = std::modf(0.5*J, &intpart);
			if (fracpart != 0) return 0;
			else
			{
				if (j1 > j2 + j3 || j2 > j1 + j3 || j3 > j1 + j2 ||
					j1 < fabs(j2 - j3) || j2 < fabs(j1 - j3) || j3 < fabs(j2 - j1)) return 0;
				std::vector<double> factor(3);
				factor[0] = J / 2. - j1;
				factor[1] = J / 2. - j2;
				factor[2] = J / 2. - j3;
				std::sort(factor.begin(), factor.end());//smallest first, largest last

				for (unsigned int i = 0; i < factor.size(); i++)
				{
					fracpart = std::modf(factor[i], &intpart);// J should be integer
					if (fracpart != 0) return 0;
				}

				Result += LogFactorialFraction(2 * factor[0], factor[2]);//taking the largest factor under the square root in denominator (factor[2]!)*(factor[2]!)
				Result += LogFactorialFraction(2 * factor[1], factor[2]);
				Result += LogFactorialFraction(2 * factor[2], J + 1);
				Result *= 0.5;// end of square root

				Result += LogFactorialFraction(J / 2, factor[1]);
				Result += LogFactorialFraction(1, factor[0]);//smallest factor goes unpaired with numerator

				fracpart = std::modf(0.25*J, &intpart);//fracpart = 0 if J/2 is even
				if (fracpart == 0) { Result = exp(Result); }
				else { Result = -exp(Result); }

				return Result;
			}
		}
		else if ((m1 == 0 || m2 == 0 || m3 == 0) && (m1 == 0.5 || m2 == 0.5 || m3 == 0.5))
		{
			// /ja   jb   J\ - Brink and Satcher, "Angular momentum", p 138
			// \0.5 -0.5  0 /
			double ja, jb;
			int columns_permutation = 1;//account for sign flip
			if (j1 > j2 + j3 || j2 > j1 + j3 || j3 > j1 + j2 ||
				j1 < fabs(j2 - j3) || j2 < fabs(j1 - j3) || j3 < fabs(j2 - j1)) return 0;

			if (m1 == 0)
			{
				J = j1;
				if (m2 == 0.5)
				{
					ja = j2;
					jb = j3;
				}
				else
				{
					ja = j3;
					jb = j2;
					columns_permutation = -1;
				}
			}
			else if (m1 == 0.5)
			{
				ja = j1;
				if (m2 == 0)
				{
					J = j2;
					jb = j3;
					columns_permutation = -1;
				}
				else
				{
					J = j3;
					jb = j2;
				}
			}
			else
			{
				jb = j1;
				if (m2 == 0)
				{
					J = j2;
					ja = j3;
				}
				else
				{
					J = j3;
					ja = j2;
					columns_permutation = -1;
				}
			}

			double K;
			fracpart = std::modf(0.5*(ja + jb + J), &intpart);
			if (fracpart == 0) { K = J; }
			else  { K = J + 1; }

			std::vector<double> factor_qsrt(3);//sort factorials under square root
			factor_qsrt[0] = ja + jb - J;
			factor_qsrt[1] = ja + J - jb;
			factor_qsrt[2] = jb + J - ja;
			std::sort(factor_qsrt.begin(), factor_qsrt.end());

			std::vector<double> factor(3);//sort factorials under square root
			factor[0] = (ja + jb - K) / 2;
			factor[1] = (ja + K - jb - 1) / 2;
			factor[2] = (jb + J - ja) / 2;
			std::sort(factor.begin(), factor.end());

			Result += LogFactorialFraction(factor_qsrt[0], factor[2]);
			Result += LogFactorialFraction(factor_qsrt[1], factor[2]);
			Result += LogFactorialFraction(factor_qsrt[2], (J + ja + jb + 1));
			Result *= 0.5;

			Result += LogFactorialFraction(1, factor[0]);
			Result += LogFactorialFraction(0.5*(K + ja + jb), factor[1]);

			fracpart = std::modf(0.25*(K + ja + jb) - 0.5, &intpart);
			if (fracpart != 0) columns_permutation *= -1;
			Result = columns_permutation*exp(Result)*2. / sqrt((2 * ja + 1)*(2 * jb + 1));

			return Result;
		}
		else
		{
			std::vector<double> Regge(9);
			std::vector<double> ReggeInt(9);
			Regge[0] = j2 + j3 - j1;//R(11)
			Regge[1] = j1 - j2 + j3;//R(12)
			Regge[2] = j1 + j2 - j3;//R(13)
			Regge[3] = j1 + m1;//R(21)
			Regge[4] = j2 + m2;//R(22)
			Regge[5] = j3 + m3;//R(23)
			Regge[6] = j1 - m1;//R(31)
			Regge[7] = j2 - m2;//R(32)
			Regge[8] = j3 - m3;//R(33)

			// if any Regge symbol is <0 or not integer, 3j symbol = 0
			for (unsigned int i = 0; i < Regge.size(); i++)
			{
				if (Regge[i] < 0.) return 0.;
				fracpart = std::modf(Regge[i], &intpart);
				if (Regge[i] != intpart) return 0.;
				ReggeInt[i] = (int)intpart;
			}

			//using 8.3.29 of Varshalovich. Outer part is everything under the square root==================NEEDS TO BE FINISHED=======================

			return 0.;
		}

	}
}

//this function evaluates log(Num!/Denom!). It is used to calculate factorial summations in Wigner3j funstion. Allows to avoid getting huge factorial products
double LogFactorialFraction(double Num, double Denom)
{
	double Result = 0.;
	if (Num > Denom)
	{
		for (int i = Denom + 1; i <= Num; i++)
		{
			Result += log((double)i);
		}
	}
	else if (Num < Denom)
	{
		for (int i = Num + 1; i <= Denom; i++)
		{
			Result -= log((double)i);
		}
	}
	return Result;
}

namespace RateData{
	std::vector<InverseEIIdata> inverse(const std::vector<EIIdata>& eiiVec)
	{
		std::vector<InverseEIIdata> tbrVec(eiiVec.size()+1);
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
	void read_vector(const std::string& s, std::vector<T>&v) {
		auto iss = std::istringstream(s);

		std::string str;
		T tmp;
		while (iss >> str) {
			auto ss = std::stringstream(str);
			ss >> tmp;
			v.push_back(tmp);
		}
	}

	std::string find_bracket_contents(std::string &src, char left, char right) {
		size_t first_idx = src.find(left);
		size_t last_idx = src.rfind(right);
		return src.substr(first_idx+1, last_idx-first_idx-1);
	}

	// Reads a ratefile input and stores the data in PutHere
	// Returns true on successful opening
	bool ReadRates(const std::string & input, std::vector<Rate>& PutHere) {
		PutHere.clear();

		Rate Tmp;
		std::ifstream infile;
		infile.open(input);
		if (!infile.is_open()) {
			return false;
		}

		while (!infile.eof())
		{
			std::string line;
			getline(infile, line);
			if (line[0] == '#') continue; // skip comments
			if (line[0] == '\0') continue; // skip null bytes

			std::stringstream stream(line);
			stream >> Tmp.val >> Tmp.from >> Tmp.to >> Tmp.energy;

			PutHere.push_back(Tmp);
		}

		infile.close();
		return true; // Returns true if all went well

	}

	// reads a "JSON-style" file and stores the data in an EIIdata structure.
	// Returns true on success
	// Possibly the worst parser ever written.
	bool ReadEIIParams(const std::string & input, std::vector<EIIdata> & PutHere) {
		PutHere.clear();
		std::ifstream infile;
		infile.open(input);
		if (!infile.is_open()) {
			return false;
		}

		std::string line;

		EIIdata tmp;
		std::string tmpstr;
		bool in_top_brace = false;
		bool in_eii_record = false;

		int i=0; // configuration counter
		while (!infile.eof())
		{

			getline(infile, line);
			#ifdef DEBUG_VERBOSE
			std::cout<<"[ DEBUG ] [ ReadEII ] "<<line<<"\n";
			#endif
			if (line[0] == '#') continue; // skip comments
			if (in_top_brace) {
				if (in_eii_record) {
					if (line[0] == '}') {
						in_eii_record = false;
						PutHere.push_back(tmp);
						continue;
					}
					std::string pref = "  \"fin\":";
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
					std::string prefix="\"configuration ";
					if(line.compare(0, prefix.size(), prefix) == 0) {
						// New entry in the array
						in_eii_record = true;
						tmpstr = line.substr(prefix.size());
						int j = stoi(tmpstr.substr(0,tmpstr.find('"')));
						tmp.init = j;
						tmp.resize(0);
						if(i != j) std::cerr<<"[ ReadEII ] Unexpected index: got "<<j<<" expected "<<i<<"\n";
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

	void WriteRates(const std::string& fname, const std::vector<Rate>& rates) {
		FILE * fl = safe_fopen(fname.c_str(), "w");
		fprintf(fl, "# val from to energy(Ha)\n");
		for (auto& R : rates) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
		fclose(fl);
	}

	void WriteEIIParams(const std::string& fname, const std::vector<EIIdata>& rates) {
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

