#include "RateData.h"

using namespace RateData::Atom = Atom;


void Atom::save_csv(const std::filesystem::path root,  char which)
{


}


void Atom::load_csv(const std::filesystem::path root,  char which)
{
    if (which & PHOTO != 0){

    }
    if (which & FLUOR != 0){

    }
    if (which & AUGER != 0){

    }
    if (which & EII != 0){

    }
}

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



// Reads a ratefile input and stores the data in PutHere
// Returns true on successful opening
bool ReadRates(const std::string & input, std::vector<Rate>& PutHere) {
    PutHere.clear();

    Rate Tmp;
    std::ifstream infile;
    infile.open(input);
    if (!infile.good()) {
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

// Serialises EII params into a tabular format
// where init, fin are indices (initial and final states),
// occ is the initial state occupancy
// ionB, kin are calculated parameters
// init, fin, occ, ionB, kin
bool ReadEIIParams(const std::string & input, std::vector<EIIdata> & PutHere) {
    PutHere.clear();
    std::ifstream infile;
    infile.open(input);
    if (!infile.good()) {
        return false;
    }

    std::string line;

    EIIdata tmp;
    std::string tmpstr;

    unsigned tmp_init;
    unsigned tmp_fin;
    unsigned tmp_occ;
    double tmp_ionB, tmp_kin;
    
    while (!infile.eof())
    {
        std::string line;
        getline(infile, line);
        if (line[0] == '#') continue; // skip comments
        if (line[0] == '\0') continue; // skip null bytes

        std::stringstream stream(line);
        
        stream >> tmp_init >> tmp_fin >> tmp_occ >> tmp_ionB >> tmp_kin;

        if (PutHere.size() == 0){
            tmp.init = tmp_init;
        }else if { tmp_init != tmp.init){
            PutHere.push_back(tmp);
        }
        

        
    }

    infile.close();
    return true; // Returns true if all went to plan

};

void WriteRates(const std::string& fname, const std::vector<Rate>& rates) {
    FILE * fl = fopen(fname.c_str(), "w");
    fprintf(fl, "# val from to energy(Ha)\n");
    for (auto& R : rates) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
    fclose(fl);
}

void WriteEIIParams(const std::string& fname, const std::vector<EIIdata>& rates) {
    FILE * fl = fopen(fname.c_str(), "w");
    fprintf(fl, "# from to occ ionB(Ha) kin(Ha)\n");
    for (auto&R : rates) {
        for (size_t j=0; j<R.fin.size(); j++){
            fprintf(fl, "%6ld %6ld %3ld %1.8e %1.8e\n", R.init, R.fin[j], R.occ[j], R.ionB[j], R.kin[j]);
        }
    }
    fclose(fl);
}
