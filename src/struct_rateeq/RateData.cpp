#include "RateData.hpp"

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

// template<typename T>
// void read_vector(const std::string& s, std::vector<T>&v) {
//     auto iss = std::istringstream(s);

//     std::string str;
//     T tmp;
//     while (iss >> str) {
//         auto ss = std::stringstream(str);
//         ss >> tmp;
//         v.push_back(tmp);
//     }
// }
/*

// DEPRECEATED
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

// DEPRECEATED
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

// DEPRECEATED
void WriteRates(const std::string& fname, const std::vector<Rate>& rates) {
    FILE * fl = fopen(fname.c_str(), "w");
    fprintf(fl, "# val from to energy(Ha)\n");
    for (auto& R : rates) fprintf(fl, "%1.8e %6ld %6ld %1.8e\n", R.val, R.from, R.to, R.energy);
    fclose(fl);
}

// DEPRECEATED
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
*/


// IO version 2
void RateData::save_csv(const std::filesystem::path& root, const std::vector<Rate>& rates){
    std::string fname(root + ".csv");
    FILE * fl = fopen(fname.c_str(), "w");
    fprintf(fl, "# from to rate(1/Ha) energy(Ha)\n");
    for (auto& R : rates) fprintf(fl, "%6ld %6ld %1.8e %1.8e\n", R.from, R.to, R.val, R.energy);
    fclose(fl);
}

void RateData::load_csv(const std::filesystem::path& root, std::vector<Rate>& rates){
    std::string fname(root + ".csv");
    std::ifstream ifs(fname);
    if (!ifs.good()){
        throw std::runtime_error("File Not Found");
    }
    std::string line;
    while (getline(fname, line)){
        // ignore comments
        if (line[0] == '#') continue;
        RateData::Rate R;
        std::istringstream iss(line);
        iss >> R.from >> R.to >> R.val >> R.energy;
        rates.push_back(R);
    }
}

void RateData::save_csv(const std::filesystem::path& root, const std::vector<EIIdata>& rates){
    std::string fname(root + ".csv");
    std::ofstream ofs(fname);
    if (!ofs.good()){
        throw std::runtime_error("File Not Found");
    }
    // format it
    ofs << std::setprecision(8);

    fprintf(fl, "# from to occ ionB(Ha) kin(Ha)\n");
    for (auto&R : rates) {
        for (size_t j=0; j<R.fin.size(); j++){
            ofs << std::setw(4) <<  R.init << R.fin[j] << R.occ[j];
            ofs << std::setw(10) << R.ionB[j] << R.kin[j] << "\n"; // fuck windows, all my homies hate windows
        }
    }
    fclose(fl);
}

void RateData::load_csv(const std::filesystem::path& root, std::vector<EIIdata>& rates){
    std::string fname(root + ".csv");
    std::ifstream ifs(fname);
    if (!ifs.good()){
        throw std::runtime_error("File Not Found");
    }

    std::string line;

    while (getline(fname, line)){
        // ignore comments
        if (line[0] == '#') continue;
        RateData::EIIdata eii;
        // temporary variables
        int init, fin, occ;
        double ionB, kin;

        std::istringstream iss(line);
        iss >> init >> fin >> occ >> ionB >> kin;

        size_t l = rates.size();
        if (l > 0 && rates[l-1].init == init) {
            // append to the existing last entry
            rates[l-1].push_back(fin, occ, ionB, kin);
        } else {
            // new EIIData enetry
            EIIdata R;
            R.init = init;
            R.push_back(fin, occ, ionB, kin);
            rates.push_back(R)
        }
        rates.push_back(R);
    }
}

void RateData::save_csv(const std::filesystem::path& root, const RateData::Atom& a){
    save_csv(root + a.name + ".photo", a.Photo);
    save_csv(root + a.name + ".auger", a.Auger);
    save_csv(root + a.name + ".fluor", a.Fluor);
    save_csv(root + a.name + ".eii", a.EIIparams);
}

void RateData::load_csv(const std::filesystem::path& root, Atom& a){
    load_csv(root + a.name + ".photo", a.Photo);
    load_csv(root + a.name + ".auger", a.Auger);
    load_csv(root + a.name + ".fluor", a.Fluor);
    load_csv(root + a.name + ".eii", a.EIIparams);
}