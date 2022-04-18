/**
 * @file calc_rates.cpp
 * @author Alaric Sanders (als217@cam.ac.uk)
 * @brief Calculates Photoionisation, Auger decay and Impact Ionisation rates
 * @version 0.1
 * @date 2022-04-15
 * 
 * @copyright Copyright (c) 2022
 * 
 */

#include <vanity.hpp>
#include <string>
#include <fstream>
#include <filesystem>

#include "numerical/Constant.h"
#include "hartreefock/HartreeFock.h"
#include "hartreefock/ComputeRateParam.h"
#include "hartreefock/HFInput.h"


// Usage: 
// calc_rates /path/to/atom.toml /path/to/outdir PhotonEnergy(eV)
int main (int argc, char**argv){
    std::cout<<"\033[1m";
    print_file("config/banner.txt");
    std::cout<<"\033[34m";
    print_file("config/version.txt");
    std::cout<<"\033[0m"<<"\n"<<"\n";

    std::cout<<"Copyright (C) 2022  Alaric Sanders and Alexander Kozlov"<<"\n";
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY."<<"\n";
    std::cout<<"This is free software, and you are welcome to redistribute it"<<"\n";
    std::cout<<"under the conditions of the GPLv3; see LICENSE for details."<<"\n";

    // check if we have the right # args
    if (argc < 3) {
        std::cerr<<"Usage: calc_rates /path/to/atom.toml /path/to/outdir wavelength(Angstrom)"
        std::cerr<<std::endl;
        throw std::runtime_error("Incorrect usage");
    }

    // check if the supplied file exists
    std::filesystem::path infile = argv[1];
    std::ifstream in(infile);
    if (!in.good())){
        std::cerr<<"Could not read file "<<argv[1]<<std::endl;
        throw std::runtime_error("Bad infile");
    }

    // check if we can write to the output directory
    std::filesystem::path stem = argv[2];
    // append to path
    stem /= infile.replace_extension() + std::string("_") + argv[3] 
    
    std::ofstream log(stem + ".log");
    if (!log.good())){
        std::cerr<<"Could not access directory "<<argv[2]<<std::endl;
        throw std::runtime_error("Bad outdir");
    }

    double wavelength = atof(argv[3])/Constant::eV_per_Ha;

    //////////////////////////////////////////////
    // All input params validated
    // Time to calculate
    std::vector<int> final_occ(Orbits[a].size(), 0);
    std::vector<int> max_occ(Orbits[a].size(), 0);
    for (size_t i = 0; i < max_occ.size(); i++) {
        if (fabs(Orbits[a][i].Energy) > Omega()) final_occ[i] = Orbits[a][i].occupancy();
        max_occ[i] = Orbits[a][i].occupancy();
    }
    RateData::Atom rates = Dynamics.SolvePlasmaBEB(max_occ, final_occ, log);

    log.close();
    
    RateData::save_csv(stem, rates);
}