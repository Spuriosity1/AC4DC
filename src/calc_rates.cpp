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
#include <HFInput.hpp>
#include <string>
#include <fstream>
#include <filesystem>

#include "numerical/Constant.hpp"
#include "hartreefock/HartreeFock.hpp"
#include "hartreefock/ComputeRateParam.hpp"
#include "hartreefock/HFInput.hpp"


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
    if (argc < 4) {
        std::cerr<<"Usage: calc_rates /path/to/atom.toml /path/to/outdir wavelength(Angstrom) [num_threads]"<<std::endl;
        throw std::runtime_error("Incorrect usage");
    }

    // check if the supplied file exists
    std::filesystem::path infile = argv[1];
    std::ifstream in(infile);
    
    if (!in.is_open()){
        std::cerr<<"Could not read file "<<argv[1]<<std::endl;
        throw std::runtime_error("Bad infile");
    }


    // check if we can write to the output directory
    std::filesystem::path stem = argv[2];
    // append to path
    stem /= infile.filename().replace_extension();
    stem += std::string("_") + argv[3] + "A";
    stem.replace_extension(".log");

    std::ofstream log(stem);
    if (!log.is_open()){
        std::cerr<<"Could not access directory "<<stem<<std::endl;
        throw std::runtime_error("Bad outdir");
    }

    double wavelength = atof(argv[3])/Constant::eV_per_Ha;

    //////////////////////////////////////////////
    // All input params validated
    // Time to calculate

    HFInput input;
    input.from_toml(in);

    std::vector<int> final_occ(input.orbitals.size(), 0);
    std::vector<int> max_occ(input.orbitals.size(), 0);
    std::vector<RadialWF> Orbitals;
    for (auto& orb: input.orbitals){
        Orbitals.push_back(RadialWF(input.grid_points));
        Orbitals.back().set_N(orb.n);
        Orbitals.back().set_L(orb.l);//setting L overwrites occupancy with 4L+2
        Orbitals.back().set_occupancy(orb.max_occ);
    }
    
    for (size_t i = 0; i < max_occ.size(); i++) {
        if (fabs(Orbitals[i].Energy) > input.omega) final_occ[i] = Orbitals[i].occupancy();
        max_occ[i] = Orbitals[i].occupancy();
    }

    // the grid to construct everyhting over
    Grid lattice(input.grid_points, input.grid_min/input.Z, input.grid_max, 4);

    Potential U(&lattice, input.Z, input.nuclear_potential);
    
    ComputeRateParam Dynamics(lattice, Orbitals, U, input, 4);

    RateData::Atom rates = Dynamics.SolvePlasmaBEB(max_occ, final_occ, log);

    log.close();
    stem.replace_extension("");
    
    RateData::save_csv(stem, rates);
}