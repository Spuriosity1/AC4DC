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

// (C) Alaric Sanders 2020

#include "ComputeRateParam.h"
#include "ElectronSolver.h"
#include "Input.h"
#include "Constant.h"
#include <iostream>

using namespace std;

// Rate system solver.
// Uses precomputed rates from AC4DC for all atomic cross-section data.
// KEEP IN MIND:
// - For every atom X listed in the .mol file, AC4DC must be run for the file X.inp
// - AC4DC has input parameters for pulse width, energy and fluence.
// - Only photon energy affects the rate calculations.
// Let scripts/run.py handle all of these details.


void print_banner(const char* fname){
    std::ifstream ifs(fname, ifstream::in);

    char c = ifs.get();
    while (ifs.good()) {
        std::cout << c;
        c = ifs.get();
    }
    ifs.close();
}

void try_mkdir(const std::string& fname) {
    if (mkdir(fname.c_str(), ACCESSPERMS) == -1) {
        if (errno != EEXIST)
            std::cerr<<"mkdir error attempting to create "<< fname << ":" << errno;
    }
}

int get_file_names(const char* infile_, std::string &tag, std::string &logfile, std::string&outdir) {
    // Takes infile of the form "DIR/Lysozyme.mol"
    // Stores "Lysozyme" in tag, "output/log/run_Lysozyme" in logfile
    std::string infile = string(infile_);
    size_t tagstart = infile.rfind('/');
    size_t tagend = infile.rfind('.');
    tagstart = (tagstart==string::npos) ? 0 : tagstart + 1;// Exclude leading slash
    tagend = (tagend==string::npos) ? infile.size() : tagend;
    tag = infile.substr(tagstart, tagend-tagstart);
    // guarantee the existence of a folder structure
    try_mkdir("output");
    try_mkdir("output/log");
    try_mkdir("output/__Molecular");
    logfile = "output/log/run_" + tag + ".log";
    // check correct format
    outdir = "output/__Molecular/"+tag+"/";
    try_mkdir(outdir);
    std::string extension = infile.substr(tagend);
    if (extension != ".mol") {
        std::cerr<<"This file is for coupled calculations. Please provide a .mol file similar to Lysozyme.mol"<<"\n";
        return 1;
    }
    return 0;
}

struct CmdParser{
    CmdParser(int argc, const char *argv[]) {
        if (argc < 2) {
            std::cout << "Usage: solver path/to/molecular/in.mol [-rh]";
            valid_input = false;
        }
        for (int a=2; a<argc; a++) {
            if (argv[a][0] != '-')
                continue;

            int i=1;
            while (argv[a][i]!='\0') {
                switch (argv[a][i]) {
                    case 's':
                        // recalculate
                        std::cout<<"\033[1;35mWarning:\033[0m Searching output folder for precalculated rates."<<"\n";
                        std::cout<<"If these were calculated for a different X-ray energy, the results will be wrong!"<<"\n";
                        recalc = false;
                        break;
                    case 'x':
                        // X - sections only
                        std::cout<<"\033[1;31mSolving for cross-sections only.\033[0m"<<"\n";
                        solve_rate_eq = false;
                        break;
                    case 'h':
                        // Usage help.
                        std::cout<<"This is physics code, were you really expecting documentation?"<<"\n";
                        std::cout<<"  -s Look for stored precalculated rate coefficients"<<"\n";
                        std::cout<<"  -x Skip rate-equaton solving"<<"\n";
                        break;
                    case 'w':
                        // Warranty.
                        std::cout<<"This program is provided as-is, with no warranty, explicit or implied."<<"\n";
                        std::cout<<"It is a simplified model of XFEL plasma dynamics, however it should not"<<"\n";
                        std::cout<<"be viewed as accurate under all conditions."<<"\n";
                        exit(0);
                        break;
                    case 'c':
                        // Copyright.
                        print_banner("LICENSE");
                        exit(0);
                        break;
                    default:
                        std::cout<<"Flag '"<<argv[a][i]<<"' is not a recognised flag."<<"\n";
                }
                i++;
            }
        }
    }
    bool recalc = true;
    bool valid_input = true;
    bool solve_rate_eq = true;
};

int main(int argc, const char *argv[]) {
    CmdParser runsettings(argc, argv);
    if (!runsettings.valid_input) {
        return 1;
    }

    std::cout<<"\033[1m";
    print_banner("config/banner.txt");
    std::cout<<"\033[34m";
    print_banner("config/version.txt");
    std::cout<<"\033[0m"<<"\n"<<"\n";

    std::string name, logname, outdir;

    std::cout<<"Copyright (C) 2020  Alaric Sanders and Alexander Kozlov"<<"\n";
    std::cout<<"This program comes with ABSOLUTELY NO WARRANTY; for details run `ac4dc -w'."<<"\n";
    std::cout<<"This is free software, and you are welcome to redistribute it"<<"\n";
    std::cout<<"under certain conditions; run `ac4dc -c' for details."<<"\n";

    if (get_file_names(argv[1], name, logname, outdir) == 1)
        return 1;

    std::cout<<"Running simulation for target "<<name<<"\n";
    std::cout << "logfile name: " << logname <<"\n";
    std::ofstream log(logname);
    std::cout << "\033[1;32mInitialising... \033[0m" <<"\n";
    ElectronSolver S(argv[1], log); // Contains all of the collision parameters.
    std::cout << "\033[1;32mComputing cross sections... \033[0m" <<"\n";
    S.compute_cross_sections(log, runsettings.recalc);
    if (runsettings.solve_rate_eq) {
        std::cout << "\033[1;32mSolving rate equations... \033[0m" <<"\n";
        S.solve();
        std::cout << "\033[1;32mDone! \033[0m" <<"\n";
        S.save(outdir);
    } else {
        std::cout << "\033[1;32mDone! \033[0m" <<"\n";
    }

    return 0;
}
