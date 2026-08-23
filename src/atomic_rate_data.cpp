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

// Stage 1 of the AC4DC suite: the Hartree-Fock atomic physics.
// Computes the photoionisation / fluorescence / Auger / electron-impact-ionisation
// rates and form factors for a single element at a single photon energy, and stores
// them in <out_dir>/<element>_<omega>eV.h5 for ac4dc (Stage 2) to consume. <out_dir>
// defaults to the directory containing the .inp input file.

#include "ComputeRateParam.h"
#include "HartreeFock.h"
#include "Input.h"
#include "Grid.h"
#include "Potential.h"
#include "RadialWF.h"
#include "Constant.h"
#include "RateHDF5.h"

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
namespace fs = std::filesystem;

static void usage() {
    cout << "Usage: atomic_rate_data <input.inp> <photon_energy_eV> [-o out_dir] [-j num_threads]" << endl;
    cout << "  <input.inp>         path to an atomic input file (e.g. input/atoms/C.inp)" << endl;
    cout << "  <photon_energy_eV>  XFEL photon energy in eV (e.g. 9000)" << endl;
    cout << "  -o out_dir          directory for the output .h5 (default: the input's directory)" << endl;
    cout << "  -j num_threads      OpenMP threads for the HF calculation (default 4)" << endl;
    cout << endl;
    cout << "Writes <out_dir>/<element>_<photon_energy>eV.h5" << endl;
}

int main(int argc, const char* argv[]) {
    if (argc < 3) { usage(); return 1; }

    const string inp_path = argv[1];
    double omega_eV;
    try {
        omega_eV = stod(argv[2]);
    } catch (const exception&) {
        cerr << "Could not parse photon energy '" << argv[2] << "' as a number." << endl;
        usage();
        return 1;
    }
    if (omega_eV <= 0) { cerr << "Photon energy must be positive." << endl; return 1; }

    int num_threads = 4;
    string out_dir; // empty => default to the input's directory
    for (int a = 3; a < argc; ++a) {
        if (string(argv[a]) == "-j" && a + 1 < argc) num_threads = atoi(argv[++a]);
        else if (string(argv[a]) == "-o" && a + 1 < argc) out_dir = argv[++a];
        else if (string(argv[a]) == "-h") { usage(); return 0; }
    }

    if (!fs::exists(inp_path)) {
        cerr << "Atomic input file not found: " << inp_path << endl;
        return 1;
    }

    // Element id = the input's basename with its extension (".inp") dropped.
    const string element = fs::path(inp_path).stem().string();
    if (out_dir.empty()) {
        out_dir = fs::path(inp_path).parent_path().string();
        if (out_dir.empty()) out_dir = ".";
    }

    const double omega_au = omega_eV / Constant::eV_per_Ha;

    std::error_code ec;
    fs::create_directories(out_dir, ec);
    if (ec) cerr << "Could not create output directory " << out_dir
                 << ": " << ec.message() << endl;

    ofstream log(out_dir + "/log_" + element + "_"
                 + to_string((long)llround(omega_eV)) + "eV.txt");

    cout << "\033[1;32m[ atomic_rate_data ]\033[0m computing rates for " << element
         << " at " << omega_eV << " eV (" << num_threads << " threads)" << endl;

    // --- Single-atom setup, mirroring the per-atom loop in MolInp's constructor ---
    Grid Latt(0);
    vector<RadialWF> Orb;
    Input in(const_cast<char*>(inp_path.c_str()), Orb, Latt, log);
    in.Set_Pulse(omega_au, /*fluence*/ 1.0, /*width*/ 1.0); // only omega affects the rates
    in.Set_Num_Threads(num_threads);
    Potential U(&Latt, in.Nuclear_Z(), in.Pot_Model());

    // --- Solve HF, then compute the plasma rate parameters ---
    HartreeFock HF(Latt, Orb, U, in, log);
    ComputeRateParam Dynamics(Latt, Orb, U, in);

    vector<int> final_occ(Orb.size(), 0);
    vector<int> max_occ(Orb.size(), 0);
    vector<bool> shell_check(Orb.size(), false);
    for (size_t i = 0; i < Orb.size(); i++) {
        // Lock in electrons bound more deeply than the photon energy.
        if (fabs(Orb[i].Energy) > in.Omega()) final_occ[i] = Orb[i].occupancy();
        max_occ[i] = Orb[i].occupancy();
        shell_check[i] = Orb[i].is_shell();
    }

    RateData::Atom atom = Dynamics.SolvePlasmaBEB(max_occ, final_occ, shell_check, log);
    atom.name = element; // element identifier used by ac4dc's lookup and stored in the file

    const string out_path = RateHDF5::filename(element, omega_eV, out_dir);
    cout << "\033[1;32m[ atomic_rate_data ]\033[0m writing " << out_path
         << " (" << atom.num_conf << " configurations)" << endl;
    RateHDF5::write(out_path, atom, Dynamics.AllFormFactors(), in.Nuclear_Z(), omega_eV);

    cout << "\033[38;5;47mDone!\033[0m" << endl;
    return 0;
}
