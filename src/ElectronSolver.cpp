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

#include "ElectronSolver.h"
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include "SplineIntegral.h"
#include <fstream>
#include <algorithm>
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <chrono>
#include <math.h>
#include <omp.h>
#include "config.h"


state_type ElectronSolver::get_ground_state() {
    state_type initial_condition;
    assert(initial_condition.atomP.size() == input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        initial_condition.atomP[a][0] = input_params.Store[a].nAtoms;
        for(size_t i=1; i<initial_condition.atomP.size(); i++) {
            initial_condition.atomP[a][i] = 0.;
        }
    }
    initial_condition.F=0;
    // std::cout<<"[ Rate Solver ] initial condition:"<<std::endl;
    // std::cout<<"[ Rate Solver ] "<<initial_condition<<std::endl;
    return initial_condition;
}

void ElectronSolver::get_energy_bounds(double& max, double& min) {
    max = 0;
    min = 1e9;
    for(auto& atom : input_params.Store) {
        for(auto& r : atom.Photo) {
            if (r.energy > max) max=r.energy;
            if (r.energy < min) min=r.energy;
        }
        for(auto& r : atom.Auger) {
            if (r.energy > max) max=r.energy;
            if (r.energy < min) min=r.energy;
        }
    }
}


void ElectronSolver::compute_cross_sections(std::ofstream& _log, bool recalc) {
    input_params.calc_rates(_log, recalc);
    hasRates = true;

    // Set up the container class to have the correct size
    Distribution::set_elec_points(input_params.Num_Elec_Points(), input_params.Min_Elec_E(), input_params.Max_Elec_E(), input_params.elec_grid_type);
    state_type::set_P_shape(input_params.Store);
    // Set up the rate equations (setup called from parent Adams_BM)
    this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, 5e-3);
    // create the tensor of coefficients
    RATE_EII.resize(input_params.Store.size());
    RATE_TBR.resize(input_params.Store.size());
    for (size_t a=0; a<input_params.Store.size(); a++) {
        size_t N = input_params.Num_Elec_Points();
        RATE_EII[a].resize(N);
        RATE_TBR[a].resize(N*(N+1)/2);
    }
    precompute_gamma_coeffs();
    Distribution::precompute_Q_coeffs(input_params.Store);
}

void ElectronSolver::solve() {
    assert (hasRates || "No rates found! Use ElectronSolver::compute_cross_sections(log)\n");
    auto start = std::chrono::system_clock::now();

    good_state = true;
    if (input_params.pulse_shape ==  PulseShape::square){
        this->iterate(-input_params.Width(), timespan_au - input_params.Width()); // Inherited from ABM
    } else {
        this->iterate(-timespan_au/2, timespan_au/2); // Inherited from ABM
    }
    
    
    cout<<"[ Rate Solver ] Using timestep "<<this->dt*Constant::fs_per_au<<" fs"<<std::endl;
    
    double time = this->t[0];
    int retries = 1;
    while (!good_state) {
        std::cerr<<"\033[93;1m[ Rate Solver ] Halving timestep...\033[0m"<<std::endl;
        good_state = true;
        input_params.num_time_steps *= 2;
        if (input_params.num_time_steps > MAX_T_PTS){
            std::cerr<<"\033[31;1m[ Rate Solver ] Exceeded maximum T size. Skipping remaining iteration."<<std::endl;
            break;
        }
        if (timestep_reached - time < 0 || retries == 0){
            std::cerr<<"\033[31;1m[ Rate Solver ] Shorter timestep failed to improve convergence. Skipping remaining iteration."<<std::endl;
            break;
        } else {
            time = timestep_reached;
        }
        retries--;
        this->setup(get_ground_state(), this->timespan_au/input_params.num_time_steps, 5e-3);
        if (input_params.pulse_shape ==  PulseShape::square){
            this->iterate(-input_params.Width(), timespan_au - input_params.Width()); // Inherited from ABM
        } else {
            this->iterate(-timespan_au/2, timespan_au/2); // Inherited from ABM
        }
    }
    
    

    auto end = chrono::system_clock::now();
    chrono::duration<double> elapsed_seconds = end-start;
    time_t end_time = std::chrono::system_clock::to_time_t(end);

    cout << "[ Solver ] finished computation at " << ctime(&end_time) << endl;
    long secs = elapsed_seconds.count();
    cout<<"[ Solver ] ODE iteration took "<< secs/60 <<"m "<< secs%60 << "s" << endl;

}


// Populates RATE_EII and RATE_TBR based on calls to Distribution::Gamma_eii/ Distribution::Gamma_tbr
// 
void ElectronSolver::precompute_gamma_coeffs() {
    std::cout<<"[ Gamma precalc ] Beginning coefficient computation..."<<std::endl;
    size_t N = Distribution::size;
    for (size_t a = 0; a < input_params.Store.size(); a++) {
        std::cout<<"\n[ Gamma precalc ] Atom "<<a+1<<"/"<<input_params.Store.size()<<std::endl;
        auto eiiVec = input_params.Store[a].EIIparams;
        vector<RateData::InverseEIIdata> tbrVec = RateData::inverse(eiiVec);
        size_t counter=1;
        #pragma omp parallel default(none) shared(a, N, counter, tbrVec, eiiVec, RATE_EII, RATE_TBR, std::cout)
		{
			#pragma omp for schedule(dynamic) nowait
            for (size_t n=0; n<N; n++) {
                #pragma omp critical
                {
                    std::cout<<"\r[ Gamma precalc ] "<<counter<<"/"<<N<<" thread "<<omp_get_thread_num()<<std::flush;
                    counter++;
                }
                #ifndef NO_EII
                Distribution::Gamma_eii(RATE_EII[a][n], eiiVec, n);
                #endif
                // Weird indexing exploits symmetry of Gamma_TBR
                // such that only half of the coefficients are stored
                #ifndef NO_TBR
                for (size_t m=n+1; m<N; m++) {
                    size_t k = N + (N*(N-1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                    // // k = N... N(N+1)/2
                    Distribution::Gamma_tbr(RATE_TBR[a][k], tbrVec, n, m);
                }
                Distribution::Gamma_tbr(RATE_TBR[a][n], tbrVec, n, n);
                #endif
                
            }
        }
    }
    std::cout<<"\n[ Gamma precalc ] Done."<<std::endl;
}


// The Big One: The global ODE functon
// Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j]
// d/dt f = Q_B[f](t)

// Non-stiff part of the system
void ElectronSolver::sys(const state_type& s, state_type& sdot, const double t) {
    sdot=0;
    Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);

    if (!good_state) return;

    for (size_t a = 0; a < s.atomP.size(); a++) {
        const bound_t& P = s.atomP[a];
        bound_t& Pdot = sdot.atomP[a];

        // PHOTOIONISATION
        double J = pf(t); // photon flux in atomic units
        for ( auto& r : input_params.Store[a].Photo) {
            double tmp = r.val*J*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(r.energy, r.val*J*P[r.from]);
            sdot.bound_charge +=  tmp;
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*J*P[r.from]);
        }

        // FLUORESCENCE
        for ( auto& r : input_params.Store[a].Fluor) {
            Pdot[r.to] += r.val*P[r.from];
            Pdot[r.from] -= r.val*P[r.from];
            // Energy from optical photon assumed lost
        }

        // AUGER
        for ( auto& r : input_params.Store[a].Auger) {
            double tmp = r.val*P[r.from];
            Pdot[r.to] += tmp;
            Pdot[r.from] -= tmp;
            sdot.F.addDeltaSpike(r.energy, r.val*P[r.from]);
            // sdot.F.add_maxwellian(r.energy*2./3., r.val*P[r.from]);
            // Distribution::addDeltaLike(vec_dqdt, r.energy, r.val*P[r.from]);
            sdot.bound_charge +=  tmp;
        }

        // EII / TBR bound-state dynamics
        size_t N = Distribution::size;
        for (size_t n=0; n<N; n++) {
            double tmp=0; // aggregator
            
            #ifndef NO_EII
            for (size_t init=0;  init<RATE_EII[a][n].size(); init++) {
                for (auto& finPair : RATE_EII[a][n][init]) {
                    tmp = finPair.val*s.F[n]*P[init];
                    Pdot[finPair.idx] += tmp;
                    Pdot[init] -= tmp;
                    sdot.bound_charge += tmp;
                }
            }
            #endif
            
            
            #ifndef NO_TBR
            // exploit the symmetry: strange indexing engineered to only store the upper triangular part.
            // Note that RATE_TBR has the same geometry as InverseEIIdata.
            for (size_t m=n+1; m<N; m++) {
                size_t k = N + (N*(N-1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                // k = N... N(N+1)/2-1
                // W += RATE_TBR[a][k]*s.F[n]*s.F[m]*2;
                for (size_t init=0;  init<RATE_TBR[a][k].size(); init++) {
                    for (auto& finPair : RATE_TBR[a][k][init]) {
                        tmp = finPair.val*s.F[n]*s.F[m]*P[init]*2;
                        Pdot[finPair.idx] += tmp;
                        Pdot[init] -= tmp;
                        sdot.bound_charge -= tmp;
                    }
                }
            }
            // the diagonal
            // W += RATE_TBR[a][n]*s.F[n]*s.F[n];
            for (size_t init=0;  init<RATE_TBR[a][n].size(); init++) {
                for (auto& finPair : RATE_TBR[a][n][init]) {
                    tmp = finPair.val*s.F[n]*s.F[n]*P[init];
                    Pdot[finPair.idx] += tmp;
                    Pdot[init] -= tmp;
                    sdot.bound_charge -= tmp;
                }
            }
            #endif
        }

        // Free-electron parts
        #ifdef NO_EII
        #warning No impact ionisation
        #else
        s.F.get_Q_eii(vec_dqdt, a, P);
        #endif
        #ifdef NO_TBR
        #warning No three-body recombination
        #else
        s.F.get_Q_tbr(vec_dqdt, a, P);
        #endif
    }


    sdot.F.applyDelta(vec_dqdt);
    // This is loss.
    sdot.F.addLoss(s.F, input_params.loss_geometry, s.bound_charge);

    if (isnan(s.norm()) || isnan(sdot.norm())) {
        cerr<<"NaN encountered in ODE iteration."<<endl;
        cerr<< "t = "<<t*Constant::fs_per_au<<"fs"<<endl;
        good_state = false;
        timestep_reached = t*Constant::fs_per_au;
    }
}


// 'badly-behaved' part of the system
void ElectronSolver::sys2(const state_type& s, state_type& sdot, const double t) {
    sdot=0;
    Eigen::VectorXd vec_dqdt = Eigen::VectorXd::Zero(Distribution::size);
    
    // compute the dfdt vector
    
    // for (size_t a = 0; a < s.atomP.size(); a++) {
    //     const bound_t& P = s.atomP[a];  
        
    // }
    
    

    #ifdef NO_EE
    #warning No electron-electron interactions
    #else
    s.F.get_Q_ee(vec_dqdt); // Electron-electon repulsions
    #endif
    sdot.F.applyDelta(vec_dqdt);
}

// IO functions
void ElectronSolver::save(const std::string& _dir) {
    string dir = _dir; // make a copy of the const value
    dir = (dir.back() == '/') ? dir : dir + "/";

    saveFree(dir+"freeDist.csv");
    saveFreeRaw(dir+"freeDistRaw.csv");
    saveBound(dir);

    std::vector<double> fake_t;
    size_t num_t_points = input_params.Out_T_size();
    if ( num_t_points >  t.size() ) num_t_points = t.size();
    size_t t_idx_step = t.size() / num_t_points;
    for (size_t i=0; i<num_t_points; i++) {
        fake_t.push_back(t[i*t_idx_step]);
    }
    pf.save(fake_t,dir+"intensity.csv");

}

void ElectronSolver::saveFree(const std::string& fname) {
    // Saves a table of free-electron dynamics to file fname
    ofstream f;
    cout << "[ Free ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Time (fs) | Density @ energy (eV):" <<endl;
    f << "#           | "<<Distribution::output_energies_eV(this->input_params.Out_F_size())<<endl;

    assert(y.size() == t.size());
    size_t num_t_points = input_params.Out_T_size();
    if ( num_t_points >  t.size() ) num_t_points = t.size();
    size_t t_idx_step = t.size() / num_t_points;
    for (size_t i=0; i<num_t_points; i++) {
        f<<t[i*t_idx_step]*Constant::fs_per_au<<" "<<y[i*t_idx_step].F.output_densities(this->input_params.Out_F_size())<<endl;
    }
    f.close();
}

void ElectronSolver::saveFreeRaw(const std::string& fname) {
    ofstream f;
    cout << "[ Free ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Energy Knot: "<< Distribution::output_knots_eV() << endl;
    f << "# Time (fs) | Expansion Coeffs"  << endl;

    assert(y.size() == t.size());
    
    for (size_t i=0; i<t.size(); i++) {
        f<<t[i]*Constant::fs_per_au<<" "<<y[i].F<<endl;
    }
    f.close();
}

void ElectronSolver::saveBound(const std::string& dir) {
    // saves a table of bound-electron dynamics , split by atom, to folder dir.
    assert(y.size() == t.size());
    // Iterate over atom types
    for (size_t a=0; a<input_params.Store.size(); a++) {
        ofstream f;
        string fname = dir+"dist_"+input_params.Store[a].name+".csv";
        cout << "[ Atom ] Saving to file "<<fname<<"..."<<endl;
        f.open(fname);
        f << "# Ionic electron dynamics"<<endl;
        f << "# Time (fs) | State occupancy (Probability times number of atoms)" <<endl;
        f << "#           | ";
        // Index, Max_occ inherited from MolInp
        for (auto& cfgname : input_params.Store[a].index_names) {
            f << cfgname << " ";
        }
        f<<endl;
        // Iterate over time.
        size_t num_t_points = input_params.Out_T_size();
        if ( num_t_points >  t.size() ) num_t_points = t.size();
        size_t t_idx_step = t.size() / num_t_points;        
        for (size_t i=0; i<num_t_points; i++) {
            // Make sure all "natom-dimensioned" objects are the size expected
            assert(input_params.Store.size() == y[i].atomP.size());
            
            f<<t[i*t_idx_step]*Constant::fs_per_au << ' ' << y[i*t_idx_step].atomP[a]<<endl;
        }
        f.close();
    }

}
