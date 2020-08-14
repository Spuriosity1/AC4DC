#include "ElectronSolver.h"
#include "HartreeFock.h"
#include "ComputeRateParam.h"
#include <fstream>
#include <sys/stat.h>
#include <algorithm>
#include <bspline.hpp>
#include "transitionrate.hpp"
#include <eigen3/Eigen/SparseCore>


void PhotonFlux::set_parameters(double fluence, double fwhm){
    // The photon flux model
    // Gaussian model A e^{-t^2/B}
    cout<<"[ Flux ] fluence="<<fluence<<", fwhm="<<fwhm<<endl;
    B = fwhm*fwhm/(4*0.6931471806); // f^2/4ln(2)
    A = fluence/pow(Constant::Pi*B,0.5);
}

inline double PhotonFlux::operator()(double t){
    // Returns flux at time t (same units as fluence)
    return A*exp(-t*t/B);
}

void PhotonFlux::save(const vector<double>& Tvec, const std::string& fname){
    ofstream f;
    cout << "[ Flux ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Time (fs) | Flux (pht/cm2/s)" <<endl;
    for (auto& t : Tvec){
        f << t*Constant::fs_per_au << " ";
        double intensity = (*this)(t)*1e15;
        intensity /= Constant::fs_per_au*Constant::cm_per_au*Constant::cm_per_au;
        f << intensity <<std::endl;
    }
    f.close();
}


ElectronSolver::ElectronSolver(const char* filename, ofstream& log) :
    MolInp(filename, log), pf(fluence, width), Adams_BM(3) // 3rd order Adams method!
{
    timespan_au = this->width*10;
}

state_type ElectronSolver::get_ground_state() {
    state_type initial_condition;
    initial_condition *= 0; // May not be necessary, but probably not a bad idea.
    assert(initial_condition.atomP.size() == Store.size());
    for (size_t a=0; a<Store.size(); a++) {
        initial_condition.atomP[a][0] = Store[a].nAtoms;
    }
    return initial_condition;
}

void ElectronSolver::compute_cross_sections(std::ofstream& _log, bool recalc) {
    this->calc_rates(_log, recalc);
    hasRates = true;
    // Set up the container class to have the correct size
    Distribution::set_elec_points(num_elec_points, min_elec_e, max_elec_e);
    state_type::set_P_shape(this->Store);
    // Set up the rate equations (setup called from parent Adams_BM)
    this->setup(get_ground_state(), this->timespan_au/num_time_steps);
    // create the tensor of coefficients
    RATE_EII.resize(Store.size());
    RATE_TBR.resize(Store.size());
    for (size_t a=0; a<Store.size(); a++){
        RATE_EII[a].resize(num_elec_points);
        RATE_TBR[a].resize(num_elec_points*(num_elec_points+1)/2);
        size_t dim = state_type::P_size(a);
        for (auto& W : RATE_EII[a]){
            W.resize(dim, dim);
        }
        for (auto& W : RATE_TBR[a]){
            W.resize(dim, dim);
        }
    }
    precompute_gamma_coeffs();
}

void ElectronSolver::solve(){
    cout<<"Using timestep "<<this->dt*Constant::fs_per_au<<" fs"<<endl;
    if (!hasRates) {
        std::cerr <<
"No rates have been loaded into the solver. \
Use ElectronSolver::compute_cross_sections(log)\n" << std::endl;
        return;
    }
    bool converged = false;
    // TODO: repeat this until convergence.
    this->iterate(-timespan_au/2, this->num_time_steps); // Inherited from ABM
}

void ElectronSolver::precompute_gamma_coeffs(){
    std::cout<<"[rate precalc] Beginning coefficient computation..."<<std::endl;
    size_t N = Distribution::size;
    for (size_t a = 0; a < Store.size(); a++) {
        std::cout<<"[rate precalc] Atom "<<a+1<<"/"<<Store.size()<<std::endl;
        for (size_t n=0; n<N; n++){
            for (auto& eii: Store[a].EIIparams){
                Eigen::SparseMatrix<double> Gamma;
                RATE_EII[a][n]= Distribution::Gamma_eii(eii, n, a);
                for (size_t m=n+1; m<N; m++){
                    size_t k = (N*(N+1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
                    // k = N-1... N(N+1)/2
                    RATE_TBR[a][k] = Distribution::Gamma_tbr(eii, n, m, a);
                }
                RATE_TBR[a][n] = Distribution::Gamma_tbr(eii, n, n, a);
            }
        }
    }
    std::cout<<"[rate precalc] Done."<<std::endl;
}

// The Big One: Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j]
// d/dt f = Q_B[f](t)
// system
void ElectronSolver::sys(const state_type& s, state_type& sdot, const double t){
    sdot=0;
    for (size_t a = 0; a < s.atomP.size(); a++) {

        // create an appropriately-sized W[i][j]
        GammaType::TransitionRate W(s.atomP[a].size());

        // PHOTOIONISATION
        double J = pf(t); // photon flux in uhhhhh [TODO: UNITS]
        for ( auto& r : Store[a].Photo) {
            W(r.to, r.from) += r.val*J;
            // sdot.F.addDeltaLike(omega - r.energy, r.val*J);
        }

        // FLUORESCENCE
        for ( auto& r : Store[a].Fluor) {
            W(r.to, r.from) += r.val;
            // Energy from optical photon assumed lost
        }

        // AUGER
        for ( auto& r : Store[a].Auger) {
            W(r.to, r.from) += r.val;
            // sdot.F.addDeltaLike(r.energy, r.val);
        }

        // // EII / TBR
        // double N = Distribution::size;
        // for (size_t n=0; n<N; n++){
        //     W += RATE_EII[a][n]*s.F[n];
        //     // exploit the symmetry: goofy indexing to only store the upper triangular part.
        //     for (size_t m=n+1; m<N; m++){
        //         size_t k = (N*(N+1)/2) - (N-n)*(N-n-1)/2 + m - n - 1;
        //         // k = N-1... N(N+1)/2
        //         W += RATE_TBR[a][k]*s.F[n]*s.F[m]*2;
        //     }
        //     // the diagonal
        //     W += RATE_TBR[a][n]*s.F[n]*s.F[n];
        // }



        const bound_t& P = s.atomP[a];
        bound_t& Pdot = sdot.atomP[a];

        Pdot = W.calc_delta(P);


        // Changes in F
        // sdot.F.add_Qeii(a, s.F, P);
        // sdot.F.add_Qtbr(a, s.F, P);

    }
    sdot.F.add_Qee(s.F); // Electron-electon repulsions

    #ifndef NDEBUG
    if (isnan(s.norm())){
        cerr<< "nan encountered in state!"<<endl;
        cerr<< "t = "<<t<<"au = "<<t*Constant::fs_per_au<<"fs"<<endl;
        cerr<< s <<endl;
        throw runtime_error("bad state");
    }
    if (isnan(sdot.norm())){
        cerr<< "nan encountered in state derivative!"<<endl;
        cerr<< "t = "<<t<<"au = "<<t*Constant::fs_per_au<<"fs"<<endl;
        cerr<< sdot <<endl;
        throw runtime_error("bad state_dot");
    }
    #endif

}

// IO functions
void ElectronSolver::save(const std::string& _dir){
    string dir = _dir; // make a copy of the const value
    dir = (dir.back() == '/') ? dir : dir + "/";

    saveFree(dir+"freeDist.csv");
    saveBound(dir);

    pf.save(this->t,dir+"intensity.csv");

}

void ElectronSolver::saveFree(const std::string& fname){
    // Saves a table of free-electron dynamics to file fname
    ofstream f;
    cout << "[ Free ] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Time (fs) | Density @ energy (eV):" <<endl;
    f << "#           | "<<Distribution::get_energies_eV()<<endl;

    assert(y.size() == t.size());
    for (int i=0; i<y.size(); i++){
        f<<t[i]*Constant::fs_per_au<<y[i].F<<endl;
    }
    f.close();
}

void ElectronSolver::saveBound(const std::string& dir){
    // saves a table of bound-electron dynamics , split by atom, to folder dir.
    assert(y.size() == t.size());
    // Iterate over atom types
    for (size_t a=0; a<Store.size(); a++) {
        ofstream f;
        string fname = dir+"dist_"+Store[a].name+".csv";
        cout << "[ Atom ] Saving to file "<<fname<<"..."<<endl;
        f.open(fname);
        f << "# Ionic electron dynamics"<<endl;
        f << "# Time (fs) | State occupancy (Probability times number of atoms)" <<endl;
        f << "#           | States: ";
        // Index, Max_occ inherited from MolInp
        for (auto& cfgname : Store[a].index_names){
            f << cfgname << " ";
        }
        f<<endl;
        // Iterate over time.
        for (size_t i=0; i<y.size(); i++){
            // Make sure all "natom-dimensioned" objects are the size expected
            assert(Store.size() == y[i].atomP.size());
            f<<t[i]*Constant::fs_per_au << ' ' << y[i].atomP[a]<<endl;
        }
        f.close();
    }

}

/*
Weight::Weight(size_t _size){
    size = _size;
    W = (double *) malloc(sizeof(double)*size*size); // allocate the memory
    std::fill(W, W+size*size, 0);
}

Weight::~Weight(){
    free(W);
}

double Weight::from(size_t idx){
    double x=0;
    for (size_t i = 0; i < size; i++) {
        x += operator()(i, idx);
    }
    return x;
}

double Weight::to(size_t idx){
    double x=0;
    for (size_t i = 0; i < size; i++) {
        x += operator()(idx, i);
    }
    return x;
}
*/
