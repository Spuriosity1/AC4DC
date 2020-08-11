#include "ElectronSolver.h"
#include "HartreeFock.h"
#include "RateEquationSolver.h"
#include <fstream>
#include <sys/stat.h>
#include <algorithm>


void PhotonFlux::set_parameters(double fluence, double fwhm){
    // The photon flux model
    // Gaussian model A e^-t^2/B
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
    f << "# Free electron dynamics"<<endl;
    f << "# Time (fs) | Intensity (J/cm2)" <<endl;
    for (auto& t : Tvec){
        f << t << " "<< this->operator()(t)<<std::endl;
    }
    f.close();
}


ElectronSolver::ElectronSolver(const char* filename, ofstream& log) :
    MolInp(filename, log), pf(width*Constant::fs_in_au, 100000*fluence*Constant::Fluence_in_au)
{
    timespan = this->width*10*Constant::fs_in_au;
    precompute_gamma_coeffs();
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
    this->setup(get_ground_state(), this->timespan/num_time_steps);
    cout<<"Using timestep "<<this->dt<<" fs"<<endl;
}

void ElectronSolver::solve(){
    if (!hasRates) {
        std::cerr <<
"No rates have been loaded into the solver. \
Use ElectronSolver::compute_cross_sections(log)\n" << endl;
        return;
    }
    bool converged = false;
    // TODO: repeat this until convergence.
    this->iterate(-timespan/2, this->num_time_steps); // Inherited from ABM
}

void ElectronSolver::precompute_gamma_coeffs(){
    size_t N = Distribution::size;
    size_t gamma_dim = Store[a].Photo.size();
    for (size_t a = 0; a < s.atomP.size(); a++) {
        for (size_t k=0; k<N*(N-1)/2; k++){
            RATE_EII[a][k] = new Eigen::SparseMatrix<double>(gamma_dim, gamma_dim);
        }
        for (size_t n=0; n<N; n++){
            RATE_EII[a][n] = new Eigen::SparseMatrix<double>(gamma_dim, gamma_dim);

            for (auto& eii: Store[a].EIIparams){
                Eigen::SparseMatrix<double> Gamma(gamma_dim, gamma_dim);
                Gamma_eii(Gamma, eii, n, a);
                *RATE_EII[a][n] += Gamma;
                for (size_t m=0; m<n; m++){
                    size_t k = (N*(N+1)/2) - (N-m)*(N-m-1)/2 + n - m - 1
                    // k = N-1... N(N+1)/2
                    Gamma_tbr(Gamma, eii, n, m, a);
                    *RATE_TBR[a][k] += Gamma;
                }
                Gamma_tbr(Gamma, eii, n, n, a);
                *RATE_TBR[a][n] += Gamma;
            }
        }
    }
}

// The Big One: Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j]
// d/dt f = Q_B[f](t)
// system
void ElectronSolver::sys(const state_type& s, state_type& sdot, const double t){
    sdot=0;
    for (size_t a = 0; a < s.atomP.size(); a++) {

        // create an appropriately-sized W[i][j]
        Eigen::MatrixXd W(s.atomP[a].size(), s.atomP[a].size());

        // PHOTOIONISATION
        double J = pf(t); // photon flux in uhhhhh [TODO: UNITS]
        for ( auto& r : Store[a].Photo) {
            W(r.to, r.from) += r.val*J;
            sdot.F.addDeltaLike(r.energy, r.val*J);
        }

        // FLUORESCENCE
        for ( auto& r : Store[a].Fluor) {
            W(r.to, r.from) += r.val;
            // Energy from optical photon assumed lost
        }

        // AUGER
        for ( auto& r : Store[a].Auger) {
            W(r.to, r.from) += r.val;
            sdot.F.addDeltaLike(r.energy, r.val);
        }

        // EII / TBR
        double N = Distribution::size;
        for (size_t n=0; n<N; n++){
            W += RATE_EII[a][n]*s.F.f[n];
            // exploit the symmetry: goofy indexing to only store the upper triangular part.
            for (size_t m=0; m<n; m++){
                size_t k = (N*(N+1)/2) - (N-m)*(N-m-1)/2 + n - m - 1
                // k = N-1... N(N+1)/2
                W += RATE_TBR[a][k]*s.F.f[n]*s.F.f[m]*2;
            }
            // the diagonal
            W += RATE_TBR[a][n]*s.F.f[n]*s.F.f[n];
        }



        const state_type::bound_t& P = s.atomP[a];
        state_type::bound_t& Pdot = sdot.atomP[a];

        // Compute the changes in P
        for (size_t i = 0; i < P.size(); i++) {
            Pdot[i] = 0;
            for (size_t j = 0; j < i; j++) {
                Pdot[i] += W(i, j) * P[j] - W(j, i) * P[i];
            }
            // Avoid j=i
            for (size_t j = i+1; j < P.size(); j++) {
                Pdot[i] += W(i, j) * P[j] - W(j, i) * P[i];
            }

        }

        // Changes in F
        sdot.F.add_Qeii(a, s.F, P);
        sdot.F.add_Qtbr(a, s.F, P);

    }
    sdot.F.add_Qee(s.F); // Electron-electon repulsions

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
    f << "# Time (fs) | Density [ UNITS ] @ energy:" <<endl;
    f << "#           | "<<Distribution::get_energies_eV()<<endl;

    assert(y.size() == t.size());
    for (int i=0; i<y.size(); i++){
        f<<t[i]<<y[i].F<<endl;
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
            f<<t[i] << ' ' << y[i].atomP[a]<<endl;
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
