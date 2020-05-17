#include "ElectronSolver.h"
#include "HartreeFock.h"
#include "RateEquationSolver.h"
#include <fstream>
#include <sys/stat.h>

void PhotonFlux::set_parameters(double fluence, double fwhm){
    // The photon flux model
    // Gaussian model A e^-t^2/B
    B = fwhm*fwhm/(4*0.6931471806); // f^2/4ln(2)
    A = fluence/pow(Constant::Pi*B,0.5);
}
inline double PhotonFlux::operator()(double t){
    // Returns flux at time t
    return A*exp(-t*t/B);
};

Weight::Weight(size_t _size){
    size = _size;
    W = (double*) malloc(sizeof(double)*size*size); // allocate the memory
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

ElectronSolver::ElectronSolver(const char* filename, ofstream& log) :
    MolInp(filename, log), pf(width, 100000*fluence)
{
    timespan = this->width*10;
}

state_type ElectronSolver::get_ground_state() {
    state_type initial_condition;
    initial_condition = 0; // May not be necessary, but probably not a bad idea.
    for (size_t a; a<initial_condition.atomP.size(); a++) {
        initial_condition.atomP[a][0] = this->Store[a].nAtoms;
    }
    return initial_condition;
}

void ElectronSolver::compute_cross_sections(std::ofstream& _log, bool recalc) {
    this->calc_rates(_log, recalc);
    hasRates = true;

    // Set up the container class to have the correct size
    state_type::set_elec_points(num_elec_points, min_elec_e, max_elec_e);
    state_type::set_P_shape(this->Store);
    // Set up the rate equations (setup called from parent Adams_BM)
    this->setup(get_ground_state(), this->timespan/num_time_steps);
    cout<<"Using timestep "<<this->dt<<" fs"<<endl;
};

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


// The Big One: Incorporates all of the right hand side to the global
// d/dt P[i] = \sum_i=1^N W_ij - W_ji P[j]
// d/dt f = Q_B[f](t)
void ElectronSolver::sys(const state_type& s, state_type& sdot, const double t){
    for (size_t a = 0; a < s.atomP.size(); a++) {

        // create an appropriately-sized W[i][j]
        Weight W(s.atomP[a].size());

        double J = pf(t); // photon flux in [TODO: UNITS] units
        for ( auto& pht : Store[a].Photo) {
            W(pht.to, pht.from) = pht.val*J;
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
    }
}

// IO functions
void ElectronSolver::save(const std::string& _dir, bool saveSeparate){
    string dir = _dir; // make a copy of the const value
    if (dir.back() == '/') dir.pop_back();

    if (saveSeparate){
        saveFree(dir);
        saveBound(dir+"/boundDist.csv");
    } else {
        saveAll(dir+"/dist.csv");
    }

}

void ElectronSolver::saveAll(const std::string& fname){
    ofstream f;
    cout << "[All] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Electron dynamics"<<endl;
    f << "# Time (fs) [ bound1 | bound2 | bound3 ] free" <<endl;
    assert(y.size() == t.size());
    for (int i=0; i<y.size(); i++){
        f<<t[i]<<" "<<y[i]<<endl;
    }
    f.close();
}

void ElectronSolver::saveFree(const std::string& fname){
    // Saves a table of free-electron dynamics to file fname
    ofstream f;
    cout << "[Free] Saving to file "<<fname<<"..."<<endl;
    f.open(fname);
    f << "# Free electron dynamics"<<endl;
    f << "# Time (fs) | Density [ UNITS ] @ energy:" <<endl;
    f << "#           | ";
    for (int i=0; i<num_elec_points; i++) {
        f<<state_type::f_grid[i]<<" ";
    }
    f<<endl;
    assert(y.size() == t.size());
    for (int i=0; i<y.size(); i++){
        f<<t[i]<<", ";
        for (size_t j = 0; j < y[i].f.size(); j++) {
            f<<y[i].f[j]<<", ";
        }
        f<<endl;
    }
    f.close();
}

void ElectronSolver::saveBound(const std::string& dir){
    // saves a table of bound-electron dynamics , split by atom, to folder dir.
    assert(y.size() == t.size());
    // Iterate over atom types
    for (size_t a=0; a<Store.size(); a++) {
        ofstream f;
        string fname = dir+"/dist_"+Store[a].name+".csv";
        cout << "[Bound] saving to file "<<fname<<"..."<<endl;
        f.open(fname);
        f << "# Ionic electron dynamics"<<endl;
        f << "# Time (fs) | State occupancy (Probability times number of atoms)" <<endl;
        f << "#           | States:";
        // Index, Max_occ inherited from MolInp
        for (auto& cfgname : Store[a].index_names){
            f << cfgname << " ";
        }
        f<<endl;
        // Iterate over time.
        for (size_t i=0; i<y.size(); i++){
            // Make sure all "natom-dimensioned" objects are the size expected
            assert(Store.size() == y[i].atomP.size());
            f<<t[i];
            for (auto& state_prob : y[i].atomP[a]) {
                f<<", "<<state_prob;
            }
            f<<endl;
        }
        f.close();
    }

}
