/**
 * @file ee_dynamics.cpp
 * @brief 
 * 
 */
#include "HybridIntegrator.hpp"
#include "RateSystem.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// INTERNAL UNITS:
// eV, fs
class EEcoll : ode::Hybrid<Distribution>{
public:
    EEcoll(double min_e, double max_e, size_t num_e_points, const GridSpacing& g) : ode::Hybrid<Distribution>(2) {
        //                            num, min, max
        Distribution::set_elec_points(num_e_points, min_e, max_e, g);
        
        // HACK: give it an empty Store
        vector<RateData::Atom> dummy(0);
        Distribution::precompute_Q_coeffs(dummy);
        
    }

    void set_initial_condition(double E, double density, double spike_density, double _dt) {
        Distribution init;
        cerr<<"Initial conditions: n="<<density<<" Maxwellian T="<<E/4.*Constant::eV_per_Ha<<" eV, ";
        cerr<<"curve with n="<<spike_density<<" delta at "<<E*Constant::eV_per_Ha<<"eV"<<endl;
        init = 0;
        init.add_maxwellian(E/4., density);
        init.addDeltaSpike(E, spike_density);
        this->setup(init, _dt);
        cout<<init.density(Distribution::size)<<endl;
        // this->y and this->t are now sized appropriately
    }

    void run_sim(double final_time){
        try
        {
            this->iterate(0, final_time);
        }
        catch(const std::exception& e)
        {
            std::cerr << e.what() << '\n';
        }
    }

    void printraw(string fname) {
        std::ofstream os(fname);
        os<<"#t | E (eV)"<<endl;
        os<<"#  |  "<<Distribution::output_knots_eV()<<endl;
        for (size_t i = 0; i < y.size(); i++) {
            os<<t[i]*Constant::fs_per_au<<' '<<y[i]<<endl;
        }
    }

    void print(string fname) {
        std::ofstream os(fname);
        os<<"#t | E (eV)"<<endl;
        os<<"#  |  "<<Distribution::output_energies_eV(Distribution::size*10)<<endl;
        for (size_t i = 0; i < y.size(); i++) {
            os<<t[i]*Constant::fs_per_au<<' '<<y[i].output_densities(Distribution::size*10)<<endl;
        }
    }
protected:
    void sys_bound(const Distribution& q, Distribution& qdot, const double t) {
        qdot=0;
    }
    void sys_ee(const Distribution& q, Distribution& qdot) {
        qdot=0;
        Eigen::VectorXd v = Eigen::VectorXd::Zero(Distribution::size);
        q.get_Q_ee(v);
        qdot.applyDeltaF(v);
        if (isnan(qdot.norm())) throw runtime_error("NaN encountered in sdot");
    }
};

int main(int argc, char const *argv[]) {

    if (argc < 10) {
        cerr<<"Usage: "<<argv[0]<<" [fname] [T (eV)] [fin_time (fs)] [num_t_points] [num_e_points] [thermal density (A^-3)] [spike density (A^-3)] [max e (eV)] [gridstyle]"<<endl;
        return 1;
    }
    string fname(argv[1]);
    double temperature = atof(argv[2]);
    double fin_time = atof(argv[3]);
    int num_t_pts = atoi(argv[4]);
    int num_e_pts = atoi(argv[5]);
    double density = atof(argv[6]);
    double sdensity = atof(argv[7]);
    double max_e = atof(argv[8]);
    GridSpacing gs;
    istringstream is(argv[9]);
    is >> gs;
    gs.num_low = num_e_pts/2;
    gs.zero_degree_0=0;
    gs.zero_degree_inf=2;
    

    cerr<<"MB Density ="<<density<<" elec per Angstrom3 @ kT = "<<temperature/4<<"eV"<<endl;
    cerr<<"Spike Density ="<<sdensity<<" elec per Angstrom3 @ E = "<<temperature<<"eV"<<endl;
    cerr<<"Final Time: "<<fin_time<<" fs in "<<num_t_pts<<" timesteps"<<endl;

    fin_time /= Constant::fs_per_au;
    temperature /= Constant::eV_per_Ha;
    max_e /= Constant::eV_per_Ha;
    density *= Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;
    sdensity *= Constant::Angs_per_au*Constant::Angs_per_au*Constant::Angs_per_au;    

    gs.transition_e = max_e/3.;

    double step = fin_time/num_t_pts;
    double min_e = 0;
    EEcoll I(min_e, max_e, num_e_pts, gs);
    cerr<<"final time: "<<fin_time<<" au"<<endl;
    I.set_initial_condition(temperature, density, sdensity, step);
    I.run_sim(fin_time);
    I.print(fname);
    I.print(fname + ".raw");
    return 0;
}
