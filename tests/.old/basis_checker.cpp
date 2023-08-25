#include <iostream>
#include "src/FreeDistribution.h"
#include "src/Constant.h"
#include <sstream>

using namespace std;


const auto& one = [](double e) -> double {return 1.;};

class BasisTester : public SplineIntegral{
    public:
    BasisTester(size_t F_size, double min_e, double max_e, GridSpacing grid_type) {
        
        Distribution::set_elec_points(F_size, min_e, max_e, grid_type);
        // HACK: There are two distinct BasisSet-inheriting things, the Distribution static BasisIntegral
        // and this object. 
        this->set_parameters(F_size, min_e, max_e, grid_type);

        cout<<"# Grid: ";
        for (auto& e : this->knot)
        {
            cout<<e<<" ";
        }
        cout<<endl;
        
    }

    void check_basis()
    {
        double min_e = this->knot[0];
        double max_e = this->knot[num_funcs];
        double e = min_e;
        double de = (max_e - min_e)/this->num_funcs/100;
        
        cerr<<"Testing basis unity partition: "<<endl;
        double avg =0;
        for (size_t i = 0; i < this->num_funcs*100; i++)
        {
            int J0 = this->i_from_e(e);
            double tmp =0;
            for (int J = J0 - BSPLINE_ORDER; J<= J0 + BSPLINE_ORDER; J++) {
                tmp += this->at(J, e);
            }
            // for (int J=0; J<num_funcs; J++) {
            //     tmp += (*this)(J, e);
            // }
            avg += tmp;
            if (fabs(tmp - 1) > 1e-8) cerr<<"[ WARN ] sum of basis funcs @ e="<<e<<" differs from 1:  "<<tmp<<endl;
            e += de;
        }
        avg /= num_funcs * 100;
        cerr<<"Average basis height: "<<avg<<endl<<endl;
    }

};


int main(int argc, char const *argv[]) {
    if (argc < 6) cerr << "usage: "<<argv[0]<<" [min_e (Ha)] [max_e (Ha)] [num_basis] [electron temperature] [grid style]";
    double density = 17.12;
    double min_e = atof(argv[1]);
    double max_e = atof(argv[2]);
    int N = atoi(argv[3]);
    double T = atof(argv[4]);
    istringstream grid(argv[5]);
    GridSpacing gt;
    grid >> gt;
    
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_basis();
    Distribution F;
    F = 0;
    F.add_maxwellian(density, T);
    Eigen::VectorXd v = Eigen::VectorXd::Zero(Distribution::size);
    F.addDeltaLike(v, T*3, density);
    F.addDeltaSpike(T*2, density);
    F.applyDeltaF(v);
    cerr<<"Integral of F: should be "<<3*density<<", is "<<F.integral(one)<<endl;
    
    cout<<"#Distribution shape: "<<endl;
    cout<<Distribution::output_energies_eV(Distribution::size*10)<<endl;
    cout<<F.output_densities(Distribution::size*10)<<endl;
    
    
    return 0;
}
