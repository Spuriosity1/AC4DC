#include "FreeDistribution.h"
#include "Constant.h"
#include <iostream>

using namespace std;

RateData::EIIdata get_fake_eii(){
    RateData::EIIdata tmp;
    tmp.init = 0;
    tmp.push_back(3, 2, 20, 1.9);
    tmp.push_back(5, 2, 44.3, 8.9);
    tmp.push_back(1, 1, 1.3, 0.9);
    tmp.push_back(7, 1, 5, 0.9);
    return tmp;
}

const auto& one = [](double e) -> double {return 1.;};

class BasisTester : public SplineIntegral{
    public:
    BasisTester(size_t F_size, double min_e, double max_e, GridSpacing grid_type){
        
        Distribution::set_elec_points(F_size, min_e, max_e, grid_type);
        // HACK: There are two distinct BasisSet-inheriting things, the Distribution static BasisIntegral
        // and this object. 
        this->set_parameters(F_size, min_e, max_e, 0, grid_type);
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
            size_t J0 = this->i_from_e(e);
            double tmp =0;
            for (int J = J0 - BSPLINE_ORDER; J<= J0 + BSPLINE_ORDER; J++){
                tmp += this->at(J, e);
            }
            avg += tmp;
            if (fabs(tmp - 1) > 1e-8) cerr<<"[ WARN ] sum of basis funcs @ e="<<e<<" differs from 1:  "<<tmp<<endl;
            e += de;
        }
        avg /= num_funcs * 100;
        cerr<<"Average basis height: "<<avg<<endl<<endl;
    }

    void check_eii()
    {
        
        RateData::EIIdata eii_process = get_fake_eii();
        std::vector<double> tot_gamma;
        int num_fin_states = eii_process.size();
        tot_gamma.resize(num_fin_states, 0);
        std::vector<SparsePair> eg;
        double dQ = 0;

        cout <<"K energy sum_eta gamma_K   sum_J Q_K "<<endl;
        for (size_t K = 0; K < Distribution::size; K++)
        {
            Gamma_eii(eg, eii_process, K);
            double tot = 0;
            for (size_t j = 0; j < num_fin_states; j++)
            {
                tot_gamma[j] += eg[j].val;
                tot += eg[j].val;
            }

            double dQ_tmp = 0;
            for (size_t J = 0; J < Distribution::size; J++)
            {
                dQ_tmp += calc_Q_eii(eii_process, J, K);
            }

            cout << K <<" "<<(supp_max(K)+supp_min(K))*0.5<<" "<<tot <<" "<< dQ_tmp <<endl;
            dQ += dQ_tmp;

        }
        cerr<<endl<<"init -> fin   gamma_total   "<<endl;
        double tot=0;
        for (size_t j = 0; j < num_fin_states; j++)
        {
            cerr<<"[ "<<eii_process.init<<"->"<<eii_process.fin[j]<<" ] "<<tot_gamma[j]<<endl;
            tot += tot_gamma[j];
        }
        cerr<<tot<<endl;
        
        
        

        cerr<<" total dQ/dP_init:"<<endl;
        cerr<< dQ <<endl;
        cerr<<" ratio: "<<endl;

        cerr<<"Fraction: " << tot/ dQ <<endl;

        
    }
};


int main(int argc, char const *argv[]) {
    if (argc < 5) cerr << "usage: "<<argv[0]<<" [min_e (Ha)] [max_e (Ha)] [num_basis] [grid type = [l]inear, [q]uadratic, [e]xponential ]";
    double density = 1;
    double min_e = atof(argv[1]);
    double max_e = atof(argv[2]);
    int N = atoi(argv[3]);
    GridSpacing gt;
    switch (argv[4][0])
    {
    case 'l':
        gt = GridSpacing::linear;
        break;
    case 'q':
        gt = GridSpacing::quadratic;
        break;
    case 'e':
        gt = GridSpacing::exponential;
        break;
    default:
        cerr<<"Unrecognised grid type: '"<<argv[4]<<"'"<<endl;
        return 1;
        break;
    }
    
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_basis();
    // Distribution F;
    // F.set_maxwellian(T, density);
    // cerr<<"Integral of F: should be "<<density<<", is "<<F.integral(one)<<endl;
    // F.addDeltaSpike(T, density);

    // cout<<"#Distribution shape: "<<endl;
    // cout<<Distribution::output_energies_eV(Distribution::size*10)<<endl;
    // cout<<F.output_densities(Distribution::size*10)<<endl;
    
    bt.check_eii();
    
    return 0;
}
