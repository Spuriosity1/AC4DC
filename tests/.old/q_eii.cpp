#include <iostream>
#include "src/FreeDistribution.h"
#include "src/Constant.h"
#include <fstream>


using namespace std;

RateData::EIIdata get_fake_eii() {
    RateData::EIIdata tmp;
    tmp.init = 0;
    tmp.push_back(9, 2, 11.3, 1.16);
    tmp.push_back(3, 2, 0.71, 1.94);
    tmp.push_back(1, 2, 1.93, 2.69);
    // tmp.push_back(1, 1, 11, 1);
    return tmp;
}

const auto& one = [](double e) -> double {return 1.;};

class BasisTester : public SplineIntegral{
    public:
    BasisTester(size_t F_size, double min_e, double max_e, GridSpacing grid_type) {
        Distribution::set_elec_points(F_size, min_e, max_e, grid_type);
        // HACK: There are two distinct BasisSet-inheriting things, the Distribution static BasisIntegral
        // and this object. 
        
        this->set_parameters(F_size, min_e, max_e, grid_type);
    }

    void check_eii(string filestem)
    {
        
        RateData::EIIdata eii_process = get_fake_eii();
        std::vector<double> tot_gamma;
        int num_fin_states = eii_process.size();
        tot_gamma.resize(num_fin_states, 0);
        std::vector<SparsePair> eg;
        double dQ = 0;
        ofstream out(filestem + "_qeii.csv");
        ofstream mat_out(filestem + "qeii_matrix.csv");

        out <<"#K energy sum_eta gamma_K   sum_J Q_K "<<endl;
        for (size_t K = 0; K < Distribution::size; K++)
        {
            Gamma_eii(eg, eii_process, K);
            double tot = 0;
            for (size_t j = 0; j < eg.size(); j++)
            {
                tot_gamma[j] += eg[j].val;
                tot += eg[j].val;
            }

            double dQ_tmp = 0;
            mat_out<<K;
            for (size_t J = 0; J < Distribution::size; J++)
            {
                double QKJ = calc_Q_eii(eii_process, J, K);
                mat_out<<" "<<QKJ;
                dQ_tmp += QKJ;
            }
            mat_out<<endl;

            out << K <<" "<<(supp_max(K)+supp_min(K))*0.5<<" "<<tot <<" "<< dQ_tmp <<endl;
            dQ += dQ_tmp;

        }
        out.close();
        mat_out.close();
        cout<<endl<<"init -> fin   gamma_total   "<<endl;
        double tot=0;
        for (size_t j = 0; j < num_fin_states; j++)
        {
            cout<<"[ "<<eii_process.init<<"->"<<eii_process.fin[j]<<" ] "<<tot_gamma[j]<<endl;
            tot += tot_gamma[j];
        }
        cout<<tot<<endl;

        cout<<" total dQ/dP_init:"<<endl;
        cout<< dQ <<endl;
        cout<<" ratio: "<<endl;

        cout<<"Fraction: " << tot/ dQ <<endl;

        
    }
};


int main(int argc, char const *argv[]) {
    if (argc < 6) cerr << "usage: "<<argv[0]<<" [filestem] [min_e (Ha)] [max_e (Ha)] [num_basis] [grid type = [l]inear, [q]uadratic, [e]xponential ]";
    double density = 1;
    string filestem(argv[1]);
    double min_e = atof(argv[2]);
    double max_e = atof(argv[3]);
    int N = atoi(argv[4]);
    GridSpacing gt;
    istringstream inp(argv[5]);
    inp >> gt;
    gt.num_low = N/3;
    gt.transition_e = max_e/5;
    gt.zero_degree_0 = 0;
    gt.zero_degree_inf = 4;
    
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_eii(filestem);
    
    return 0;
}
