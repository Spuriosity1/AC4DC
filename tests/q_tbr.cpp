#include <iostream>
#include <fstream>
#include "src/FreeDistribution.h"
#include "src/Constant.h"


using namespace std;

RateData::EIIdata get_fake_eii(){
    RateData::EIIdata tmp;
    tmp.init = 0;
    tmp.push_back(5, 2, 40, 8.9);
    // tmp.push_back(1, 1, 1.3, 0.9);
    // tmp.push_back(3, 2, 20, 1.9);
    // tmp.push_back(7, 1, 5, 0.9);
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

    void check_eii(string fstem)
    {
        
        RateData::InverseEIIdata tbr_process = get_fake_eii();
        // Index runs along final states
        std::vector<double> tot_gamma;
        int num_fin_states = tbr_process.size();
        tot_gamma.resize(num_fin_states, 0);
        std::vector<SparsePair> eg;
        double dQ = 0;

        string gamma_loc = fstem + "_gamma.csv";
        string Q_loc = fstem + "_Q.csv";

        ofstream gout(gamma_loc);
        ofstream qout(Q_loc);

        std::vector<std::vector<double>> Q_mat;
        Q_mat.resize(Distribution::size);
        for (auto&& v : Q_mat){
            v.resize(Distribution::size, 0);
        }
        
        // Aggregate charges, summed over J.
        SplineIntegral::sparse_matrix sm;
        for (size_t J = 0; J < Distribution::size; J++)
        {
            sm = calc_Q_tbr(tbr_process, J);
            for (auto&& entry : sm){
                Q_mat[entry.K][entry.L] += entry.val;
            }
        }

        gout<<"#L= ";
        qout<<"#L= ";
        for (size_t L=0; L<Distribution::size; L++){
            gout <<(supp_max(L)+supp_min(L))*0.5<<" ";
            qout <<(supp_max(L)+supp_min(L))*0.5<<" ";
        }
        gout<<endl;
        qout<<endl;

        for (size_t K = 0; K < Distribution::size; K++)
        {
            gout <<(supp_max(K)+supp_min(K))*0.5;
            qout <<(supp_max(K)+supp_min(K))*0.5;
            for (size_t L=0; L<Distribution::size; L++){
                Gamma_tbr(eg, tbr_process, K, L);
                double tot = 0;
                for (size_t j = 0; j < num_fin_states; j++)
                {
                    tot_gamma[j] += eg[j].val;
                    tot += eg[j].val;
                }

                double dQ_tmp = Q_mat[K][L];
                gout<<" "<<tot;
                qout<<" "<<dQ_tmp;
                dQ += dQ_tmp;
            }
            gout<<endl;
            qout<<endl;
        }
        gout.close();
        qout.close();

        cerr<<endl<<"init -> fin   gamma_total   "<<endl;
        double tot=0;
        for (size_t j = 0; j < num_fin_states; j++)
        {
            cerr<<"[ "<<tbr_process.init<<"->"<<tbr_process.fin[j]<<" ] "<<tot_gamma[j]<<endl;
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
    if (argc < 6) cerr << "usage: "<<argv[0]<<" [filestem] [min_e (Ha)] [max_e (Ha)] [num_basis] [ grid type = [l]inear, [q]uadratic, [e]xponential ]";
    double density = 1;
    string fstem(argv[1]);
    double min_e = atof(argv[2]);
    double max_e = atof(argv[3]);
    int N = atoi(argv[4]);
    GridSpacing gt;
    istringstream is(argv[5]);
    is >> gt;
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_eii(fstem);
    
    return 0;
}