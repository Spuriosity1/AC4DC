#include <iostream>
#include <fstream>
#include "src/FreeDistribution.h"
#include "src/Constant.h"


using namespace std;

class BasisTester : public SplineIntegral{
    public:
    BasisTester(size_t F_size, double min_e, double max_e, GridSpacing grid_type) {
        
        Distribution::set_elec_points(F_size, min_e, max_e, grid_type);
        // HACK: There are two distinct BasisSet-inheriting things, the Distribution static BasisIntegral
        // and this object. 
        this->set_parameters(F_size, min_e, max_e, grid_type);
    }

    void check_ee(string fname)
    {
        
        // Index runs along final states
        double dQ = 0;        

        ofstream qout(fname);

        // Make Q_mat the correct size to store all the summed Q^J_{KL} coefficients 
        std::vector<std::vector<std::vector<double>>> Q_mat;
        
        Q_mat.resize(Distribution::size);
        for (auto&& v : Q_mat) {
            v.resize(Distribution::size);
            for (auto& vi : v) {
                vi.resize(Distribution::size, 0);
            }
        }
        
        // Read all QEE params into file
        SplineIntegral::pair_list ls;
        for (size_t J = 0; J < Distribution::size; J++) {
            for (size_t K=0; K < Distribution::size; K++) {
                ls = calc_Q_ee(J, K);
                for (auto&& entry : ls) {
                    Q_mat[K][entry.idx][J] = entry.val;
                }
            }
        } 

        // Print some helpful labels
        
        cout<<"#L= ";
        for (size_t L=0; L<Distribution::size; L++) {
            cout <<(supp_max(L)+supp_min(L))*0.5<<" ";
        }
        cout<<endl;
        
        // Print the whole QEE tensor to the supplied file
        for (size_t J = 0; J < Distribution::size; J++)
        {
            qout<<"J = "<<J<<endl;
            for (size_t K = 0; K < Distribution::size; K++)
            {
                for (size_t L = 0; L < Distribution::size; L++)
                {
                    qout<<(Q_mat[K][L][J]+ Q_mat[L][K][J])/2.<<" ";
                }
                qout<<endl;
            }
        }
        qout.close();
        

        for (size_t K = 0; K < Distribution::size; K++)
        {
            for (size_t L=0; L < Distribution::size; L++) {

                // Aggregate charges in each (K,L) combination, summed over J.
                // This should sum to zero.
                // = J_KL (0) - J_KL(max)
                double dQ_tmp=0;
                Eigen::VectorXd v(Distribution::size);
                Eigen::VectorXd w(Distribution::size);
                for (size_t J = 0; J < Distribution::size; J++)
                {
                    v[J] = (Q_mat[K][L][J] + Q_mat[L][K][J])/ 2;
                    w[J] = areas[J];
                }
                cout<<" "<<this->Sinv(v).dot(w);
            }
            cout<<endl;
        }
    }
};


int main(int argc, char const *argv[]) {
    if (argc < 6) cerr << "usage: "<<argv[0]<<" [filestem] [min_e (Ha)] [max_e (Ha)] [num_basis] [ grid type = [l]inear, [q]uadratic, [e]xponential ]";
    string fstem(argv[1]);
    double min_e = atof(argv[2]);
    double max_e = atof(argv[3]);
    int N = atoi(argv[4]);
    GridSpacing gt;
    gt.zero_degree_0=0;
    gt.zero_degree_inf = 3;
    istringstream is(argv[5]);
    is >> gt;
    BasisTester bt(N, min_e, max_e, gt);
    bt.check_ee(fstem);
    
    return 0;
}
