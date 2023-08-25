#include "src/SplineBasis.h"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc < 5) {
        cout<<"Usage: "<<argv[0]<<" [fstem] [num_funcs] [min_e (Ha)] [max_e (Ha)] (linear|quadratic|exponential) "<<endl;
    }
    std::string dummy(argv[1]);
    ofstream fout(dummy+"_vals.csv");
    ofstream dfout(dummy+"_derivative.csv");
    int num_funcs = atoi(argv[2]);
    double min = atof(argv[3]);
    double max = atof(argv[4]);
    istringstream is(argv[5]);
    GridSpacing gs;
    is >> gs;
    BasisSet basis;
    gs.zero_degree_0 = 1;
    gs.zero_degree_inf = 0;
    basis.set_parameters(num_funcs, min, max, gs);
    cout<<"Testing partition of unity"<<endl;
    double step = (max - min) / num_funcs / 10;
    double x = 0;
    for (size_t i = 0; i < basis.num_funcs*10; i++)
    {
        double tmp=0;
        for (size_t j = 0; j < basis.num_funcs; j++)
        {
            tmp += basis(j, x);
        }
        x += step;
        cout<<" "<<tmp<<endl;
    }
    cout<<"\nAreas ";
    for (size_t i = 0; i < basis.num_funcs; i++)
    {
        cout<<basis.areas[i]<<" ";
    }
    cout<<"\navg_e ";
    for (size_t i = 0; i < basis.num_funcs; i++)
    {
        cout<<basis.avg_e[i]<<" ";
    }
    
    double de = max/(num_funcs*20);
    double e=0;
    for (size_t i=0; i<num_funcs*24; i++) {
        fout<<e<<" ";
        dfout<<e<<" ";
        for (size_t j=0; j<num_funcs; j++) {
            fout<<basis(j, e)<<" ";
            dfout<<basis.D(j, e)<<" ";
        }
        fout<<endl;
        dfout<<endl;
        e += de;
    }
    fout.close();
    dfout.close();
    return 0;
}
