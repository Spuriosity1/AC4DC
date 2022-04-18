#include "src/SplineBasis.hpp"
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

int main(int argc, char const *argv[]) {
    if (argc < 5) {
        std::cout<<"Usage: "<<argv[0]<<" [fstem] [num_funcs] [min_e (Ha)] [max_e (Ha)] (linear|quadratic|exponential) "<<"\n";
    }
    std::string dummy(argv[1]);
    std::ofstream fout(dummy+"_vals.csv");
    std::ofstream dfout(dummy+"_derivative.csv");
    int num_funcs = atoi(argv[2]);
    double min = atof(argv[3]);
    double max = atof(argv[4]);
    istd::stringstream is(argv[5]);
    GridSpacing gs;
    is >> gs;
    BasisSet basis;
    gs.zero_degree_0 = 1;
    gs.zero_degree_inf = 0;
    basis.set_parameters(num_funcs, min, max, gs);
    std::cout<<"Testing partition of unity"<<"\n";
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
        std::cout<<" "<<tmp<<"\n";
    }
    std::cout<<"\nAreas ";
    for (size_t i = 0; i < basis.num_funcs; i++)
    {
        std::cout<<basis.areas[i]<<" ";
    }
    std::cout<<"\navg_e ";
    for (size_t i = 0; i < basis.num_funcs; i++)
    {
        std::cout<<basis.avg_e[i]<<" ";
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
        fout<<"\n";
        dfout<<"\n";
        e += de;
    }
    fout.close();
    dfout.close();
    return 0;
}
