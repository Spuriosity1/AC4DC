#include "src/ImpactCrossSection.hpp"
#include "src/Constant.hpp"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    if(argc < 4) {
        std::cerr<<"Usage: "<<argv[0]<<" T B N "<<"\n";
        return 1;
    }
    double T = atof(argv[1]);
    double B = atof(argv[2]);
    int num_integral = atoi(argv[3]);
    double U = B/7;
    double occ = 2;

    double tot_sigma = ImpactCrossSection::sigmaBEB(T, B, U, occ);

    std::cout<<"Total cross-section = "<<tot_sigma<<"\n";
     
    double max_e = T-B;
    // double min_e = 0;
    double tmp =0;
    // Naive integration
    double de = max_e/num_integral;
    for (size_t i = 0; i < num_integral; i++)
    {
        double W = i*de;
        tmp += ImpactCrossSection::DsigmaBEB(T, W, B, U, occ);
    }
    tmp *= de;
    std::cout<<"Integrated differential cross-section over [0, T-B] = "<<tmp<<"\n";
    double err = 2*tot_sigma - tmp;
    std::cout<<"2sigma - ∫dsigma = "<<err<<"\n";
    std::cout<<"Discrepancy: "<<100*err/tot_sigma<<"%"<<"\n";
    std::cerr<<"Testing symmetry:"<<"\n";
    std::cerr<<"W, sigma(W), sigma(T-B-W)"<<"\n";
    int num_points = 30;
    de = T/num_points/2;
    double W1 = 0;
    double W2 = T - B;
    for(size_t i=0; i<num_points; i++) {
        W1 += de;
        W2 -= de;
        std::cerr<<W1 <<" "<< ImpactCrossSection::DsigmaBEB(T, W1, B, U, occ)<<" "<<ImpactCrossSection::DsigmaBEB(T, W2, B, U, occ)<<"\n";
    }
    
    
    return 0;
}
