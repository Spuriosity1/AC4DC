#include "src/Dipole.h"
#include "src/Constant.h"
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    if(argc < 4) {
        cerr<<"Usage: "<<argv[0]<<" T B N "<<endl;
        return 1;
    }
    double T = atof(argv[1]);
    double B = atof(argv[2]);
    int num_integral = atoi(argv[3]);
    double U = B/7;
    double occ = 2;

    double tot_sigma = Dipole::sigmaBEB(T, B, U, occ);

    cout<<"Total cross-section = "<<tot_sigma<<endl;
     
    double max_e = T-B;
    // double min_e = 0;
    double tmp =0;
    // Naive integration
    double de = max_e/num_integral;
    for (size_t i = 0; i < num_integral; i++)
    {
        double W = i*de;
        tmp += Dipole::DsigmaBEB(T, W, B, U, occ);
    }
    tmp *= de;
    cout<<"Integrated differential cross-section over [0, T-B] = "<<tmp<<endl;
    double err = 2*tot_sigma - tmp;
    cout<<"2sigma - ∫dsigma = "<<err<<endl;
    cout<<"Discrepancy: "<<100*err/tot_sigma<<"%"<<endl;
    cerr<<"Testing symmetry:"<<endl;
    cerr<<"W, sigma(W), sigma(T-B-W)"<<endl;
    int num_points = 30;
    de = T/num_points/2;
    double W1 = 0;
    double W2 = T - B;
    for(size_t i=0; i<num_points; i++) {
        W1 += de;
        W2 -= de;
        cerr<<W1 <<" "<< Dipole::DsigmaBEB(T, W1, B, U, occ)<<" "<<Dipole::DsigmaBEB(T, W2, B, U, occ)<<endl;
    }
    
    
    return 0;
}
