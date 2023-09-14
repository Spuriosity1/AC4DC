#include "RateData.h"
#include <math.h>
#include <iostream>
#include <fstream>
#include <random>

using namespace RateData;
using namespace std;

void rand_rate(Rate& rate){
    rate.energy = (double) rand()/RAND_MAX;
    rate.from = rand() % 100;
    rate.to = rand() % 100;
    rate.val = 1 + (double) rand()/RAND_MAX;
}

void rand_eii(EIIdata& eii){
    eii.init = rand();
    eii.resize(0);
    for (unsigned i=0; i< rand() % 20; i++){
        eii.push_back(rand()%100,rand()%100,(double) rand()/RAND_MAX, 10 + (double) rand()/RAND_MAX);
    }
}


int main(int argc, char const *argv[])
{
    if (argc > 1){
        srand(atoi(argv[1]));
    } else {
        srand(0);
    }

    std::vector<Rate> rates(100);
    for (auto& rate : rates){
        rand_rate(rate);
    }

    json j1 = rates[0];

    cout << "Serialising Rate ...\n";
    cout << std::scientific << j1;
    cout<<"\n\n";

    std::vector<EIIdata> eiidata(100);
    for (auto& eii : eiidata){
        rand_eii(eii);
    }

    json j2 = eiidata[0];
    cout << "Serialising EII ...\n";
    cout << j2;
    cout<<"\n\n";

    Atom a;
    a.atomic_density = 1.9;
    a.Auger.resize(100);
    a.EIIparams.resize(100);
    a.Fluor.resize(100);
    a.Photo.resize(100);


    return 0;
}
