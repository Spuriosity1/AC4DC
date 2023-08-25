#include "src/Constant.h"
#include <vector>
#include <iostream>

using namespace std;

int main(int argc, char const *argv[])
{
    if (argc < 3){
        cerr<<"Usage: "<<argv[0]<<" [infile] [outfile]\n";
        return 1;
    }
    std::vector<RateData::EIIdata> eiiVec;
    RateData::ReadEIIParams(argv[1], eiiVec);
    RateData::WriteEIIParams(argv[2], RateData::inverse(eiiVec));
    
    return 0;
}
