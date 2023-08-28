#include "src/ComputeRateParam.h"

// This file reads a rate from file and outputs what it thinks is the same thing.
// Diagnostic tool for parser errors.

using namespace std;
int main(int argc, char const *argv[])
{
    if (argc < 3) {
        cerr<<"Need to specify an input and output file."<<endl;
        cerr<<"Usage: rate_io_test output/C/Xsections/EII.json testOutput/C/EII.json"<<endl;
        cerr<<"       rate_io_test output/C/Xsections/Auger.txt"<<endl;
    }
    string instring(argv[1]);
    string outstring(argv[2]);
    string ext = instring.substr(instring.rfind('.'));
    if (ext == ".json") {
        // EII data
        std::vector<RateData::EIIdata> data;
        RateData::ReadEIIParams(instring, data);
        RateData::WriteEIIParams(outstring, data);
    } else if (ext == ".txt") {
        // rate
        std::vector<RateData::Rate> data;
        RateData::ReadRates(instring, data);
        RateData::WriteRates(outstring, data);
    } else {
        cerr<<"Unrecognised file type"<<endl;
    }
    // RateData::ReadRates(argv[2])
    return 0;
}

