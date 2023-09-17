#pragma once

#include "RadialWF.h"
#include "Grid.h"
#include "Potential.h"
#include <vector>
#include "Constant.h"
#include <fstream>
#include "HFInputParam.h"
#include <RateData.h>

using namespace RateData;

class BoundTransport{
public:
    /**
    * @brief // Finds the lowest energy configuration allowed 
    * via transport of valence electrons to the next unoccupied orbital that is at or above the outer orbital.
    * @param light_orbitals 
    * @param heavy_orbitals 
    * @param Input 
    * @return 
    */
    BoundTransport(vector<RadialWF> &light_orbitals, vector<RadialWF>&heavy_orbitals, InputData::HFInputParam & Input);
    
    bound_transport bbe_transport;

};
