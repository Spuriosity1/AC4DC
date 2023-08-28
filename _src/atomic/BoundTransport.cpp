#include "BoundTransport.h"
#include "DecayRates.h"
#include "NumericalIntegral.h"
#include <wignerSymbols.h>
#include <fstream>
#include "ComputeRateParam.h"



BoundTransport::BoundTransport(vector<RadialWF> &light_orbitals, vector<RadialWF>&heavy_orbitals, Input & Input){

    // // Find the lowest energy state for configuration combinations.
	// for (int i = 0; i < light_orbitals.size(); i++){
    //     for (int j = 0; j < heavy_orbitals.size(); j++){

    //     }
    // }

    // for (size_t i = 0;i < dimension - 1; i++)//last configuration is lowest electron count state//dimension-1
    // {
    //     vector<RadialWF> Orbitals = orbitals;
    //     for (size_t j = 0;j < Orbitals.size(); j++) {
    //         Orbitals[j].set_occupancy(orbitals[j].occupancy() - ComputeRateParam::Index[i][j]);
    //         N_elec += Orbitals[j].occupancy();
    //         // Store shell flag if orbital corresponds to shell
    //         if(shell_check[j]) Orbitals[j].flag_shell(true);
    //     }
    // }

	// {
	// 	Orbitals[i].set_N(orbitals[i].N());
	// 	Orbitals[i].set_L(orbitals[i].L());        // Resets occupancy
	// 	Orbitals[i].Energy = orbitals[i].Energy;
	// 	if (orbitals[i].occupancy() != 0)
	// 	{
	// 		N_elec += orbitals[i].occupancy();
	// 		if (j < orbitals[i].pract_infinity()) { j = orbitals[i].pract_infinity(); Infinity = lattice.R(j); }
	// 		if (omega + orbitals[i].Energy >= 0.25)
	// 		{
	// 			if (0.5*k_min*k_min > omega + orbitals[i].Energy) k_min = sqrt(2 * (omega + orbitals[i].Energy));
	// 			if (0.5*k_max*k_max < omega + orbitals[i].Energy) k_max = sqrt(2 * (omega + orbitals[i].Energy));
	// 		}
	// 	}
	// 	Orbitals[i].set_occupancy(orbitals[i].occupancy());
	// }  

// light_orbitals[i].energy;

}
