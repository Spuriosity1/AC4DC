/**
 * @file DynamicRegions.cpp
 * @author Spencer Passmore
 * @brief 
 * @note may want to change to have the raw files save with not much more fineness than the load ones.
 */
/*===========================================================================
This file is part of AC4DC.

    AC4DC is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    AC4DC is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with AC4DC.  If not, see <https://www.gnu.org/licenses/>.
===========================================================================*/
// (C) Spencer Passmore 2023

#include "DynamicRegions.h"
#include <iostream>
#include <vector>
#include <variant>
#include <math.h>
#include <execution>


/**
 * @brief 
 * @details  
 * Static regions (named for the feature they are fitting when dynamic region not present):
 * - Low E divergent region
 * - Auger region
 * - Tail region
 * Dynamic regions:
 * - MB region
 * - Photoelectron region
 * transition region has more knots for higher energies (as it becomes quite wide and we already are using a low number of knots).
 * @todo 
 * 1. Add a config file for users to provide custom grids. 
 * 2. Add a dynamic region that will automatically provide support for peaks higher than the photoelectron peaks.
 */

GridRegions::GridRegions(){
    int pts_per_dirac = 10; 
    // Defaults for Maxwell-Boltzmann (thermalised electrons) region bounds, but overwritten in some presets.
    mb_min_over_kT =  0.2922; //  90% of electrons above this point
    mb_max_over_kT =  2.3208; //  80% of electrons below this point (lower since not as sharp)
} 
void GridRegions::initialise_regions(DynamicGridPreset preset){
    // Initialise regions
    regions = {};
    // Carbon target (a relatively divergent-prone example) good grid for the unstable early times:
    // indices  1|30|55|65|85|105|145|150 
    // energies 4|10|50|200|500|4500|6500|10000    
    int pts_per_dirac = 1;   
    double trans_scaling = max(1.,preset.pulse_omega/6000);
    string preset_name = "None";
    switch (preset.selected){
      case DynamicGridPreset::high_acc:
        preset_name = "High accuracy";
        mb_max_over_kT = 2.3208*4/3; 
        pts_per_dirac = 35; // (minimum) per region, will have less contributed per photopeak when region is wide enough that point density is lower than some overlapping static regions. 
        regions = {
            Region(10,10,50, "static"), // low support (MB is very fine early on, grid does not do well with sudden transitions in knot density.)
            Region(12,50,200,"static"), 
            Region(25,200,600,"static"), // auger
            Region((int)(0.5+ 8*trans_scaling),600,preset.pulse_omega/4,"static"), // transition
            Region(35,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            Region(7,preset.pulse_omega*6/4,preset.pulse_omega*2.5,"static"), // high tail. An energy ceiling at least twice the photopeak ensures all non-negligible grid points will obey charge conservation (Sanders).
            Region(40,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };
        break;
      case DynamicGridPreset::medium_acc:
      preset_name = "Medium accuracy";
      mb_max_over_kT = 2.3208*4/3;  
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, "static"), Region(8,20,50, "static"),  // low support
            Region(10,50,200,"static"), 
            Region(15,200,600,"static"), // auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,"static"), // transition
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            Region(7,preset.pulse_omega*6/4,preset.pulse_omega*2.5,"static"), // high tail
            Region(25,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };        
        break;
      case DynamicGridPreset::low_acc:
        preset_name = "Low accuracy";
        pts_per_dirac = 10;
        regions = {
            Region(5,10,20, "static"), Region(5,20,50, "static"),  // low support
            Region(7,50,200,"static"), 
            Region(7,200,600,"static"), // auger
            Region((int)(0.5+ 6*trans_scaling),600,preset.pulse_omega/4,"static"), // transition
            Region(15,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            Region(4,preset.pulse_omega*6/4,preset.pulse_omega*2.5,"static"), // high tail
            Region(15,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };    
        break;  
      case DynamicGridPreset::heavy_support:
        preset_name = "Heavy atom support";  // TODO change name this was just the K-shell photoelectrons in Al.
        mb_max_over_kT = 2.3208*4/3;  
        pts_per_dirac = 10;
        regions = {
            //Region(1,1.2,5,"static"),
            Region(5,10,20, "static"), Region(5,20,50, "static"),  // low support
            Region(12,50,200,"static"),  // auger
            Region(12,200,600,"static"), 
            Region((int)(0.5+ 6*trans_scaling),600,preset.pulse_omega/4,"static"), // transition
            Region(15,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            Region(4,preset.pulse_omega*6/4,preset.pulse_omega*2.5,"static"), // high tail
            Region(20,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };    
        break;          
      case DynamicGridPreset::dismal_acc:
        preset_name = "Dismal accuracy";
        pts_per_dirac = 5;
        regions = {
            Region(4,10,20, "static"), Region(4,20,50, "static"),  // low support
            Region(5,50,200,"static"), 
            Region(5,200,600,"static"), // auger
            Region((int)(0.5+ 3*trans_scaling),600,preset.pulse_omega/4,"static"), // transition
            Region(10,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            Region(4,preset.pulse_omega*6/4,preset.pulse_omega*2.5,"static"), // high tail
            Region(10,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };
        break;
      // Effectively replaces dynamic dirac regions with static 30 points around dirac region
      case DynamicGridPreset::no_dirac:
        preset_name = "No dynamic dirac";
        pts_per_dirac = 1;
        regions = {
            Region(15,10,50, "static"), // low support
            Region(8,50,200,"static"), 
            Region(8,200,600,"static"), // auger
            Region(30+(int)(0.5+ 4*trans_scaling),600,preset.pulse_omega*1.1,"static"),  
            Region(10,preset.pulse_omega*1.1,preset.pulse_omega*2.5,"static"), // high tail
            // Region(4,600,preset.pulse_omega/4,"static"), // transition
            // Region(20,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            // Region(5,preset.pulse_omega*6/4,preset.pulse_omega*2,"static"), // high tail            
            Region(25,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };     
        break;     
      // Idea: optimal grid dynamics should look something like noticeably more grid points early than late, as runge's phenomenon is most significant early on.
      // setting dynamic grid point count such that after a few fs their density is comparable to overlapping static region densities forces this.
      case DynamicGridPreset::training_wheels:
        preset_name = "No dynamic dirac";
        pts_per_dirac = 15;
        regions = {
            Region(10,10,50, "static"), // low support
            Region(8,50,200,"static"), 
            Region(8,200,600,"static"), // auger
            Region((int)(0.5+ 5*trans_scaling),600,preset.pulse_omega/4,"static"), // transition
            Region(30,preset.pulse_omega/4,preset.pulse_omega*6/4,"static"),  // photo
            Region(5,preset.pulse_omega*6/4,preset.pulse_omega*2.5,"static"), // high tail
            Region(15,-1,-1,"mb"), // Maxwell-boltzmann distribution
        };        
        break;              
    }
    std::vector<Region> common_regions = {
        Region(1,5,10,"static"),  // low divergent (purpose is just placing a point at 5 eV. Below this is unnecessarily costly, and sometimes breaks - either because I have brittle code or it's fundamentally untenable.)
        // 4 Photoelectron peaks. need to include num regions defined here in FeatureRegimes. TODO fix this issue.
        Region(pts_per_dirac,-1,-1,"dirac"), Region(pts_per_dirac,-1,-1,"dirac"), 
        Region(pts_per_dirac,-1,-1,"dirac"), Region(pts_per_dirac,-1,-1,"dirac")
    };     
    regions.insert(regions.end(), common_regions.begin(), common_regions.end() );
    std::cout << "[ Dynamic Grid ] '"<< preset_name << "' preset used to initialise dynamic regions." << endl;
}


/**
 * @brief Update dynamic regions' energy bounds. (Using linear regions).
 * @details currently places centre of region on peak.
 */
void GridRegions::update_regions(FeatureRegimes rgm){
    std::cout << "[Dynamic Grid] Updating regions " << std::endl;
    size_t peak_idx = 0;
    for (size_t r = 0; r < regions.size(); r ++){
        switch(regions[r].get_type()[0]){
            case 's': // static
            break; 
            case 'd':
                regions[r].update_region(rgm.dirac_peaks[peak_idx],rgm.dirac_minimums[peak_idx],rgm.dirac_maximums[peak_idx]);
                peak_idx++;
            break;
            case 'm':
                regions[r].update_region(rgm.mb_peak,rgm.mb_min,rgm.mb_max);
            break;
            default:
                std::cout <<"Error, unrecognised region type" <<regions[r].get_type() << std::endl;
        }
    }
    assert(rgm.dirac_peaks.size() == peak_idx);
}


double GridRegions::dynamic_min_inner_knot(){
    double min_E = INFINITY;
    for(Region reg : regions){
        min_E = min(reg.get_E_min(),min_E);
    }
    return min_E;
 
}
double GridRegions::dynamic_max_inner_knot(){
    double max_E = -INFINITY;
    for(Region reg : regions){
        max_E = max(reg.get_E_max(),max_E);
    }
    return max_E;
}
// void GridRegions::set_static_region_energies(vector<double> energy_boundaries){

// }

void Region::update_region(double new_centre, double new_min, double new_max){
    assert(type != "static");    
    if (E_min != new_min || E_max != new_max){
        std::cout <<"energy range of region of type '"<<type<<"' updated from:\n" 
        <<E_min*Constant::eV_per_Ha<<" - "<<E_max*Constant::eV_per_Ha << std::endl;    
        E_min = new_min;
        E_max = new_max;
        std::cout << "to:\n" 
        <<E_min*Constant::eV_per_Ha<<" - "<<E_max*Constant::eV_per_Ha << std::endl;
    }
    else std::cout<<"Energy range of region of type '"<<type<<"' had NO update." <<std::endl;
}

// powers not implemented yet since doesn't seem necessary
/**
 * @brief Returns next knot, or -1 if it's outside the region's (energy) range. 
 * 
 * @param previous_knot 
 * @return double 
 */
double Region::get_next_knot(double previous_knot, bool first_non_zero){
    if (first_non_zero){
        return E_min;
    }
    if(power != 1){std::cout << "WARNING powers not implemented for dynamic grid" <<std::endl;}
    // We get the knot that is lowest relative to the previous knot + the knot density and above E_min
    double next_knot = previous_knot + get_point_density();
    if ((E_min <= next_knot && next_knot <= E_max)){
        return next_knot;
    }
    // ... Or if the next knot would be below the region, the lowest knot that is the first point of this region instead.
    if(E_min > previous_knot) return E_min;
    return -1;
}

