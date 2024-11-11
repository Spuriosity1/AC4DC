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
 * @note Code chooses max density region, so overlapping regions are allowed.
 * @todo 
 * 1. Add a config file for users to provide custom grids. 
 * 2. Add a dynamic region that will automatically provide support for peaks higher than the photoelectron peaks.
 */

GridRegions::GridRegions(){
    int pts_per_dirac = 10; 
    // Defaults for Maxwell-Boltzmann (thermalised electrons) region bounds, but overwritten in some presets.
    mb_min_over_kT =  0.2922; //  90% of electrons above this point
    mb_max_over_kT =  2.3208; //  80% of electrons below this point (lower since not as sharp)
    first_gp_min_E = 4/Constant::eV_per_Ha;
} 
void GridRegions::initialise_regions(DynamicGridPreset& preset){
    first_gp_min_E = first_gp_min_E*Constant::eV_per_Ha;  // Converted to eV for consistency. Convert to Ha at end.
    preset.min_dirac_region_peak_energy = 1500;  // eV. Convert to Ha at end.
    // Initialise regions
    size_t num_dirac_regions = 4;
    char dirac_region = Region::dirac;
    regions = {};
    // Carbon target (a relatively divergent-prone example) good grid for the unstable early times:
    // indices  1|30|55|65|85|105|145|150 
    // energies 4|10|50|200|500|4500|6500|10000    
    int pts_per_dirac = 1;   
    double max_knot_energy = 2.5*max(preset.pulse_omega,preset.electron_source_energy);
    double trans_scaling = max(1.,preset.pulse_omega/6000);
    string preset_name = "None";
    switch (preset.selected){
      case DynamicGridPreset::high_acc:
        preset_name = "High accuracy";
        mb_max_over_kT = 2.3208*4/3; 
        pts_per_dirac = 35; // (minimum) per region, will have less contributed per photopeak when region is wide enough that point density is lower than some overlapping static regions. 
        regions = {
            Region(10,10,50, Region::fixed), // low support (MB is very fine early on, grid does not do well with sudden transitions in knot density.)
            Region(12,50,200,Region::fixed), 
            Region(25,200,600,Region::fixed), // auger
            Region((int)(0.5+ 8*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(35,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail. An energy ceiling at least twice the photopeak ensures all non-negligible grid points will obey charge conservation (Sanders).
            Region(40,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };
        break;
      case DynamicGridPreset::medium_acc:
      preset_name = "Medium accuracy";
      mb_max_over_kT = 2.3208*4/3;  
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, Region::fixed), Region(8,20,50, Region::fixed),  // low support
            Region(10,50,200,Region::fixed), 
            Region(15,200,600,Region::fixed), // auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(25,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };        
        break;
      case DynamicGridPreset::low_acc:
        preset_name = "Low accuracy";
        pts_per_dirac = 10;
        regions = {
            Region(5,10,20, Region::fixed), Region(5,20,50, Region::fixed),  // low support
            Region(7,50,200,Region::fixed), 
            Region(7,200,600,Region::fixed), // auger
            Region((int)(0.5+ 6*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(15,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(4,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(15,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };    
        break;  
      case DynamicGridPreset::heavy_support:
        preset_name = "Heavy atom support";  // TODO change name this was just the K-shell photoelectrons in Al.
        mb_max_over_kT = 2.3208*4/3;  
        pts_per_dirac = 10;
        regions = {
            //Region(1,1.2,5,Region::fixed),
            Region(7,10,20, Region::fixed), Region(15,20,50, Region::fixed),  // low support
            Region(25,50,200,Region::fixed),  // Slow photoelectron Support 
            Region(20,200,600,Region::fixed), 
            Region((int)(0.5+ 6*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(15,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(4,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(20,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };    
        break;          
      case DynamicGridPreset::Zr_support:
        preset_name = "Zr_fast Auger support (1/2 L hole(s), M->L|M->free (1st order): 1589/1615 eV, M->L|N->free (2nd order): 1918/1970 eV)";  // TODO change name this was just the K-shell photoelectrons in Al.
        mb_max_over_kT = 2.3208*4/3;  
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, Region::fixed), Region(8,20,50, Region::fixed),  // low support
            Region(10,50,200,Region::fixed), 
            Region(10,200,600,Region::fixed), // light auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition (Overlaps with Zr auger support but that's fine since code chooses max density region)
            Region(20,1200,2400,Region::fixed), // Zr Auger support
            Region(4,800,1200,Region::fixed), // Zr Auger lower support
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(25,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };   
        break;     
      case DynamicGridPreset::lower_dirac_support:
        preset_name = "Lower dirac regions";  // Support extends below bottom of transition region.
        preset.min_dirac_region_peak_energy  = 350;
        mb_max_over_kT = 2.3208*4/3;  
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, Region::fixed), Region(8,20,50, Region::fixed),  // low support
            Region(10,50,200,Region::fixed), 
            Region(15,200,600,Region::fixed), // auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(25,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };        
        break;      
      case DynamicGridPreset::Galli_support:
        first_gp_min_E = 2; // Lower min E makes sim. much slower. Unknown if the solution converges better (let alone by how much) when adding knots for these low energies. 
        preset_name = "Galli support";  
        preset.min_dirac_region_peak_energy  = 350;
        mb_max_over_kT = 2.3208*4/3;  
        mb_min_over_kT =  0.2922*1/4;
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, Region::fixed), Region(8,20,50, Region::fixed),  // low support
            Region(10,50,200,Region::fixed), 
            Region(15,200,600,Region::fixed), // auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(40,-1,-1,Region::mb_log), // Maxwell-boltzmann distribution
        };        
        break;      
      case DynamicGridPreset::all_log_grid:  // Similar to preset "dirac", with lower min energy like "Galli_support". However regions are all logarithmic.
        preset_name = "Lower dirac regions, log all";  // Support extends below bottom of transition region.
        preset.min_dirac_region_peak_energy  = 350;
        mb_max_over_kT = 2.3208*4/3;  
        mb_min_over_kT = 0.2922*1/4; 
        pts_per_dirac = 15;  
        ////
        dirac_region = Region::dirac_log;
        regions = {
            Region(7,10,20, Region::fixed_log), Region(8,20,50, Region::fixed_log),  // low support
            Region(10,50,200,Region::fixed_log), 
            Region(15,200,600,Region::fixed_log), // auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,Region::fixed_log), // transition
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed_log),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed_log), // high tail
            Region(25,-1,-1,Region::mb_log), // Maxwell-boltzmann distribution
        };        
        break;       
      case DynamicGridPreset::mb_log_grid:  // Similar to preset "dirac", with lower min energy like "Galli_support". However MB region is logarithmic.
        preset_name = "Lower dirac regions, log MB";  // Support extends below bottom of transition region.
        preset.min_dirac_region_peak_energy  = 350;
        mb_max_over_kT = 2.3208*4/3;  
        mb_min_over_kT = 0.2922*1/4; 
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, Region::fixed), Region(8,20,50, Region::fixed),  // low support
            Region(10,50,200,Region::fixed), 
            Region(15,200,600,Region::fixed), // auger
            Region((int)(0.5+ 7*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(25,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(7,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(25,-1,-1,Region::mb_log), // Maxwell-boltzmann distribution
        };        
        break;      
      case DynamicGridPreset::static_high_energy:  // i.e. Static dirac. Same as all_log_grid, except one big logarithmic static region past 1/4 photon energy.
        preset_name = "Static high energy regions, log all";  // Support extends below bottom of transition region.
        preset.min_dirac_region_peak_energy  = 350;
        mb_max_over_kT = 2.3208*4/3;  
        mb_min_over_kT = 0.2922*1/4; 
        pts_per_dirac = 15;  
        regions = {
            Region(7,10,20, Region::fixed_log), Region(8,20,50, Region::fixed_log),  // low support
            Region(10,50,200,Region::fixed_log), 
            Region(15,200,600,Region::fixed_log), // auger
            Region(30,600,preset.pulse_omega/3,Region::fixed_log), // transition
            Region(50,preset.pulse_omega/3,max_knot_energy/2.5*6/4,Region::fixed_log),  // high-energy
            Region(10,max_knot_energy/2.5*6/4,max_knot_energy,Region::fixed_log),  // tail
            Region(25,-1,-1,Region::mb_log), // Maxwell-boltzmann distribution
        };        
        break;                       
      case DynamicGridPreset::dismal_acc:
        preset_name = "Dismal accuracy";
        pts_per_dirac = 5;
        regions = {
            Region(4,10,20, Region::fixed), Region(4,20,50, Region::fixed),  // low support
            Region(5,50,200,Region::fixed), 
            Region(5,200,600,Region::fixed), // auger
            Region((int)(0.5+ 3*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(10,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(4,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(10,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };
        break;
      // Effectively replaces dynamic dirac regions with static 30 points around dirac region
      case DynamicGridPreset::no_dirac:
        preset_name = "No dynamic dirac";
        pts_per_dirac = 1;
        regions = {
            Region(15,10,50, Region::fixed), // low support
            Region(8,50,200,Region::fixed), 
            Region(8,200,600,Region::fixed), // auger
            Region(30+(int)(0.5+ 4*trans_scaling),600,preset.pulse_omega*1.1,Region::fixed),  
            Region(10,preset.pulse_omega*1.1,max_knot_energy,Region::fixed), // high tail
            // Region(4,600,preset.pulse_omega/4,Region::fixed), // transition
            // Region(20,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            // Region(5,preset.pulse_omega*6/4,preset.pulse_omega*2,Region::fixed), // high tail            
            Region(25,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };     
        break;     
      // Idea: optimal grid dynamics should look something like noticeably more grid points early than late, as runge's phenomenon is most significant early on.
      // setting dynamic grid point count such that after a few fs their density is comparable to overlapping static region densities forces this.
      case DynamicGridPreset::training_wheels:
        preset_name = "Training wheels";
        pts_per_dirac = 15;
        regions = {
            Region(10,10,50, Region::fixed), // low support
            Region(8,50,200,Region::fixed), 
            Region(8,200,600,Region::fixed), // auger
            Region((int)(0.5+ 5*trans_scaling),600,preset.pulse_omega/4,Region::fixed), // transition
            Region(30,preset.pulse_omega/4,preset.pulse_omega*6/4,Region::fixed),  // photo
            Region(5,preset.pulse_omega*6/4,max_knot_energy,Region::fixed), // high tail
            Region(15,-1,-1,Region::mb), // Maxwell-boltzmann distribution
        };        
        break;              
    }
    // Region for electron source.
    if(preset.electron_source_energy != -1){
        char static_region = Region::fixed;
        if (preset.selected == DynamicGridPreset::all_log_grid || preset.selected == DynamicGridPreset::static_high_energy)
            static_region = Region::fixed_log;
        // Inner region with high density of knots, outer region of support to ensure not too high a difference with surroundings.
        regions.push_back(Region(int(pts_per_dirac),preset.electron_source_energy*3/4,preset.electron_source_energy*4/3,static_region));
        regions.push_back(Region(int(pts_per_dirac),preset.electron_source_energy*7/8,preset.electron_source_energy*8/7,static_region));
        regions.push_back(Region(int(pts_per_dirac),preset.electron_source_energy*1/2,preset.electron_source_energy*2,static_region));
    }    
    std::vector<Region> common_regions = {
        Region(1,5,10,Region::fixed),  // low divergent (purpose is just placing a point at 5 eV. Below this is unnecessarily costly, and sometimes breaks - either because I have brittle code or it's fundamentally untenable.)
        // 4 Photoelectron peaks. need to include num regions defined here in FeatureRegimes. TODO fix this issue.
    };     
    for(size_t i = 0; i < num_dirac_regions; i++){
                common_regions.push_back(Region(pts_per_dirac,-1,-1,dirac_region));
    }
    regions.insert(regions.end(), common_regions.begin(), common_regions.end() );
    std::cout <<"\033[38:5:208m"<<"[ Dynamic Grid ] '"<< preset_name << "' preset used to initialise dynamic regions."<<"\033[0m"<< endl;
    
    preset.name = preset_name;
    // Convert to atomic units
    preset.min_dirac_region_peak_energy  /= Constant::eV_per_Ha;
    first_gp_min_E /= Constant::eV_per_Ha;
}


/**
 * @brief Update dynamic regions' energy bounds. (Using linear regions).
 * @details currently places centre of region on peak.
 */
void GridRegions::update_regions(FeatureRegimes rgm){
    std::cout << "[ Dynamic Grid ] Updating regions" << std::endl;
    size_t peak_idx = 0;
    for (size_t r = 0; r < regions.size(); r ++){
        switch(regions[r].get_type()){
            case Region::fixed:
            case Region::fixed_log:
            break; 
            case Region::dirac:
            case Region::dirac_log:
                regions[r].update_region(rgm.dirac_peaks[peak_idx],rgm.dirac_minimums[peak_idx],rgm.dirac_maximums[peak_idx]);
                peak_idx++;
            break;
            case Region::mb:
            case Region::mb_log:
                regions[r].update_region(rgm.mb_peak,rgm.mb_min,rgm.mb_max);
            break;
            default:
                std::cout <<"Error, unrecognised region type: " <<regions[r].get_type() << std::endl;
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
    assert(type != Region::fixed);    
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
    double next_knot = previous_knot + get_inv_point_density(previous_knot);
    if ((E_min <= next_knot && next_knot <= E_max)){
        return next_knot;
    }
    // ... Or if the next knot would be below the region, the lowest knot that is the first point of this region instead.
    if(E_min > previous_knot) return E_min;

    // Otherwise the next knot is outside this region!
    return -1;
}
double Region::get_inv_point_density(double previous_knot){
    switch(get_type()){
        case Region::fixed:
        case Region::dirac:
        case Region::mb:
        {
            return (E_max-E_min)/num_points;
        }
        break;
        case Region::fixed_log:
        case Region::dirac_log:            
        case Region::mb_log:
        {   
            double a = get_E_min(); double  b = get_E_max(); 
            double x = previous_knot; double N = get_num_points();


            double last_point = a;
            double next_point = a;
            double last_size = -1;           
            // Don't take derivative, construct dummy knot to get spacing between discretised grid points.
            for(size_t n=2; n<= N; n++) {
                 last_point = next_point; 
                 next_point = a*pow(b/a,(n-1)/(N-1));
                 if(next_point > x)
                    return next_point-last_point;
            }                        
            // double p = 1.05;
            // if (x < 0)
            //     return INFINITY;
            // assert(N > 1); // -> x > a - A
            
            // double x_test = x*Constant::eV_per_Ha;
            // double K = (pow(p,N-1)*a-b)/(1-pow(p,N-1));
            // double A = p/(a+K);
            
            // // Construct dummy knot to get spacing between discretised grid points.
            // std::vector<double> discrete_points(N);
            // double last_point = a;
            // double next_point = a;
            // double last_size = -1;
            // for(size_t n=2; n<= N; n++) {
            //      last_point = next_point;
            //      next_point = pow(p,n)/A-K;
            //      double test_point = Constant::eV_per_Ha*last_point;
            //      double test_point2 = Constant::eV_per_Ha*next_point;
            //      double test_x = Constant::eV_per_Ha*x;
            //      if(next_point > x)
            //         return next_point-last_point;
            // }
        }
        break;
        default:
            std::cout <<"Error, unrecognised region type: " <<get_type()<< std::endl;
    }    
    return -1;
}

