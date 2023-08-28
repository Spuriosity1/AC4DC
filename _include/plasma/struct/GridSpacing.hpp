/**
 * @file GridSpacing.hpp
 * 
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

#include <vector>
#include <iostream>
#include <algorithm>
#ifndef GRIDSPACING_CXX_H
#define GRIDSPACING_CXX_H

// Determines GridRegions.regions
struct DynamicGridPreset{
    const static char dismal_acc = 0;
    const static char low_acc = 1;
    const static char medium_acc = 2;
    const static char high_acc = 3;
    const static char no_dirac = 4;
    const static char training_wheels = 5;
    const static char heavy_support = 6;
    const static char unknown = 101;
    char selected = unknown;  
    double pulse_omega = -1;  // Photon energy [eV]
};

struct GridSpacing {
    const static char manual = 0;
    const static char dynamic = 1;
    const static char unknown = 101;
    char mode = dynamic;  // default is dynamic, since manual grid true/false is contained in manual grid header, which we want to be able to exclude from input file.
    unsigned zero_degree_0 = 0; // Number of derivatives to set to zero at the origin
    unsigned zero_degree_inf = 0; // Number of derivatives to set to zero at infinity
};

struct Cutoffs{
    double transition_e = 0; //used to estimate Coulomb logarithm. 
    double min_coulomb_density=0; // Density below which to ignore Coulomb interactions 
};

struct FeatureRegimes{
    double mb_peak=0, mb_min=0, mb_max=0;  
    size_t num_dirac_peaks = 4;
    std::vector<double> dirac_peaks = std::vector<double>(num_dirac_peaks); 
    std::vector<double> dirac_minimums = std::vector<double>(num_dirac_peaks); 
    std::vector<double> dirac_maximums = std::vector<double>(num_dirac_peaks); 
};

struct ManualGridBoundaries {
    // Note this current implementation includes the energy/index of last grid point.
    std::vector<int> bndry_idx ={-1};  
    std::vector<double> bndry_E ={-1.};  // [eV]
    std::vector<double> powers ={1.};
    int parsed_count = 0; 
};

namespace {
    [[maybe_unused]] std::ostream& operator<<(std::ostream& os, GridSpacing gs) {
        switch (gs.mode)
        {
        case GridSpacing::manual:
            os << "manual";
            break;
        case GridSpacing::dynamic:
            os << "dynamic";
            break;
        default:
            os << "Unknown grid type";
            break;
        }
        return os;
    }

    /// Sets grid style/type to dynamic or manual
    [[maybe_unused]] std::istream& operator>>(std::istream& is, GridSpacing& gs) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No grid mode provided, defaulting to dynamic mode..."<<std::endl;
            gs.mode = GridSpacing::dynamic;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'T':
            gs.mode = GridSpacing::manual;
            break;
        case 't':
            gs.mode = GridSpacing::manual;
            break;
        case 'F':
            gs.mode = GridSpacing::dynamic;
            break;
        case 'f':
            gs.mode = GridSpacing::dynamic;
            break;          
        default:
            std::cerr<<"Unrecognised grid mode \""<<tmp<<"\", defaulting to dynamic mode..."<<std::endl;
            gs.mode = GridSpacing::dynamic;
            break;
        }
        return is;
    }
    /// Sets dynamic grid regions' preset. TODO replace input with full string.
    [[maybe_unused]] std::istream& operator>>(std::istream& is, DynamicGridPreset& preset) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No dynamic grid preset provided, defaulting to medium accuracy..."<<std::endl;
            preset.selected = DynamicGridPreset::medium_acc;
            return is;
        }
        switch ((char) tmp[0])
        {
        case 'd':
            preset.selected = DynamicGridPreset::dismal_acc;
            break;
        case 'l':
            preset.selected = DynamicGridPreset::low_acc;
            break;
        case 'm':
            preset.selected = DynamicGridPreset::medium_acc;
            break;
        case 'h':
            preset.selected = DynamicGridPreset::high_acc;
            break;    
        case 'n':
            preset.selected = DynamicGridPreset::no_dirac;
            break;         
        case 't':
            preset.selected = DynamicGridPreset::training_wheels;
            break;                      
        case 'A':
            preset.selected = DynamicGridPreset::heavy_support;
            break;                      
        default:
            std::cerr<<"Unrecognised grid preset \""<<tmp<<"\", defaulting to medium accuracy..."<<std::endl;
            preset.selected = DynamicGridPreset::medium_acc;
            break;
        }
        return is;
    }    
    /// Sets region boundaries
    [[maybe_unused]] std::istream& operator>>(std::istream& is, ManualGridBoundaries& gb) {
        std::string tmp;
        is >> tmp;
        if (tmp.length() == 0) {
            std::cerr<<"No boundaries provided, default case not implemented yet...!"<<std::endl;
            return is;
        }

        // Get vector from string. I'm not proud of this.
        std::vector<double> target_vector;
        std::string delimiter_char = "|";
        std::string skip_char = " ";
        size_t pos = 0;
        //size_t space_pos = 0;
        double token;
        std::string token_str;  
        while ((pos = tmp.find(delimiter_char)) != std::string::npos) {
            token_str = tmp.substr(0, pos);
            // // get rid of spaces 
            // std::string::iterator end_pos = std::remove(token_str.begin(), token_str.end(), ' ');
            // token_str.erase(end_pos, token_str.end());     
            // space_pos = token_str.find(' ');
            // token_str = tmp.substr(0, space_pos);

            // Add num to vector
            token =  stod(token_str);
            target_vector.push_back(token);
            // Erase input string up to next delimiter
            tmp.erase(0, pos + delimiter_char.length());
        }        
        if(tmp.size() > 0){
            std::string::iterator end_pos = std::remove(tmp.begin(), tmp.end(), ' ');
            tmp.erase(end_pos, tmp.end());
            target_vector.push_back(stod(tmp));
        }

        if(gb.parsed_count == 0){
            gb.bndry_idx.resize(target_vector.size());
            for(size_t i = 0; i < target_vector.size(); i++){
                target_vector[i] += 0.5;
                gb.bndry_idx[i] = (int)target_vector[i];
            }
        }
        else if (gb.parsed_count == 1) gb.bndry_E = target_vector;
        else if (gb.parsed_count == 2) gb.powers = target_vector;
        else{
            std::cerr << "Too many grid region variables parsed.";}
        gb.parsed_count++;
        return is;
    }    
}

#endif