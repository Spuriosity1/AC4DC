/**
 * @file DynamicRegions.h
 * @brief Defines the classes that handle the dynamic grid mechanics
 */

#include "Constant.h"
#include "GridSpacing.hpp"
#include <vector>
#include <iostream>

using namespace std;

class Region
{
public:
    // // Dummy region
    // Region(): type{"blank"},power{1}{};
    /**
     * @brief Construct a new Region object
     * 
     * @param num_points 
     * @param E_min in eV!
     * @param E_max in eV!
     * @param type one of: "static", "dirac", "mb"
     */
    Region(int num_points, double E_min, double E_max, const string type) :
        num_points{num_points},
        E_min{E_min/Constant::eV_per_Ha},
        E_max{E_max/Constant::eV_per_Ha},
        type{type},
        power{1}   
        {
            num_points_left = num_points;
        }  
    double get_point_density(){return (E_max-E_min)/num_points;};
    double get_E_min(){return E_min;}
    double get_E_max(){return E_max;}
    const string get_type(){return type;}
    // (just used for sorting)
    bool operator < (Region& other){
        return E_max  < (other.get_E_min());
    }    
    double get_next_knot(double previous_knot, bool first_non_zero);
    void update_region(double new_centre,double new_min, double new_max);
    void set_num_points(int n){num_points = n;}
    int get_num_points(){return num_points;}

private:
    int num_points;
    int num_points_left;
    double E_min;
    double E_max;
    string type;
    int power;
};

class GridRegions
{
public:
    GridRegions();
    GridRegions(size_t num_static);
    void update_regions(FeatureRegimes rgm);
    std::vector<Region> regions;
    double mb_min_over_kT; // min of MB region = mb_min_over_kT*kT;
    double mb_max_over_kT; // max of MB region = mb_min_over_kT*kT;

protected:
    void initialise_regions(DynamicGridPreset preset);
    void set_static_region_energies(vector<double> energy_boundaries);
    double dynamic_min_inner_knot();
    double dynamic_max_inner_knot();
    

    // GridRegions& operator+=(const GridRegions& other_GR){
    //     regions.insert(regions.end(), other_GR.regions.begin(), other_GR.regions.end());
    //     return *this;
    // }
    // GridRegions& operator=(const GridRegions& other_GR){
    //     regions = other_GR.regions;
    //     return *this;
    // }    

private:
    int NUM_STATIC_REGIONS = 3;
    int NUM_DYNAMIC_REGIONS = 2;
};

