/**
 * @file SpineBasis.cpp
 * @brief @copybrief SplineBasis.h
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

#include "SplineBasis.h"
#include <assert.h>
#include "Constant.h"
#include "BSpline.h"
#include <algorithm>

// classic false position method
double find_root(std::function<double(double)> func, double a, double b, double tolerance=1e-8){
    assert(b > a);
    assert(tolerance > 0);
    double x = 0;
    while (fabs(func(x)) > tolerance){
        // x = (a*func(b) - b*func(a))/(func(b) - func(a));
        x = (a+b)/2;
        // split interval into [a, x]u[x, b]
        double fx = func(x);
        if (fx*func(a) < 0){
            // root in interval [x, a]
            b = x;
        } else if (fx*func(b) < 0){
            a = x;
        } else {
            throw runtime_error("Root Find Error: f does not change sign");
        }
    }
    return x;
}


//                  value to find, sorted array, length of array, index that *arr refers to
// Returns an index idx such that val lies in the interval [arr[idx], arr[idx+1])
int r_binsearch(double val, const double* arr, int len, int idx) {
    if (len==1) return idx;
    int mid = len/2;
    // arr has indices i i+1 ... i+len-1
    if (*(arr + mid) > val) {
        return r_binsearch(val, arr, mid, idx);
    } else {
        return r_binsearch(val, arr+mid, len - mid, idx+mid);
    }
}

/// @todo make variable names clearer
int lower_idx(double val, const std::vector<double>& arr){
    return r_binsearch(val, arr.data(), arr.size(), 0);
}

/// @todo make variable names clearer
int BasisSet::i_from_e(double e) {
    // size_t idx = r_binsearch(e, avg_e.data(), avg_e.size(), 0);
    size_t idx = lower_idx(e, avg_e);
    if (idx == avg_e.size()) return 0;
    // returns the closer one
    return (avg_e[idx+1] - e > e - avg_e[idx]) ? idx : idx +1;
}

int BasisSet::lower_i_from_e(double e) {
    // size_t idx = r_binsearch(e, avg_e.data(), avg_e.size(), 0);
    return lower_idx(e, avg_e);
}

/**
 * @brief Constructs the grid over which the electron distribution is solved. 
 * @details Idea: Want num_funcs usable B-splines
 *  If grid has num_funcs+k+1 points, num_funcs of these are usable splines
 *  grid layout for open boundary:
 *  t0=t1=t2=...=tk-1,   tn+1 = tn+2 =... =tn+k
 *  homogeneous Dirichlet boundary:
 *  t0=t1=t2=...=tk-2,   tn+2 =... =tn+k
 *  Neumann boundary:
 *  t0=t1=...=tk-3, tn+3=...=tn+k
 *     
 *  boundary at minimm energy enforces energy conservation
 * 
 * @param gt grid type, also referred to as grid_style. Only used for zero_degree_0 and zero_degree_inf. TODO really need to refactor around how we are using gt.
 */

std::vector<double> BasisSet::set_knot(const GridSpacing& gt, FeatureRegimes& rgm, bool trial, bool do_not_update_regions){ // 'trial' isn't used by program currently
    if (!do_not_update_regions){
        update_regions(rgm);
    }
    Z_0 = gt.zero_degree_0;
    Z_inf = gt.zero_degree_inf;
    if( BSPLINE_ORDER - Z_0 < 0){
        std::cerr <<"Failed to set boundary condition at 0, defaulting to "<<0<<std::endl;
        Z_0 = 0;
    }
    if (BSPLINE_ORDER - Z_inf < 0){
        std::cerr <<"Failed to set boundary condition at infinity, defaulting to "<<BSPLINE_ORDER<<std::endl;
        Z_inf = BSPLINE_ORDER;
    } 

    size_t start = BSPLINE_ORDER-Z_0-1;

    std::sort(regions.begin(),regions.end()); // Sort by min energy.
    // Set first energy to minimum of all grids present.
    std::vector<double> boundary_E = {regions[0].get_E_min()};        

    std::cout<<"[ Dynamic Knot ] Using splines of "<<BSPLINE_ORDER<<"th order "<<std::endl;
    std::cout<<"[ Dynamic Knot ] (i.e. piecewise polynomial has leading order x^"<<BSPLINE_ORDER-1<<")"<<std::endl;

    std::vector<double> new_knots;
    new_knots.reserve(400+BSPLINE_ORDER+1);

    // Set the k - zero_deg repeated knots
    for (size_t i=0; i<start; i++) {
        new_knots.push_back(_min);
    }
    // minimum (0)
    new_knots.push_back(_min);

    size_t i=start;
    // Dynamically choose next grid point. Idea: We move to next point based on whichever region's rule gives the closest one. 
    while(true){
        i++;
        // Find the smallest grid point that regions that the point is within want.
        double next_point = INFINITY;
        for (size_t r = 0; r < regions.size(); r ++){
            double point =  regions[r].get_next_knot(new_knots[i-1],i==start+1); // If first non_zero point, we get the minimum of the regions since it is safe to do so (we don't do it at other points so that no knots are closer together than the regions it is within/neighbouring.)
            if (point > new_knots[i-1]){
                next_point = min(next_point,point);
            }
        }
        if(next_point == INFINITY){
            // We are outside the max, add max point + boundary thing and finish.
            new_knots.push_back(_max);
            new_knots.push_back(_max+500/Constant::eV_per_Ha);
            break;
        }
        assert(next_point > new_knots[i-1]);
        assert(i < 999);
        new_knots.push_back(next_point);
    }
    size_t tmp_num_funcs = new_knots.size()-(start+1)-1;
    // t_{n+1+z_infinity} has been set now. Repeat it until the end.
    for (size_t i= tmp_num_funcs + 1 + Z_inf; i< tmp_num_funcs+BSPLINE_ORDER; i++) {
        new_knots.push_back(new_knots[tmp_num_funcs + Z_inf]);
    }
    
    if (!trial){
        num_funcs = tmp_num_funcs;
        knot = new_knots;
    }
    #ifdef DEBUG
    std::cerr<<"Knot: [ ";
    for (auto& ki : new_knots) {
        std::cerr<<ki<<" ";
    }
    std::cerr<<"]\n";
    #endif    
    return new_knots;
}



void BasisSet::manual_set_knot(const GridSpacing& gt){
    Z_0 = gt.zero_degree_0;
    Z_inf = gt.zero_degree_inf;
    if( BSPLINE_ORDER - Z_0 < 0){
        std::cerr <<"Failed to set boundary condition at 0, defaulting to "<<0<<std::endl;
        Z_0 = 0;
    }
    if (BSPLINE_ORDER - Z_inf < 0){
        std::cerr <<"Failed to set boundary condition at infinity, defaulting to "<<BSPLINE_ORDER<<std::endl;
        Z_inf = BSPLINE_ORDER;
    } 

    size_t start = BSPLINE_ORDER-Z_0-1;
    //size_t num_int = num_funcs + Z_inf - start;

    /// Wiggle destroyer.  Continuous piecewise.
    // e.g. 6k eV photon beam. Have delta-like stuff in MB and dirac.
    // Idea: Define e.g. four energy regions: (eV) <-0--low--10--mid--200-trans.--2000--high--10000--surplus->
    // We want to have higher density in the mid and high energy regions, and low in the low, transition, and surplus regions.
    // Form of E = An^p + B

    // Force first point to be 0
    _manual_region_bndry_index.insert(_manual_region_bndry_index.begin(), 0);
    _manual_region_bndry_energy.insert(_manual_region_bndry_energy.begin(), 0);
    _region_powers.insert(_region_powers.begin(), 0);
    assert(_manual_region_bndry_energy.size() == _manual_region_bndry_index.size());
    assert(_region_powers.size() == _manual_region_bndry_index.size() - 1);

    std::vector<double> hyb_powlaw_factor (_manual_region_bndry_index.size() - 1,0.);
    // Params that define region boundaries:
    int n_M, n_N; //Index of first point in region/next region.  
    double E_M, E_N; // Corresponding energies.
    double p;        // power law that sparseness of grid points in each region follows.
    // R region boundaries and 1 power law for each region --> R - 1 power laws.
    // 
    for (size_t rgn=1; rgn< _region_powers.size(); rgn++){   // (0 to pow of 0 gives 1 so skip n=0 point.).
        n_M = _manual_region_bndry_index[rgn];
        E_M = _manual_region_bndry_energy[rgn];
        n_N = _manual_region_bndry_index[rgn+1];
        E_N = _manual_region_bndry_energy[rgn+1];
        p = _region_powers[rgn];
        
        hyb_powlaw_factor[rgn] = (E_N - E_M)/pow(n_N-n_M, p);
    }
    // classic power law
    bool classic_mode = false;
    #ifdef CLASSIC_MANUAL_GRID
        classic_mode = true;
    #endif
    double p_powlaw;
    double A_powlaw;
    if (classic_mode){
        double transition_e  = 2000/Constant::eV_per_Ha;
        size_t num_low = 35;
        size_t M = num_low + 1;
        p_powlaw = (log(_max-_min) - log(transition_e - _min))/(log(1.*num_funcs/M));
        A_powlaw = (_max - _min)/pow(num_funcs, p_powlaw);
    }

    std::cout<<"[ Manual Knot ] Using splines of "<<BSPLINE_ORDER<<"th order "<<std::endl;
    std::cout<<"[ Manual Knot ] (i.e. piecewise polynomial has leading order x^"<<BSPLINE_ORDER-1<<")"<<std::endl;

    std::cout<<"[ Manual Knot ] Using power-law exponents: "; 
    for (size_t j = 0; j < _region_powers.size(); j++){
        cout << _region_powers[j] << ", ";
    }
    std::cout << std::endl; 
    std::cout<<"[ Manual Knot ] Using power-law factors: ";
    for (size_t j = 0; j < hyb_powlaw_factor.size(); j++){
        cout << hyb_powlaw_factor[j] << ", ";
    }       
    std::cout << std::endl;         


    knot.resize(num_funcs+BSPLINE_ORDER+1); // num_funcs == n + 1

    
    // Set the k - zero_deg repeated knots
    for (size_t i=0; i<start; i++) {
        knot[i] = _min;
    }

    // TODO:
    // - Move this grid-construction code to GridSpacing.hpp
    // - Refactor to incorporate boundaries max and min into the GridSpacing object
    // - Consider better method to deal with MB divergence at low T.

    // Keep in mind: At i=0, the value of knot will still be _min
    bool used_good_grid = false;

    for(size_t i=start; i<=num_funcs + Z_inf; i++) {
        if (!classic_mode){
            used_good_grid = true;        
            // Generalised custom spacing
            //  Get the index of the region (rgn) that this point is part of.
            size_t rgn = 0;
            for( ; rgn < _region_powers.size(); rgn++){
                if( i - start < static_cast<size_t>(_manual_region_bndry_index[rgn+1])
                    || rgn == _region_powers.size() - 1){
                    break;
                    }
            }
            n_M = _manual_region_bndry_index[rgn];
            E_M = _manual_region_bndry_energy[rgn];
            p = _region_powers[rgn];            
            // Calculate the knot energy
            knot[i] = hyb_powlaw_factor[rgn] * pow(i - n_M - start, p) + E_M;
        }
        else{
            knot[i] = A_powlaw * pow(i-start, p_powlaw) + _min;
        }        
    }

    if(!used_good_grid) std::cout << "WARNING, this grid type is obsolete, if a dynamic grid is not possible, using a targeted grid is favourable. See README." <<std::endl; 
    // t_{n+1+z_infinity} has been set now. Repeat it until the end.
    for (size_t i= num_funcs + 1 + Z_inf; i< num_funcs+BSPLINE_ORDER; i++) {
        knot[i] = knot[num_funcs + Z_inf];
    }
    

    #ifdef DEBUG
    std::cerr<<"Knot: [ ";
    for (auto& ki : knot) {
        std::cerr<<ki<<" ";
    }
    std::cerr<<"]\n";
    #endif
}

/**
 * @brief Sets up the B-spline knot to have the appropriate shape (respecting boundary conditions)
 * @details 
 * 
 * @param _num_funcs the number of B-splines to use in the basis. 
 * @param min minimum energy
 * @param max maximum energy
 * @param gt gt.zero_degree_0: The number of derivatives to set to zero: 0 = open conditions, 1=impose f(0)=0, 2=impose f(0)=f'(0) =0
 * @param  
 */
void BasisSet::set_parameters(const GridSpacing& gt, ManualGridBoundaries& manual_elec_grid_regions, FeatureRegimes& regimes) {
    this->_manual_region_bndry_index = manual_elec_grid_regions.bndry_idx;   // TODO: refactor, can replace min and max with array of region start and region ends. -S.P.
    this -> _manual_region_bndry_energy = manual_elec_grid_regions.bndry_E;   
    this -> _region_powers = manual_elec_grid_regions.powers;   

    // Idea: Want num_funcs usable B-splines
    // If grid has num_funcs+k+1 points, num_funcs of these are usable splines (num_int = num_funcs + Z_inf - start)
    // grid layout for open boundary:
    // t0=t1=t2=...=tk-1,   tn+1 = tn+2 =... =tn+k
    // homogeneous Dirichlet boundary:
    // t0=t1=t2=...=tk-2,   tn+2 =... =tn+k
    // Neumann boundary:
    // t0=t1=...=tk-3, tn+3=...=tn+k
    
    // boundary at minimum energy enforces energy conservation 
    this->_min = 0;
    if (gt.mode == GridSpacing::dynamic){
        assert(regions.size() > 0);
        this->_max = dynamic_max_inner_knot();  
        set_knot(gt,regimes);        
    }
    else{    
        assert(gt.mode == GridSpacing::manual);
        this->_max = manual_elec_grid_regions.bndry_E.back();
        num_funcs = manual_elec_grid_regions.bndry_idx.back();
        manual_set_knot(gt);
    }    
    // std::cout << "Knots set to: ";
    // for (double e: knot)
    //     std::cout << e*Constant::eV_per_Ha << ' ';
    // cout << std::endl;
    compute_overlap(num_funcs);
}


void BasisSet::compute_overlap(size_t num_funcs){
    std::vector<Eigen::Triplet<double>> tripletList;
    tripletList.reserve(BSPLINE_ORDER*2*num_funcs);
    
    // NB: S is -> Eigen::MatrixXd S(num_funcs, num_funcs);
    for (size_t i=0; i<num_funcs; i++) {
        for (size_t j=i+1; j<num_funcs; j++) {
            double tmp=overlap(i, j);
            if (tmp != 0) {
                tripletList.push_back(Eigen::Triplet<double>(i,j,tmp));
                tripletList.push_back(Eigen::Triplet<double>(j,i,tmp));
                // S(i,j) = tmp;
                // S(j,i) = tmp;
            }
        }
        // S.insert(i,i) = overlap(i,i);
        tripletList.push_back(Eigen::Triplet<double>(i,i,overlap(i, i)));
        // S(i,i) = overlap(i,i);
    }
    /// Compute overlap matrix   //s\/ quick and fast note - if we use the grid basis (with num_func splines) the splines overlap. S takes us to a spline-independent basis I believe -S.P.
    Eigen::SparseMatrix<double> S(num_funcs, num_funcs);   // (num_funcs rows and cols). 
    S.setFromTriplets(tripletList.begin(), tripletList.end());
    // Precompute the LU decomposition
    linsolver.analyzePattern(S);
    linsolver.isSymmetric(true);
    linsolver.factorize(S);
    if(linsolver.info()!=Eigen::Success) {
        throw runtime_error("Factorisation of overlap matrix failed!");
    }

    avg_e.resize(num_funcs);
    log_avg_e.resize(num_funcs);
    areas.resize(num_funcs);
    for (size_t i = 0; i < num_funcs; i++) {
        // Chooses the 'center' of the B-spline
        avg_e[i] = (this->supp_max(i) + this->supp_min(i))/2 ;
        log_avg_e[i] = log(avg_e[i]);
        double diff = (this->supp_max(i) - this->supp_min(i))/2;
        // Widths stores the integral of the j^th B-spline
        areas[i] = 0;
        for (size_t j=0; j<10; j++){
            areas[i] += gaussW_10[j]*(*this)(i, gaussX_10[j]*diff+ avg_e[i]);
        }
        areas[i] *= diff;
    }
}


void BasisSet::set_parameters(FeatureRegimes& regimes, std::vector<double> new_grid_knots) {
    assert(Z_0 !=-1 && "Z_0 hasn't been set");
    assert( BSPLINE_ORDER - Z_0 >= 0);      
    update_regions(regimes);
    knot = new_grid_knots;
    num_funcs = knot.size()-(BSPLINE_ORDER-Z_0)-1; //knot.size()- (BSPLINE_ORDER+1); 
    compute_overlap(num_funcs);
}

/**
 * @brief Sinv = S inverse. Changes basis of vector deltaf from spline basis to grid basis and vice versa I think? -S.P.
 * @details
 */
Eigen::VectorXd BasisSet::Sinv(const Eigen::VectorXd& deltaf) {
    // Solves the linear system S fdot = deltaf   (returning f dot -S.P.)
    auto tmp = linsolver.solve(deltaf);
    return tmp;
}

/**
 * @brief Sinv = S inverse.
 * @details 
 */
Eigen::MatrixXd BasisSet::Sinv(const Eigen::MatrixXd& J) {
    // Solves the linear system S fdot = deltaf
    return linsolver.solve(J);
}

double BasisSet::raw_bspline(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
}

double BasisSet::raw_Dbspline(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    return BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]);
}

double BasisSet::operator()(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    // Returns the i^th B-spline of order BSPLINE_ORDER
    
    static_assert(USING_SQRTE_PREFACTOR);
    return pow(x,0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    // return BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    
}

// Returns the i^th B-spline of order BSPLINE_ORDER evaluated at x (safe version)
double BasisSet::at(size_t i, double x) const {
    if(i >= num_funcs || i<0) return 0;
    return operator()(i, x);
}

// Returns the first derivative of the i^th B-spline of order BSPLINE_ORDER evaluated at x
double BasisSet::D(size_t i, double x) const {
    assert(i < num_funcs && i >= 0);
    static_assert(BSPLINE_ORDER > 1);
    
    static_assert(USING_SQRTE_PREFACTOR);
    return pow(x,0.5)*BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]) + 0.5*pow(x,-0.5)*BSpline::BSpline<BSPLINE_ORDER>(x, &knot[i]);
    
    // return BSpline::DBSpline<BSPLINE_ORDER>(x, &knot[i]);
    
}

double BasisSet::overlap(size_t j,size_t i) const{
    // computes the overlap of the jth raw bspline with the ith modified bspline
    double b = min(supp_max(i), supp_max(j));
    double a = max(supp_min(i), supp_min(j));
    if (a >= b) return 0;
    double tmp=0;
    for (int m=0; m<10; m++) {
        double e = gaussX_10[m]*(b-a)/2 + (b+a)/2;
        tmp += gaussW_10[m]*raw_bspline(j, e)*(*this)(i, e);
        // tmp += gaussW_10[m]*(*this)(j, e)*(*this)(i, e);
    }
    tmp *= (b-a)/2;
    return tmp;
}

// returns a vector of intervals (bottom, K1), (K1, K2), (K2, top)
// such that the splines are polynomial on these intervals
std::vector<std::pair<double,double>> BasisSet::knots_between(double bottom, double top) const {
    std::vector<std::pair<double,double> > intervals(0);
    // NOTE:
    // These calls can and should should be replaced with std::lower_bound calls.
    if (bottom >= top) return intervals;
    size_t bottomidx = lower_idx(bottom, knot);
    size_t topidx = lower_idx(top, knot) + 1;
    if (topidx >= knot.size()) topidx = knot.size() -1;
    for (size_t i = bottomidx; i<= topidx; i++){
        intervals.push_back(std::pair<double,double>(knot[i],knot[i+1]));
    }
    return intervals;
}