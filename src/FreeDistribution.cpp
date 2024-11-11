/**
 * @file FreeDistribution.cpp
 * @brief @copybrief FreeDistribution.h
 * @details Grid points are synonymous w/ knot points for the B-splines. Note that this code expands 
 * on the original AC4DC model by introducing these grid points, allowing for an approximation of a 
 * continuous electron density distribution along the energy basis. 
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

#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include "SplineIntegral.h"
#include <Eigen/StdVector>

// #define NDEBUG
// to remove asserts

// Initialise static things

size_t Distribution::CoulombLog_cutoff=0;
double Distribution::CoulombDens_min=0;
std::vector<indexed_knot> Distribution::knots_history;
// These variables are modified by set_distribution_STATIC_ONLY or by dynamic grid updates when set_basis is called.
SplineIntegral Distribution::basis; 
size_t Distribution::size=0;  
size_t Distribution::num_continuums = 1;

#ifdef FIND_INITIAL_DIRAC 
    bool Distribution::reset_on_next_grid_update = true;  // TODO duct tape implementation...
#else
    bool Distribution::reset_on_next_grid_update = false;
#endif     


// Psuedo-constructor thing (Maybe not anymore... -S.P.)
void Distribution::set_basis(size_t step, GridSpacing grid_style, Cutoffs param_cutoffs, FeatureRegimes regimes, ManualGridBoundaries elec_grid_regions){
    // Defines a grid of num_funcs points (if manual, thsi is the total number of free-electron grid points specified in the .mol file.)
    // where num_funcs is the number of non-boundary (i.e. "usable") splines/knots.
    basis.set_parameters(grid_style, elec_grid_regions,regimes);
    knots_history.push_back(indexed_knot{step,get_knot_energies()});
    Distribution::size=basis.num_funcs;
    Distribution::CoulombLog_cutoff = basis.i_from_e(param_cutoffs.transition_e);
    #ifdef INFINITE_COULOMBLOG_CUTOFF
    Distribution::CoulombLog_cutoff = basis.i_from_e(INFINITY);
    #endif

    Distribution::CoulombDens_min = param_cutoffs.min_coulomb_density;
    cout<<"[ Free ] Estimating lnLambda for energies below ";
    cout<<param_cutoffs.transition_e<<" Ha"<<endl;
    cout<<"[ Free ] Neglecting electron-electron below density of n = "<<CoulombDens_min<<"au^-3"<<endl;
}


void Distribution::set_basis(size_t step, Cutoffs param_cutoffs, FeatureRegimes regimes, std::vector<double> knots, bool update_knot_history){
    // Defines a grid of num_funcs points (if manual, this is the total number of free-electron grid points specified in the .mol file.)
    // where num_funcs is the number of non-boundary (i.e. "usable") splines/knots.
    basis.set_parameters(regimes, knots);
    if(update_knot_history){knots_history.push_back(indexed_knot{step,get_knot_energies()});}
    Distribution::size=basis.num_funcs;
    Distribution::CoulombLog_cutoff = basis.i_from_e(param_cutoffs.transition_e);
    #ifdef INFINITE_COULOMBLOG_CUTOFF
    Distribution::CoulombLog_cutoff = basis.i_from_e(INFINITY);
    #endif    
    Distribution::CoulombDens_min = param_cutoffs.min_coulomb_density;
    cout<<"[ Free ] Estimating lnLambda for energies below ";
    cout<<param_cutoffs.transition_e<<" Ha"<<endl;
    cout<<"[ Free ] Neglecting electron-electron below density of n = "<<CoulombDens_min<<"au^-3"<<endl;
}

void Distribution::set_spline_factors(size_t _c, vector<double> new_spline_factors){
    f_array[_c] = new_spline_factors; 
}
void Distribution::set_distribution_STATIC_ONLY(size_t _c, vector<double> new_knot, vector<double> new_spline_factors) {
    f_array[_c] = new_spline_factors; 
    load_knot(new_knot);
    f_array[_c].resize(size);
}
void Distribution::load_knot(vector<double> loaded_knot) {
    loaded_knot = get_trimmed_knots(loaded_knot); // Remove boundary knots
    size = loaded_knot.size();
    basis = loaded_knot; // note the overload!
}

// Adds Q_eii to the parent Distribution
void Distribution::get_Q_eii (size_t _c, Eigen::VectorXd& v, size_t a, const bound_t& P, const int & threads) const {
    assert(basis.has_Qeii());
    assert(P.size() == basis.Q_EII[a].size());
    assert((unsigned) v.size() == size);
    
    //double v_copy [size] = {0};
    
    for (size_t J=0; J<size; J++) {
        double vj=0;
        #pragma omp parallel for num_threads(threads) reduction(+ : vj) // Do NOT use collapse(2), it's about twice as slow.
        for (size_t xi=0; xi<P.size(); xi++) {
        // Loop over configurations that P refers to
            for (size_t K=0; K<size; K++) {
                v_copy[J] += P[xi]*f_array[_c][K]*basis.Q_EII[a][xi][J][K];
            }
        }
        v[J] += vj;
        //v += Eigen::Map<Eigen::VectorXd>(v_copy,size);
    }
}

/**
 * @brief 
 * @details d/dt f = Q_B[f](t)  
 * 
 * @param v 
 * @param a 
 * @param P probabilities of each atomic state;  d/dt P[i] = \sum_i=1^N W_ij - W_ji ~~~~ P[j] = d/dt(average-atomic-state)
 */
// Puts the Q_TBR changes in the supplied vector v
void Distribution::get_Q_tbr (size_t _c, Eigen::VectorXd& v, size_t a, const bound_t& P, const int & threads) const {
    assert(basis.has_Qtbr());
    assert(P.size() == basis.Q_TBR[a].size());
    double v_copy [size];
    #pragma omp parallel for num_threads(threads) reduction(+ : v_copy) collapse(2)
    for (size_t eta=0; eta<P.size(); eta++) {          // eta -> configuration
        // Loop over configurations that P refers to
        for (size_t J=0; J<size; J++) {                   // J -> grid point
            for (auto& q : basis.Q_TBR[a][eta][J]) {   // Thousands of iterations for each J - S.P.
                v_copy[J] += q.val * P[eta] * (f_array[_c][q.K] * f_array[0][q.L] + f_array[0][q.K] * f_array[_c][q.L])*0.5; //Correct?
                //v_copy[J] += q.val * P[eta] * f_array[_c][q.K] * f_array[_c][q.L];

            }
        }
    }
    v += Eigen::Map<Eigen::VectorXd>(v_copy,size);
}

// Puts the Q_EE changes in the supplied vector v
void Distribution::get_Q_ee(size_t _c, Eigen::VectorXd& v, const int & threads) const {
    assert(basis.has_Qee());
    double CoulombLog = this->CoulombLogarithm();
    // double CoulombLog = 3.4;
    if (CoulombLog <= 0) return; // skip calculation if LnLambda vanishes

    if (isnan(CoulombLog) || CoulombLog > 11.5) CoulombLog = 11.5;
    // double LnLambdaD = 0.5*log(this->k_temperature()/4/Constant::Pi/this->density());
    // if (isnan(LnLambdaD)) LnLambdaD = 11;
    // cerr<<"LnDebLen = "<<LnLambdaD<<endl;
    // A guess. This should only happen when density is zero, so Debye length is infinity.
    // Guess the sample size is about 10^5 Bohr. This shouldn't ultimately matter much.   /// Attention - S.P. // Actually it seems this isn't active? Something something fences on roads.
    double v_copy [size] = {0}; 
    #pragma omp parallel for num_threads(threads) reduction(+ : v_copy)  collapse(2)       
    for (size_t J=0; J<size; J++) {
        for (size_t K=0; K<size; K++) {
            for (auto& q : basis.Q_EE[J][K]) {
                 //v_copy[J] += q.val * f_array[0][K] * f_array[0][q.idx] * CoulombLog;  
                 v_copy[J] += q.val * (f_array[_c][K] * f_array[0][q.idx] + f_array[0][K] * f_array[_c][q.idx])*0.5 * CoulombLog; 
            }
        }
    }
    v += Eigen::Map<Eigen::VectorXd>(v_copy,size);
}




// ============================
// utility

// See add_density_destribution comments
// Expects T to have units of Ha
void Distribution::add_maxwellian(size_t _c, double T, double N) {
    Eigen::VectorXd v(size);
    for (size_t i=0; i<size; i++) {
        v[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        for (size_t j=0; j<10; j++) {
            double e = (b-a)/2 *gaussX_10[j] + (a+b)/2;                 
            v[i] +=  gaussW_10[j]*exp(-e/T)*basis(i, e)*pow(e,0.5);; 
        }
        v[i] *= (b-a)/2;
        v[i] *= N*pow(T, -1.5)*2/pow(Constant::Pi,0.5);
    }
    Eigen::VectorXd u = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++) {
        this->f_array[_c][i] += u[i];
    }
}

// Assumes new basis has same B spline order.
void Distribution::transform_basis(std::vector<double> new_knots){
    int new_basis_order = basis.BSPLINE_ORDER;
    //// Get knots that have densities
    int num_new_splines = static_cast<int>(get_trimmed_knots(new_knots).size());   // TODO replace get_trimmed_knots with get_num_funcs?. 
    std::vector<std::vector<std::vector<double>>> new_densities_container;
    for (size_t _c = 0; _c < num_continuums; _c++){
        //// Compute densities for knots
        std::vector<std::vector<double>> new_densities(num_new_splines, std::vector<double>(64, 0));

        // Cackle and iterate through each new spline.
        for (size_t i=0; static_cast<int>(i)<num_new_splines; i++){
            // Use current basis to generate the density terms for gaussian integration at for each basis point.   
            // Black magic. ଘ(੭ˊᵕˋ)੭.*･｡ﾟ
            double a = new_knots[i];                  // i.e. <new_basis>.supp_min(i);
            double b = new_knots[i+new_basis_order];  // i.e. <new_basis>.supp_max(i);        
            for(size_t j=0; j < 64; j++){
                double e = (b-a)/2 *gaussX_64[j] + (a+b)/2;
                new_densities[i][j] = (*this)(_c,e);  
            }
        }
        // Change distribution to empty one in new basis.    
        new_densities_container.push_back(new_densities);
    }
    for (size_t _c = 0; _c < num_continuums; _c++){
        vector<double> new_f(num_new_splines,0);
        f_array[_c] = new_f; 
    }
    load_knot(new_knots);
    for (size_t _c = 0; _c < num_continuums; _c++){
        f_array[_c].resize(size);
        // Add densities in new basis.
        add_density_distribution(_c, new_densities_container[_c]);

    }
}

void Distribution::add_density_distribution(size_t _c, vector<vector<double>> densities){

    assert(densities.size() == size);
    Eigen::VectorXd v(size);
    std::vector<double> energies = get_trimmed_knots(get_knot_energies());//basis.avg_e;//get_trimmed_knots(get_knot_energies());
    for (size_t i=0; i<size; i++) {
        v[i] = 0;
        double a = basis.supp_min(i);
        double b = basis.supp_max(i);
        // This uses gaussian quadrature (https://pomax.github.io/bezierinfo/legendre-gauss.html)
        // to approximate the integral of bspline_i.density_i between the boundaries of the splines, but with some basis shenanigans to spice it up.
        // Thus v[i] would correspond to the total num electrons contributed by spline i. -S.P.
        for (size_t j=0; j<64; j++) {
            double e = (b-a)/2 *gaussX_64[j] + (a+b)/2;    // e = element E_i of gaussian quadrature sum
            v[i] += gaussW_64[j]*basis.raw_bspline(i, e)*densities[i][j]; 
        }
        v[i] *= (b-a)/2;//*=densities[i]*(b-a)/2;
    }
    // Su = v, where u[i] is the electron density (in the grid basis?), S is the (sparse) overlap matrix, and v[i] is the num electrons the spline represents. 
    Eigen::VectorXd u = this->basis.Sinv(v);
    for (size_t i=0; i<size; i++) {
        this->f_array[_c][i] += u[i];
    }
}

std::vector<double> Distribution::get_trimmed_knots(std::vector<double> knots){
    // Remove boundary knots  // TODO CRITICAL this is possibly BAD as it might not work with non-default Z_0 and Z_inf?
    while(knots[0] <= basis.min_elec_e() && knots.size() > 0){
        knots.erase(knots.begin());
    }        
    while (knots.back() > basis.max_elec_e() && knots.size() > 0){
        knots.erase(knots.end() - 1);
    }   
    return knots;
}


/**
 * @brief Loads the static variables relating to the basis that are changed by grid updates.
 * @param step_idx Step to load.
 * @return returns the knot that was loaded for convenience.
 */
std::vector<double> Distribution::load_knots_from_history(size_t step_idx){
    std::vector<double> loaded_knot = get_knots_from_history(step_idx);
    load_knot(loaded_knot);
    return loaded_knot;
}

size_t Distribution::most_recent_knot_change_idx(size_t step_idx){
    size_t most_recent_update = 0;
    // Find the step of the most recent knot update as of step_idx. (will return same step if given the step of the change.)
    for(auto elem: knots_history){
        if (elem.step > step_idx) 
            break; 
        most_recent_update = elem.step;
    }
    return most_recent_update;
}


size_t Distribution::next_knot_change_idx(size_t step_idx){
    size_t next_knot_update;
    // Find the step of the next knot update as of step_idx. (will return next update idx if given the step of the change.)
    for(auto elem: knots_history){
        if (elem.step > step_idx){
            next_knot_update = elem.step;
            break; 
        }
    } 
    return next_knot_update;
}


std::vector<double> Distribution::get_knots_from_history(size_t step_idx){
    vector<double> loaded_knot;
    // Find the most recent grid update's knots as of step_idx.  (if knot was changed at step idx will return the knot loaded at that step.)
    for(auto elem: knots_history){
        if (elem.step > step_idx) 
            break;
        loaded_knot = elem.energy;
    }
    return loaded_knot;
}
// ==========================
// IO

// Returns knot energies in eV
std::string Distribution::output_knots_eV() {
    std::stringstream ss;
    for (size_t i=0;i<basis.gridlen(); i++) {
        ss << basis.grid(i)*Constant::eV_per_Ha<<" ";
    }
    return ss.str();
}

// Returns energies in eV
std::string Distribution::output_energies_eV(size_t num_pts) {
    std::stringstream ss;
    size_t pts_per_knot = num_pts / basis.gridlen();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    for (size_t i=0; i<basis.gridlen()-1; i++){
        double e = basis.grid(i);
        double de = (basis.grid(i+1) - e)/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            ss << e*Constant::eV_per_Ha<<" ";
            e += de;
        }
    }
    
    return ss.str();
}

/**
 * @brief 
 * @param num_pts 
 * @param reference_knots knots to base the outputted energies off.
 * @return 
 */
std::string Distribution::output_densities(size_t _c, size_t num_pts,std::vector<double> reference_knots) const {    
    size_t pts_per_knot = num_pts / reference_knots.size();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    // Iterate through each knot
    std::stringstream ss;
    const double units = 1./Constant::eV_per_Ha/Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;
    for (size_t i=0; i<reference_knots.size()-1; i++){
        // output pts_per_knot points between this knot to the next.
        double e = reference_knots[i];
        double de = (reference_knots[i+1] - reference_knots[i])/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            ss << (*this)(_c,e)*units<<" ";
        e += de;
        }
    }
    return ss.str();
}


// for live plotting 
////
std::vector<double> Distribution::get_energies_eV(size_t num_pts) {
    std::vector<double> energies(num_pts);
    size_t pts_per_knot = num_pts / basis.gridlen();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        energies.resize(basis.gridlen());
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    for (size_t i=0; i<basis.gridlen()-1; i++){
        double e = basis.grid(i);
        double de = (basis.grid(i+1) - e)/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            energies[i*pts_per_knot + j]=e*Constant::eV_per_Ha;
            e += de;
        }
    }
    return energies;
}
std::vector<double> Distribution::get_densities(size_t _c, size_t num_pts,std::vector<double> reference_knots) const {    
    std::vector<double> densities(num_pts);
    size_t pts_per_knot = num_pts / reference_knots.size();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        densities.resize(basis.gridlen());
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    // Iterate through each knot
    const double units = 1./Constant::eV_per_Ha/Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;
    for (size_t i=0; i<reference_knots.size()-1; i++){
        // output pts_per_knot points between this knot to the next.
        double e = reference_knots[i];
        double de = (reference_knots[i+1] - reference_knots[i])/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            densities[i*pts_per_knot + j] = (*this)(_c,e)*units;
        e += de;
        }
    }
    return densities;
}
///////



/**
 * @brief Returns the spline-interpolated electron density of states at any energy (In atomic units, i.e. 1/(Ha.a0^3)).
 * Multiplying by E returns the electron energy density.
 * @details This does a dot product of the "Frobenius-treated" basis with each spline's expansion coefficient.
 * @param e energy to get density from.
 * @return double 
 */
double Distribution::operator()(size_t _c, double e) const{
    double tmp=0;
    for (size_t j = 0; j < size; j++) {
        tmp += basis(j, e)*f_array[_c][j];
    }
    return tmp;
}
/*
void Distribution::addDeltaSpike(double e, double N) {
    int idx = basis.lower_i_from_e(e);
    assert(idx >= 0 && (unsigned) idx+1 < size);
    double A1 = basis.areas[idx];
    double A2 = basis.areas[idx+1];
    double E1 = basis.avg_e[idx];
    double E2 = basis.avg_e[idx+1];
    // Solves the matrix equation
    // A1 * a + A2 * b = N
    // A1 E1 a + A2 E2 b = (A1 + A2) e
    
    double det = (E2 - E1)*A1*A2;
    // Inverse matrix is 
    // 1/det * [ A2 E2   -A2 ]
    //         [-A1 E1    A1 ]
    
    f[idx] +=     ( A2 * E2 * N   -  A2 * (A1 + A2)* e) / det;
    f[idx + 1] += (- A1 * E1 * N  +  A1 * (A1 + A2)* e) / det;
}*/

void Distribution::addDeltaSpike(size_t a, double e, double N) {
    int idx = basis.i_from_e(e);
    f_array[0][idx] += N/basis.areas[idx];
    #ifndef TRACK_SINGLE_CONTINUUM
    f_array[a+1][idx] += N/basis.areas[idx];
    #endif
}

/**
 * @brief f <-- f + df/dt|basis  
 * @details NOT applying a dirac delta. 
 * @param v 
 */
void Distribution::applyDeltaF(size_t a,const Eigen::VectorXd& v,const int & threads) {
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v)); 
    #pragma omp for schedule(dynamic) nowait
    for (size_t i=0; i<size; i++) {
        f_array[0][i] += u[i];
        #ifndef TRACK_SINGLE_CONTINUUM
        f_array[a+1][i] += u[i];
        #endif
    }
}

void Distribution::applyDeltaF_element_scaled(size_t a,const Eigen::VectorXd& v,const int & threads){

    #ifdef TRACK_SINGLE_CONTINUUM
    throw std::runtime_error("Function was called that is incompatible with single continuum tracking.");
    #endif
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v)); 
    #pragma omp for schedule(dynamic) nowait
    for (size_t i=0; i<size; i++) {
        u[i] *= f_array[a+1][i];
        f_array[0][i] += u[i];
        f_array[a+1][i] += u[i];
    }   
}


// - 3/sqrt(2) * 3 sqrt(e) * f(e) / R_
// Very rough approximation used here
void Distribution::addLoss(size_t a, const Distribution& d, const LossGeometry &l, double rho) {
    // f += "|   i|   ||   |_"

    double escape_e = 1.333333333*Constant::Pi*l.L0*l.L0*rho;
    for (size_t i=basis.i_from_e(escape_e); i<size; i++) {
        
        f_array[0][i] -= d[0][i] * sqrt(basis.avg_e[i]/2) * l.factor();
        #ifndef TRACK_SINGLE_CONTINUUM
        f_array[a+1][i] -= d[a+1][i] * sqrt(basis.avg_e[i]/2) * l.factor();
        #endif
    }
}

// Calculating the filtration rate of photoelectrons
// The way that the photoelectrons leave will be somewhere between the rate corresponding to the 
// limit of non-interacting, and the limit of Sanders' implementation, where the velocities of electrons
// between steps are independent due to a high collision rate/low mean free path .
// The difference being, one is a flat rate, whereas one is a percentage rate.
// for a 1000eV temperature *equilibrium* plasma with carbon density and all electrons ionised, 
// the collision rate is < 0.6785 /fs. Electrons are travelling at ~ 50-100 nm/fs, in a < 10-100 nm crystal...
// This discussion is all likely very negligible.
/**
 * @brief 
 * 
 * @param d  distribution on previous step 
 * @param bg background distribution
 * @param l 
 */
void Distribution :: addFiltration(size_t a, const Distribution& d, const Distribution& bg,const LossGeometry &l){
    // first order approximation, valid if step size * vel. much smaller than volume
    // For crystals our surface is flat, so the issue is 
    //return 0.5*v;  
    
    // For a given particle within some volume and velocity v, the time for it to leave the volume is
    // equivalent if we consider the volume to have -v instead. Indeed, we can do this with all particles
    // in the volume with the same v. Assuming homogeneity and a spherical symmetry (which it roughly is for a small crystal anyway), we can
    // calculate the rate that particles leave the volume by considering an equivalent volume passing over 
    // an equivalent number of stationary particles. 
    // slice of sphere:
    // L0 is the characteristic length scale given in input file.
    // v = sqrt(2*E) as m_e is 1 in a.u..
    // dN/dt = V_sphere*v/(2*R) = v*(2/3)pi*R^2 * f[i]  
    double c = 137.036;
    double m = 1;       
    double factor = 2./3*Constant::Pi*pow(l.L0,2);
    for (size_t i=0; i < size; i++){
         double v = c*sqrt( 1 - pow((1+basis.avg_e[i]/(m*c*c)), -2) );
        f_array[0][i] += factor*v*(bg[0][i] - d[0][i]);
        #ifndef TRACK_SINGLE_CONTINUUM
        f_array[a+1][i] += factor*v*(bg[a+1][i] - d[a+1][i]);
        #endif
    }
}

// Distribution Distribution :: relativistic_loss(const Distribution& d, const LossGeometry &l, double rho) {
//     std::vector<double> knots = basis.get_knot();
//     std::vector<std::vector<double>> lost_densities(basis.num_funcs, std::vector<double>(64, 0));
//     for (size_t i=0; i<basis.num_funcs; i++){
//         // Generate the density terms for gaussian integration at each basis point.   
//         double a = knots[i];                  // i.e. <new_basis>.supp_min(i);
//         double b = knots[i+basis.BSPLINE_ORDER];  
//          for(size_t j=0; j < 64; j++){
//             double e = (b-a)/2 *gaussX_64[j] + (a+b)/2;
//                 // Get the density multiplied by the filtration rate.
//                 lost_densities[i][j] = (*this)(e) * filtration_rates(e);  
//         }
//     }    
//     Distribution loss_distribution;
//     loss_distribution = 0; 
//     loss_distribution.add_density_distribution(lost_densities);    
//     return loss_distribution;
// }

void Distribution::addDeltaLike(Eigen::VectorXd& v, double e, double height) {
    assert(e>basis.supp_min(0));
    for (size_t i=0;i<size; i++) {
        v[i] += basis(i, e) * height;
    }

}

double Distribution::norm(size_t _c) const{
    double x=0;
    //for (size_t _c = 0; _c < num_continuums; _c++){ 
    assert(f_array[_c].size() == Distribution::size);  // A very blessed check -S.P.
    for (auto&fi : f_array[_c]) {
        x+= fabs(fi);
    }
    //}
    return x;
}

double Distribution::density() const {
    double tmp=0;
    for (size_t i = 0; i < size; i++)
    {
        tmp += basis.areas[i]*f_array[0][i]; // index 0 as we are getting the density of the combined continuum  (we don't physically care about the density of the split continuums).
    }
    return tmp;
}

double Distribution::density(size_t cutoff) const {
    double tmp=0;
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.areas[i]*f_array[0][i];
    }
    return tmp;
}

// Returns an estimate of k*T
// based on kinetic energy density of the plasma below the index given by 'cutoff'
// unused
double Distribution::k_temperature(size_t cutoff) const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = this->density(cutoff);
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.avg_e[i]*f_array[0][i]*basis.areas[i];
    }
    return tmp*2./3./n;
}

// Returns an estimate of ln(N_D) based on density and low-energy arguments
double Distribution::CoulombLogarithm() const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = 0;
    for (size_t i = 0; i < CoulombLog_cutoff; i++) {
        tmp += basis.avg_e[i]*f_array[0][i]*basis.areas[i]; 
        n += f_array[0][i]*basis.areas[i];
    }
    double kT = tmp*2./3./n;
    double DebyeLength3 = pow(kT/4/Constant::Pi/n,1.5);
    n = this->density();
    if (n < this->CoulombDens_min) return 0;
    return log(4.*Constant::Pi* n*DebyeLength3);
}


// double Distribution::integral(size_t _c, double (g)(double)) {
//     double retval = 0;
//     for (size_t J = 0; J < Distribution::size; J++)
//     {
//         double a = basis.supp_min(J);
//         double b = basis.supp_max(J);
//         double tmp = 0;
//         for (size_t i=0;i<10; i++)
//         {
//             double e = gaussX_10[i]*(b-a)/2 + (b+a)/2;
//             tmp += gaussW_10[i]*basis(J, e)*g(e);
//         }
//         retval +=  tmp * f_array[_c][J] * (b-a)/2;
//     }
//     return retval;
// }

// Outputs the combined continuum.
ostream& operator<<(ostream& os, const Distribution& dist) {
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[0][i]/Constant::eV_per_Ha;
    }
    return os;
}



