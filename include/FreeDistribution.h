/**
 * @author Alaric Sanders & Spencer Passmore
 * @file FreeDistribution.h
 * @brief Defines the class Distribution that represents the energy distribution of free electrons.
 * @details Expansion of original plasma code as part of Sanders' continuum plasma extension.
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

#ifndef FREEDISTRIBUTION_CXX_H
#define FREEDISTRIBUTION_CXX_H


#include <sstream>
#include <assert.h>
#include <math.h>
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include "Constant.h"
#include "SplineIntegral.h"
#include "Dipole.h"
#include "GridSpacing.hpp"
#include "LossGeometry.hpp"
#include "config.h"

struct indexed_knot
{
    size_t step; // Step that the grid of B-spline knots was set/updated
    std::vector<double> energy; // knots
};

/**
 * @brief Electron distribution class.
 * @details Represents a statistical distribution of electron density. Internal units are atomic units.
 * @note F is referred to as the distribution throughout the code. F[i] returns the i'th spline factor, but
 * F.f_array[0], (in prior versions F.f), is the container that holds the spline factors for the total continuum (later indices are for tracking cascades initiated by each atomic species). F(e) expands out the basis to return the density at energy e divided by local knot width.
 */
class Distribution
{
public:
    Distribution() {
        f_array.resize(num_continuums);
        for (auto& f : f_array){
            f.resize(size);}   
    }


    inline std::vector<double>& operator[](size_t _c) {
        return this->f_array[_c];
    }

    inline std::vector<double> operator[](size_t _c) const{
        return this->f_array[_c];
    }

    // /**
    //  * @brief Get spline factors at 
    //  * @param n 
    //  * @return 
    //  */
    // inline double& operator()(size_t _c, size_t n) {
    //     return this->f_array[_c][n];
    // }
    // inline double operator()(size_t _c, size_t n) const{
    //     return this->f_array[_c][n];
    // }
    
    // // Accesses the continuum, not atom-split parts of it.
    // inline double& operator[](size_t n) {
    //     return this->f_array[0][n];
    // }

    // inline double operator[](size_t n) const{
    //     return this->f_array[0][n];
    // }

    // vector-space algebra

    /**
     * @brief Adds the densities of a distribution to the calling distribution.
     * @param d The Distribution object whose spline factors d.f are added to the caller
     * @return Distribution& 
     */
    Distribution& operator+=(const Distribution& d) {
        for (size_t _c = 0; _c < num_continuums; _c++){  // Iterate through each continuum
            for (size_t i=0; i<size; i++) {
                f_array[_c][i] += d.f_array[_c][i];
            }
        }
        // total += d.total;
        return *this;
    }

    /**
     * @brief Scales electron density by given double
     * @param x Factor to multiply the electron densities at each grid point by
     * @return Distribution& 
     */
    Distribution& operator*=(double x) {
        for (auto& f : this->f_array){
            for (size_t i=0; i<size; i++) {
                f[i] *= x;
            }
        }
        // total *= d.total;
        return *this;
    }

    /**
     * @brief Sets the electron density distribution to equal that of the given Distribution. 
     * @param d Distribution object that the caller sets its electron density to.
     * @return Distribution& 
     */
    Distribution& operator=(const Distribution& d) {
        f_array = d.f_array;
        // for (size_t _c = 0; _c < num_continuums; _c++){ 
        //     f_array[_c] = d.f_array[_c];
        // }
        // total = d.total;
        return *this;
    }

    /**
     * @brief Sets the spline factors at each grid point (energy) to the given double.
     * @param y Value for the electron density to be set to. 
     * @return Distribution& 
     */
    Distribution& operator=(double y) {
        // total = y;
        for (auto& f : this->f_array){
            f.resize(size);
            for (size_t i=0; i<size; i++) {
                f[i] = y;
            }
        }
        return *this;
    }

    static vector<double> get_knot_energies(){return basis.get_knot();}
    static double num_basis_funcs(){return basis.num_funcs;}
    static void initialise_dynamic_regions(DynamicGridPreset& preset){basis.initialise_regions(preset);}


    /**
     * @brief Returns an order-preserved copy of knots with points that overlap with the basis's boundary removed.
     * 
     * @param knots Knot energies ordered from lowest to highest.
     * @return std::vector<double> 
    */
    static std::vector<double> get_trimmed_knots(std::vector<double> knots);

    /**
     * @brief Replaces non-boundary grid points and their associated density spline basis expansion factors with the ones provided. 
     * @details Assumes that these knots are the knots of the simulation and have not changed, in which case it is an alternative to set_basis() that allows for a specific distribution state to be set.
     * @param new_knot new knots (Attention: all knots that are at or below basis._min, or above basis._max, are ignored)
     * @param new_f new density
     * @attention note due to code updates this function is now misleading, as is described in details.
     * @return Distribution& 
     * @todo make clearer this is only for static grid case or just put into one function.
     */
    void set_distribution_STATIC_ONLY(size_t _c, vector<double> new_knot, vector<double> new_f);
    void set_spline_factors(size_t _c, vector<double> new_f);
    // loads the knot and appropriately updates params like the overlap matrix as done in set_parameters.
    static void load_knot(vector<double> loaded_knot);

    double norm(size_t _c) const;
    //double integral(size_t _c) const; //unused
    // modifiers

    /**
     * @brief Updates the grid's spline factors (f) with a given deltaf in the spline basis.
     * @details Adds to f, where f is the electron density vector. 
     * @todo Original comment said v was column vector df/dt, but likely mean the change in f for a time step. 
     * This seems correct, given that dfdt is inputted for the spline basis's Sinv (S_inverse * <VectorXd input>) function, which names the input deltaf.
     * @param v deltaf
     */
    void applyDeltaF(size_t a, const Eigen::VectorXd& dfdt, const int & threads);
    void applyDeltaF_element_scaled(size_t a, const Eigen::VectorXd& dfdt, const int & threads);

    // Q functions
    // Computes the dfdt vector v based on internal f
    // e.g. dfdt v; F.calc_Qee(v);
    void get_Q_eii (size_t _c, Eigen::VectorXd& v, size_t a, const bound_t& P, const int & threads) const;
    void get_Q_tbr (size_t _c, Eigen::VectorXd& v, size_t a, const bound_t& P, const int & threads) const;
    void get_Q_ee  (size_t _c, Eigen::VectorXd& v, const int & threads) const;

    void get_Jac_ee (Eigen::MatrixXd& J) const; // Returns the Jacobian of Qee
    
    /// N is the Number density (inverse au^3) of particles to be added at energy e.
    static void addDeltaLike(Eigen::VectorXd& v, double e, double N);
    /// Adds a Dirac delta to the distribution
    void addDeltaSpike(size_t a, double N, double e);
    /// Applies the loss term to the distribution 
    void addLoss(size_t _c, const Distribution& d, const LossGeometry& l, double charge_density);
    void addFiltration(size_t _c, const Distribution& d, const Distribution& bg,const LossGeometry &l);
    
    /// Sets the object to have a MB distribution
    void add_maxwellian(size_t _c, double N, double T);

    /**
     * @brief 
     * 
     * @param densities Densities corresponding to each energy.
     * @param energies Energies in ascending order.
     */
    void add_density_distribution(size_t _c, std::vector<vector<double>>);

    // Precalculators
    static void Gamma_eii( eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eii, size_t J) {
        return basis.Gamma_eii(Gamma, eii, J);
    }
    static void Gamma_tbr( eiiGraph& Gamma, const std::vector<RateData::InverseEIIdata>& tbr, size_t J, size_t K) {
        return basis.Gamma_tbr(Gamma, tbr, J, K);
    }
    static void precompute_Q_coeffs(vector<RateData::Atom>& Store) {
        #ifndef NO_EII
        basis.precompute_QEII_coeffs(Store);   
        #endif
        #ifndef NO_TBR
        basis.precompute_QTBR_coeffs(Store);
        #endif
        #ifndef NO_EE
        basis.precompute_QEE_coeffs();     
        #endif
    }

    double integral(double (f)(double));
    double density() const;
    double density(size_t cutoff) const;
    /// Returns an estimate of the plasma temperature based on all entries below cutoff (in energy units)
    double k_temperature(size_t cutoff = size) const;
    double CoulombLogarithm() const;


    static std::string output_energies_eV(size_t num_pts);
    static std::vector<double> get_energies_eV(size_t num_pts);
    /**
     * @brief Returns approx. num_pts densities by outputting the same number of points per knot. For a dynamic grid, the energies are determined based on the final energy grid used.
     * @details Always outputs minimum of 1 point per knot.
     * @param num_pts 
     * @return 
     */
    std::string output_densities(size_t _c, size_t num_pts, std::vector<double> reference_knots) const;

    std::vector<double> get_densities(size_t _c, size_t num_pts, std::vector<double> reference_knots) const;
    static std::string output_knots_eV();
    // 
    /**
     * @brief Switch grid points to the knot_energies provided, and replace densities with the ones interpolated to by the splines. 
     * 
     * @param all_knot_energies All knots including basis._max and basis._min. (Currently does not support replacing boundaries). 
     * @note if this is to be used for more than sim loading, need to incorporate running initialise_grid_with_computed_cross_sections again
     */
    void transform_basis(std::vector<double> new_knots);

    // (Unused) This does electron-electron because it is CURSED
    void from_backwards_Euler(double dt, const Distribution& prev_step, double tolerance, unsigned maxiter);

    double operator()(size_t _c, double e) const;

    /**
     * @brief The setup function.
     * @details Grants the distribution its energy basis, which serves as the knot points for the spline. Assigns CoulombLog_cutoff and Distribution::CoulombDens_min.
     * @param n  Num elec points
     * @param min_e min elec energy
     * @param max_e max elec energy
     * @param grid_style 
     */
    static void set_basis(size_t step, GridSpacing grid_style, Cutoffs param_cutoffs, FeatureRegimes regimes, ManualGridBoundaries elec_grid_regions);
    static void set_basis(size_t step, Cutoffs param_cutoffs, FeatureRegimes regimes, std::vector<double> knots,bool update_knot_history = true);
    
    //void prep_adapt_knots(const FeatureRegimes& regimes);
    int adapt_knots(bool init, GridSpacing gt, Distribution original_distribution, std::vector<double> original_knots, std::vector<double> reference_energies = {}, std::vector<double> prev_Edens = {},size_t iteration = 0);

    static double get_min_E(){
        return basis.min_elec_e();
    }
    static double get_lowest_allowed_knot(){
        return basis.first_gp_min_E;
    }

    static double get_max_E(){
        return basis.max_elec_e();
    }
    static double get_mb_min(){
        return basis.mb_min_over_kT;
    }
    static double get_mb_max(){
        return basis.mb_max_over_kT;
    }
    static size_t size;
    static size_t num_continuums;
    static std::vector<double> load_knots_from_history(size_t step_idx);
    static size_t most_recent_knot_change_idx(size_t step_idx);
    static size_t next_knot_change_idx(size_t step_idx);
    static std::vector<double> get_knots_from_history(size_t step_idx);
    static void set_knot_history(size_t i, std::vector<double> replacement_knot){knots_history[i]={i,replacement_knot};}
    // history of grid points (for dynamic grid)
    static std::vector<indexed_knot> knots_history;  
    int container_size(){return f_array[0].size();}  
    
    static bool reset_on_next_grid_update; //TODO part of a duct tape implementation of FIND_INITIAL_DIRAC

private:
    std::vector<std::vector<double>> f_array;  // Spline expansion factors, placed in vector so have option to track a continuum for each element's cascades.
    static SplineIntegral basis;
    static size_t CoulombLog_cutoff;
    /// Coulomb repulsion is ignored if density is below this threshold
    static double CoulombDens_min;
};
ostream& operator<<(ostream& os, const Distribution& dist);

#endif /* end of include guard: RATESYSTEM_CXX_H */
