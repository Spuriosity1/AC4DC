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
#include <eigen3/Eigen/SparseCore>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/SparseCholesky>
#include "Constant.h"
#include "SplineIntegral.h"
#include "Dipole.h"
#include "GridSpacing.hpp"
#include "LossGeometry.hpp"
#include "config.h"


// Represents a statistical distribution of electrons. Internal units are atomic units.
class Distribution
{
public:
    Distribution() {
        f.resize(size);
    }

    inline double& operator[](size_t n) {
        return this->f[n];
    }

    inline double operator[](size_t n) const{
        return this->f[n];
    }

    // vector-space algebra
    Distribution& operator+=(const Distribution& d) {
        for (size_t i=0; i<size; i++) {
            f[i] += d.f[i];
        }
        // total += d.total;
        return *this;
    }

    Distribution& operator*=(double x) {
        for (size_t i=0; i<size; i++) {
            f[i] *= x;
        }
        // total *= d.total;
        return *this;
    }

    Distribution& operator=(const Distribution& d) {
        f = d.f;
        // total = d.total;
        return *this;
    }

    Distribution& operator=(double y) {
        // total = y;
        for (size_t i=0; i<size; i++) {
            f[i] = y;
        }
        return *this;
    }

    double norm() const;
    double integral() const;
    // modifiers

    // applies the df/dt vector v to the overall distribution
    void applyDelta(const Eigen::VectorXd& dfdt);
    

    // Q functions
    // Computes the dfdt vector v based on internal f
    // e.g. dfdt v; F.calc_Qee(v);
    void get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const;
    void get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const;
    void get_Q_ee  (Eigen::VectorXd& v) const;

    void get_Jac_ee (Eigen::MatrixXd& J) const; // Returns the Jacobian of Qee
    
    // N is the Number density (inverse au^3) of particles to be added at energy e.
    static void addDeltaLike(Eigen::VectorXd& v, double e, double N);
    // Adds a Dirac delta to the distribution
    void addDeltaSpike(double N, double e);
    // Applies the loss term to the distribution 
    void addLoss(const Distribution& d, const LossGeometry& l, double charge_density);
    
    // Sets the object to have a MB distribution
    void add_maxwellian(double N, double T);

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
    // Returns an estimate of the plasma temperature based on all entries below cutoff (in energy units)
    double k_temperature(size_t cutoff = size) const;
    double CoulombLogarithm() const;

    static std::string output_energies_eV(size_t num_pts);
    std::string output_densities(size_t num_pts) const;
    
    // This does electron-electron because it is CURSED
    void from_backwards_Euler(double dt, const Distribution& prev_step, double tolerance, unsigned maxiter);

    double operator()(double e) const;

    // The setup function
    static void set_elec_points(size_t n, double min_e, double max_e, GridSpacing grid_style);
    static std::string output_knots_eV();


    static size_t size;
private:
    // double total;
    std::vector<double> f;
    static SplineIntegral basis;
    static size_t CoulombLog_cutoff;
    static double CoulombDens_min; // ignores Coulomb repulsion if sensity is below this threhsold

};

ostream& operator<<(ostream& os, const Distribution& dist);

#endif /* end of include guard: RATESYSTEM_CXX_H */
