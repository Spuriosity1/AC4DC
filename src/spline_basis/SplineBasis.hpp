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

#ifndef SPLINEBASIS_CXX_H
#define SPLINEBASIS_CXX_H

#include <vector>
#include <eigen3/Eigen/SparseLU>
// #include <eigen3/Eigen/LU>
#include <iostream>
#include "GridSpacing.hpp"

static constexpr bool USING_SQRTE_PREFACTOR = true;

class BasisSet{
public:
    BasisSet() {};
    void set_parameters(size_t nfuncs, double min, double max, const GridSpacing& gt);
    Eigen::VectorXd Sinv(const Eigen::VectorXd& deltaf);
    Eigen::MatrixXd Sinv(const Eigen::MatrixXd& J);

    double raw_bspline(size_t i, double x) const;
    double raw_Dbspline(size_t i, double x) const;

    // Returns the ith basis function evaluated at point x, premultiplied by the Frobenius sqrt(x) value.
    double operator()(size_t i, double x) const;
    // Returns the first derivative of the ith basis function at point x
    double D(size_t i, double x) const;
    double at(size_t i, double x) const;
    inline double supp_max(unsigned i) const{
        return knot[i+BSPLINE_ORDER];
    }
    inline double supp_min(unsigned i) const{
        return knot[i];
    }
    size_t gridlen() const{
        return knot.size();
    }
    double grid(size_t i) const{
        return knot[i];
    }
    std::vector<std::pair<double,double>> knots_between(double bottom, double top) const;
    double min_elec_e() {return _min;};
    double max_elec_e() {return _max;};
    size_t num_funcs;
    const static int BSPLINE_ORDER = 3; // 1 = rectangles, 2=linear, 3=quadratic
    std::vector<double> avg_e;
    std::vector<double> log_avg_e;
    std::vector<double> areas;
    int i_from_e(double e);
    int lower_i_from_e(double e);  
protected:
    // Eigen::PartialPivLU<Eigen::MatrixXd > linsolver;
    Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> >  linsolver;
    std::vector<double> knot;
    double overlap(size_t j, size_t k) const;
    void set_knot(const GridSpacing& gt);
    
    double _min;
    double _max;
};


#endif /* end of include guard: SPLINEBASIS_CXX_H */
