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
#include <eigen3/Eigen/StdVector>

// #define NDEBUG
// to remove asserts

// Initialise static things
size_t Distribution::size=0;
size_t Distribution::CoulombLog_cutoff=0;
double Distribution::CoulombDens_min=0;
SplineIntegral Distribution::basis;

// Psuedo-constructor thing
void Distribution::set_elec_points(size_t n, double min_e, double max_e, GridSpacing grid_style) {
    // Defines a grid of n points
    basis.set_parameters(n, min_e, max_e, grid_style);
    Distribution::size=n;
    Distribution::CoulombLog_cutoff = basis.i_from_e(grid_style.transition_e);
    Distribution::CoulombDens_min = grid_style.min_coulomb_density;
    cout<<"[ Free ] Estimating lnLambda based on first ";
    cout<<CoulombLog_cutoff<<" points, up to "<<grid_style.transition_e<<" Ha"<<endl;
    cout<<"[ Free ] Neglecting electron-electron below density of n = "<<CoulombDens_min<<"au^-3"<<endl;
}

// Adds Q_eii to the parent Distribution
void Distribution::get_Q_eii (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(basis.has_Qeii());
    assert(P.size() == basis.Q_EII[a].size());
    assert((unsigned) v.size() == size);
    
    for (size_t xi=0; xi<P.size(); xi++) {
        // Loop over configurations that P refers to
        for (size_t J=0; J<size; J++) {
            for (size_t K=0; K<size; K++) {
                v[J] += P[xi]*f[K]*basis.Q_EII[a][xi][J][K];
            }
        }
    }
}

// Puts the Q_TBR changes in the supplied vector v
void Distribution::get_Q_tbr (Eigen::VectorXd& v, size_t a, const bound_t& P) const {
    assert(basis.has_Qtbr());
    assert(P.size() == basis.Q_TBR[a].size());
    for (size_t eta=0; eta<P.size(); eta++) {
        // Loop over configurations that P refers to
        for (size_t J=0; J<size; J++) {
            for (auto& q : basis.Q_TBR[a][eta][J]) {
                 v[J] += q.val * P[eta] * f[q.K] * f[q.L];
            }
        }
    }
}


// Puts the Q_EE changes into v
void Distribution::get_Q_ee(Eigen::VectorXd& v) const {
    assert(basis.has_Qee());
    double CoulombLog = this->CoulombLogarithm();
    // double CoulombLog = 3.4;
    if (CoulombLog <= 0) return; // skip calculation if LnLambda vanishes

    if (isnan(CoulombLog) || CoulombLog > 11.5) CoulombLog = 11.5;
    // double LnLambdaD = 0.5*log(this->k_temperature()/4/Constant::Pi/this->density());
    // if (isnan(LnLambdaD)) LnLambdaD = 11;
    // cerr<<"LnDebLen = "<<LnLambdaD<<endl;
    // A guess. This should only happen when density is zero, so Debye length is infinity.
    // Guess the sample size is about 10^5 Bohr. This shouldn't ultimately matter much.
    for (size_t J=0; J<size; J++) {
        for (size_t K=0; K<size; K++) {
            for (auto& q : basis.Q_EE[J][K]) {
                 v[J] += q.val * f[K] * f[q.idx] * CoulombLog;
            }
        }
    }
}

void Distribution::get_Jac_ee(Eigen::MatrixXd& M) const{
    // Returns Q^p_qjc^q + Q^p_jrc^r
    assert(basis.has_Qee());
    M = Eigen::MatrixXd::Zero(size,size);
    double CoulombLog=3;
    // double CoulombLog = CoulombLogarithm(size/3);
    if (isnan(CoulombLog) || CoulombLog <= 0) return;
    for (size_t P=0; P<size; P++) {
        for (size_t Q=0; Q<size; Q++) {
            for (auto& q : basis.Q_EE[P][Q]) {
                // q.idx is  L
                M(P,q.idx) += q.val * f[Q] * CoulombLog;
                M(P, Q) += q.val * f[q.idx] * CoulombLog;
            }
        }
    }
}

/*
// Sets f_n+1 based on using a Newton-Rhapson-Euler stiff solver
// Not currently used
void Distribution::from_backwards_Euler(double dt, const Distribution& prev_step_dist, double tolerance, unsigned maxiter){

    Eigen::MatrixXd M(size, size);
    Eigen::MatrixXd A(size, size);
    Eigen::VectorXd v(size);
    Eigen::VectorXd old_curr_step(size);
    
    unsigned idx = 0;
    const unsigned MAX_ITER = 200;

    bool converged = false;

    Eigen::Map<const Eigen::VectorXd > prev_step (prev_step_dist.f.data(), size);
    Eigen::Map<Eigen::VectorXd > curr_step(f.data(), size);

    while (!converged && idx < MAX_ITER){
        
        old_curr_step = curr_step;
        v = Eigen::VectorXd::Zero(size);
        this->get_Q_ee(v, size/3); // calculates v based on curr_step
        
        // Newton iterator
        get_Jac_ee(M);
        A = dt * basis.Sinv(M) - Eigen::MatrixXd::Identity(size,size);
        
        curr_step = -A.fullPivLu().solve(dt*v + prev_step - curr_step);
        curr_step += prev_step;
        
        // Picard iterator
        // curr_step = prev_step + dt*v;

        double delta = fabs((curr_step-old_curr_step).sum());

        if (delta/curr_step.sum() < tolerance) converged = true;
        idx++;
    }
    if (idx == MAX_ITER) std::cerr<<"Max stiff-solver iterations exceeded"<<std::endl;
}
*/

// // Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// void Distribution::add_Q_ee(const Distribution& d, double kT) {
//     double density=0;
//     double CoulombLog = log(kT/(4*Constant::Pi*density));
//     double alpha = 2*Constant::Pi*sqrt(2)/3*CoulombLog;


// }

// ============================
// utility

// Expects T to have units of Ha
void Distribution::add_maxwellian(double T, double N) {
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
        this->f[i] += u[i];
    }
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

std::string Distribution::output_densities(size_t num_pts) const {
    std::stringstream ss;
    const double units = 1./Constant::eV_per_Ha/Constant::Angs_per_au/Constant::Angs_per_au/Constant::Angs_per_au;
    size_t pts_per_knot = num_pts / basis.gridlen();
    if (pts_per_knot == 0){
        pts_per_knot = 1;
        cerr<<"[ warn ] Outputting a small number of points per knot."<<endl;
    }
    for (size_t i=0; i<basis.gridlen()-1; i++){
        double e = basis.grid(i);
        double de = (basis.grid(i+1) - e)/(pts_per_knot-1);
        for (size_t j=0; j<pts_per_knot; j++) {
            ss << (*this)(e)*units<<" ";
        e += de;
        }
    }
    return ss.str();
}

double Distribution::operator()(double e) const{
    double tmp=0;
    for (size_t j = 0; j < size; j++) {
        tmp += basis(j, e)*f[j];
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

void Distribution::addDeltaSpike(double e, double N) {
    int idx = basis.i_from_e(e);
    f[idx] += N/basis.areas[idx];
}


void Distribution::applyDelta(const Eigen::VectorXd& v) {
    Eigen::VectorXd u(size);
    u= (this->basis.Sinv(v));
    for (size_t i=0; i<size; i++) {
        f[i] += u[i];
    }
}


// - 3/sqrt(2) * 3 sqrt(e) * f(e) / R_
// Very rough approximation used here
void Distribution::addLoss(const Distribution& d, const LossGeometry &l, double rho) {
    // f += "|   i|   ||   |_"

    double escape_e = 1.333333333*Constant::Pi*l.L0*l.L0*rho;
    for (size_t i=basis.i_from_e(escape_e); i<size; i++) {
        
        f[i] -= d[i] * sqrt(basis.avg_e[i]/2) * l.factor();
    }
}

void Distribution::addDeltaLike(Eigen::VectorXd& v, double e, double height) {
    assert(e>basis.supp_min(0));
    for (size_t i=0;i<size; i++) {
        v[i] += basis(i, e) * height;
    }

}

double Distribution::norm() const{
    assert(f.size() == Distribution::size);
    double x=0;
    for (auto&fi : f) {
        x+= fabs(fi);
    }
    return x;
}

double Distribution::density() const {
    double tmp=0;
    for (size_t i = 0; i < size; i++)
    {
        tmp += basis.areas[i]*f[i];
    }
    return tmp;
}

double Distribution::density(size_t cutoff) const {
    double tmp=0;
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.areas[i]*f[i];
    }
    return tmp;
}

// Returns an estimate of k*T
// based on kinetic energy density of the plasma below the index given by 'cutoff'
double Distribution::k_temperature(size_t cutoff) const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = this->density(cutoff);
    for (size_t i = 0; i < cutoff; i++)
    {
        tmp += basis.avg_e[i]*f[i]*basis.areas[i];
    }
    return tmp*2./3./n;
}

// Returns an estimate of ln(N_D) based on density and low-energy arguments
double Distribution::CoulombLogarithm() const {
    double tmp=0;
    // Dodgy integral of e * f(e) de
    double n = 0;
    for (size_t i = 0; i < CoulombLog_cutoff; i++) {
        tmp += basis.avg_e[i]*f[i]*basis.areas[i];
        n += f[i]*basis.areas[i];
    }
    double kT = tmp*2./3./n;
    double DebyeLength3 = pow(kT/4/Constant::Pi/n,1.5);
    n = this->density();
    if (n < this->CoulombDens_min) return 0;
    return log(4.*Constant::Pi* n*DebyeLength3);
}


double Distribution::integral(double (g)(double)) {
    double retval = 0;
    for (size_t J = 0; J < Distribution::size; J++)
    {
        double a = basis.supp_min(J);
        double b = basis.supp_max(J);
        double tmp = 0;
        for (size_t i=0;i<10; i++)
        {
            double e = gaussX_10[i]*(b-a)/2 + (b+a)/2;
            tmp += gaussW_10[i]*basis(J, e)*g(e);
        }
        retval +=  tmp * f[J] * (b-a)/2;
    }
    return retval;
}


ostream& operator<<(ostream& os, const Distribution& dist) {
    for (size_t i= 0; i < Distribution::size; i++) {
        os<<" "<<dist[i]/Constant::eV_per_Ha;
    }
    return os;
}
