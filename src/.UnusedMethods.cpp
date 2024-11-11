#include "FreeDistribution.h"
#include "Dipole.h"
#include "Constant.h"
#include "SplineIntegral.h"
#include <Eigen/StdVector>

//////////////////////////////////////////////// UNUSED //////////////////////////////////////////////////////////////
///////Distribution///////
// Sets f_n+1 based on using a Newton-Rhapson-Euler stiff solver
void Distribution::from_backwards_Euler(double dt, const Distribution& prev_step_dist, double tolerance, unsigned maxiter){

    Eigen::MatrixXd M(size, size);
    Eigen::MatrixXd A(size, size);
    Eigen::VectorXd v(size);
    Eigen::VectorXd old_curr_step(size);
    
    unsigned idx = 0;
    const unsigned MAX_ITER = 200;

    bool converged = false;

    Eigen::Map<const Eigen::VectorXd > prev_step (prev_step_dist.f_array[0].data(), size);
    Eigen::Map<Eigen::VectorXd > curr_step(f_array[0].data(), size);

    while (!converged && idx < MAX_ITER){
        
        old_curr_step = curr_step;
        v = Eigen::VectorXd::Zero(size);
        this->get_Q_ee(0, v, size/3); // calculates v based on curr_step
        
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

// // Taken verbatim from Rockwood as quoted by Morgan and Penetrante in ELENDIF
// void Distribution::add_Q_ee(const Distribution& d, double kT) {
//     double density=0;
//     double CoulombLog = log(kT/(4*Constant::Pi*density));
//     double alpha = 2*Constant::Pi*sqrt(2)/3*CoulombLog;


// }

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
                M(P,q.idx) += q.val * f_array[0][Q] * CoulombLog;
                M(P, Q) += q.val * f_array[0][q.idx] * CoulombLog;
            }
        }
    }
}
////////////////////////////////////////////////////////////////////////////////////////////////////