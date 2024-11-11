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

#include "SplineIntegral.h"
#include "config.h"

// Resizes containers and fills them with the appropriate values
void SplineIntegral::precompute_QEII_coeffs(vector<RateData::Atom>& Atoms) {
    std::cout<<"[ Q precalc ] Beginning Q_eii computation...";
    // Size Q_EII appropriately
    Q_EII.resize(Atoms.size());
    
    for (size_t a=0; a<Atoms.size(); a++) {
        if (Atoms[a].bound_free_excluded){
            Q_EII[a].resize(0);
            continue;
        }
        std::cout<<"\n[ Q precalc ] Atom "<<a+1<<"/"<<Atoms.size()<<std::endl;
        Q_EII[a].resize(Atoms[a].num_conf);
        for (size_t eta=0; eta<Atoms[a].num_conf; eta++) {
            Q_EII[a][eta].resize(num_funcs);
            for (auto& QaJ : Q_EII[a][eta]) {
                    QaJ.resize(num_funcs, 0);
            }
        }
        // NOTE: length of EII vector is generally shorter than length of P
        // (usually by 1 or 2, so dense matrices are fine)

        size_t counter=1;
        #pragma omp parallel default(none) shared(a, counter, Atoms, std::cout)
		{
			#pragma omp for schedule(dynamic) nowait
            for (size_t J=0; J<num_funcs; J++) {
                #pragma omp critical
                {
                    std::cout<<"\r[ Q precalc ] "<<counter<<"/"<<num_funcs<<" thread "<<omp_get_thread_num()<<std::flush;
                    counter++;
                }
                for (size_t K=0; K<num_funcs; K++) {
                    for (auto eii : Atoms[a].EIIparams) {
                        Q_EII[a][eii.init][J][K] = calc_Q_eii(eii, J, K);
                    }
                }
            }
        }
    }
    _has_Qeii = true;

    std::cout<<"\n[ Q precalc ] Done."<<std::endl;
}

// Resizes containers and fills them with the appropriate values
void SplineIntegral::precompute_QTBR_coeffs(vector<RateData::Atom>& Atoms) {
    std::cout<<"[ Q precalc ] Beginning Q_tbr computation..."<<std::endl;
    // Size Q_TBR appropriately
    Q_TBR.resize(Atoms.size());
    for (size_t a=0; a<Atoms.size(); a++) {
        if (Atoms[a].bound_free_excluded){
            Q_TBR[a].resize(0);
            continue;
        }
        std::cout<<"\n[ Q precalc ] Atom "<<a+1<<"/"<<Atoms.size()<<std::endl;
        Q_TBR[a].resize(Atoms[a].num_conf);
        for (size_t eta=0; eta<Atoms[a].num_conf; eta++) {
            Q_TBR[a][eta].resize(num_funcs);
            for (auto& QaJ : Q_TBR[a][eta]) {
                    //QaJ.resize(0); // thse are the sparse lists
                    //QaJ = sparse_matrix(0);
                    QaJ.clear(); // these are the sparse lists
            }
        }
        // NOTE: length of EII vector is generally shorter than length of P
        // (usually by 1 or 2, so dense matrices are fine)

        size_t counter=1;
        #pragma omp parallel default(none) shared(a, counter, Atoms, std::cout)
		{
			#pragma omp for schedule(dynamic) nowait
            for (size_t J=0; J<num_funcs; J++) {
                #pragma omp critical
                {
                    std::cout<<"\r[ Q precalc ] "<<counter<<"/"<<num_funcs<<" thread "<<omp_get_thread_num()<<std::flush;
                    counter++;
                }
                
                for (auto& tbr : RateData::inverse(Atoms[a].EIIparams)) {
                    Q_TBR[a][tbr.init][J] = calc_Q_tbr(tbr, J);
                }
            }
        }
    }
    
    _has_Qtbr = true;
    std::cout<<"\n[ Q precalc ] Done."<<std::endl;
}

// Resizes containers and fills them with the appropriate values
void SplineIntegral::precompute_QEE_coeffs() {
    std::cout<<"[ Q precalc ] Beginning Q_ee computation..."<<std::endl;
    Q_EE.resize(num_funcs);
    for(size_t J=0; J<num_funcs; J++) {
        Q_EE[J].resize(num_funcs);
        for (size_t K = 0; K < num_funcs; K++) {
            Q_EE[J][K] = calc_Q_ee(J, K);
        }
    }
    _has_Qee = true;
    std::cout<<"\n[ Q precalc ] Done."<<std::endl;
}

void SplineIntegral::Gamma_eii( std::vector<SparsePair>& Gamma_xi, const RateData::EIIdata& eii, size_t K) const{
    //Gamma_xi.resize(0);
    Gamma_xi.clear();
    for (size_t eta = 0; eta < eii.fin.size(); eta++) {
        // double a = max((float) this->supp_min(K), eii.ionB[eta]);
        double a = this->supp_min(K);
        double b = this->supp_max(K);
        double tmp=0;
        for (auto&& boundary : this->knots_between(a, b)) {
            double diff = (boundary.second - boundary.first)/2.;
            double avg = (boundary.second + boundary.first)/2.;
            double tmp2 = 0;
            for (int i=0; i<GAUSS_ORDER_EII; i++) {
                double e = gaussX_EII[i]*diff + avg;
                tmp2 += gaussW_EII[i]* (*this)(K, e)*pow(e,0.5)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            tmp += tmp2*diff*1.4142; // Electron mass = 1 in atomic units
        }
        if (tmp > DBL_CUTOFF_QEE){
            SparsePair rate;
            rate.idx=eii.fin[eta];
            rate.val= tmp;
            Gamma_xi.push_back(rate);
        }
        
        // double e = knot[K];
        // double tmp = knot[K+1]-knot[K];
        // assert(this->supp_min(K) == knot[K]);
        // assert(this->supp_max(K) == knot[K+1]);
        // tmp *= pow(e,0.5) * Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
        // tmp *= 1.4142135624;

        
    }
}
/*
void SplineIntegral::Gamma_eii( std::vector<SparsePair>& Gamma_xi, const RateData::EIIdata& eii, size_t K) const{
    Gamma_xi.resize(0);
    for (size_t eta = 0; eta < eii.fin.size(); eta++) {
        // double a = max((float) this->supp_min(K), eii.ionB[eta]);
        double a = this->supp_min(K);
        double b = this->supp_max(K);
        double tmp=0;
        if (b > a) {
            for (int i=0; i<GAUSS_ORDER_EII; i++) {
                double e = gaussX_EII[i]*(b-a)/2 + (a+b)/2;
                tmp += gaussW_EII[i]* (*this)(K, e)*pow(e,0.5)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            tmp *= (b-a)/2;
            tmp *= 1.4142; // Electron mass = 1 in atomic units
            SparsePair rate;
            rate.idx=eii.fin[eta];
            rate.val= tmp;
            Gamma_xi.push_back(rate);
        }        
    }
}*/

void SplineIntegral::Gamma_tbr( std::vector<SparsePair>& Gamma_xi, const RateData::InverseEIIdata& tbr, size_t K, size_t L) const {
    // Naming convention is unavoidably confusing.
    // init and fin correspond to initial and final states for the three-body process, i.e. 
    //Gamma_xi.resize(0);
    Gamma_xi.clear();
    Gamma_xi.reserve(tbr.fin.size());
    for (size_t eta = 0; eta<tbr.fin.size(); eta++) {
        double tmp=0;
        double B=tbr.ionB[eta];

        double aK = this->supp_min(K);
        double bK = this->supp_max(K);
        double aL = this->supp_min(L);
        double bL = this->supp_max(L);

        for (int i=0; i<GAUSS_ORDER_TBR; i++) {
            double e = gaussX_TBR[i]*(bK-aK)/2 + (aK+bK)/2;
            double tmp2=0;

            for (int j=0; j<GAUSS_ORDER_TBR; j++) {
                double s = gaussX_TBR[j]*(bL-aL)/2 + (aL+bL)/2;
                double ep = s + e + B;
                tmp2 += gaussW_TBR[j]*(*this)(L, s)*ep*pow(s,-0.5)*
                    Dipole::DsigmaBEB(ep, e, B, tbr.kin[eta], tbr.occ[eta]);
            }
            tmp2 *= (bL-aL)/2;
            tmp += gaussW_TBR[i]*(*this)(K, e)*tmp2*pow(e,-0.5);
        }
        tmp *= (bK-aK)/2;

        tmp *= Constant::Pi*Constant::Pi;

        if (tmp > DBL_CUTOFF_TBR){
            SparsePair rate(tbr.fin[eta], tmp);
            Gamma_xi.push_back(rate);
        }
    }
}

void SplineIntegral::Gamma_eii(eiiGraph& Gamma, const std::vector<RateData::EIIdata>& eiiVec, size_t K) const {
    // 4pi*sqrt(2) \int_0^\inf e^1/2 phi_k(e) sigma(e) de
    //Gamma.resize(0);
    Gamma.clear();
    Gamma.reserve(eiiVec.size());
    for (size_t init=0; init<eiiVec.size(); init++) {
        assert((unsigned) eiiVec[init].init == init);
        std::vector<SparsePair> Gamma_xi(0);
        this->Gamma_eii(Gamma_xi, eiiVec[init], K);
        Gamma.push_back(Gamma_xi);
    }
}

void SplineIntegral::Gamma_tbr(eiiGraph& Gamma, const std::vector<RateData::InverseEIIdata>& eiiVec, size_t K, size_t L) const {
    //Gamma.resize(0);
    //Gamma = eiiGraph(0);
    //eiiGraph tmp;
    //Gamma.swap(tmp);
    Gamma.clear();
    Gamma.reserve(eiiVec.size());
    for (size_t init=0; init<eiiVec.size(); init++) {
        assert((unsigned) eiiVec[init].init == init);
        std::vector<SparsePair> Gamma_xi(0);
        this->Gamma_tbr(Gamma_xi, eiiVec[init], K, L);
        Gamma.push_back(Gamma_xi);
    }
}

/*
double SplineIntegral::calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K) const {
    // computes sum of all Q_eii integrals away from eii.init, integrated against basis functions f_J, f_K
    // J is the 'free' index
    // sqrt(2) sum_eta \int dep sqrt(ep) f(ep) dsigma/de (e | ep) - sqrt(e) sigma * f_J(e)
    // 'missing' factor of 2 in front of RHS is not a mistake! (arises from definition of dsigma)

    // AD HOC TIME
    // if ( J > num_funcs - BSPLINE_ORDER) return 0;
    // if ( K > num_funcs - BSPLINE_ORDER) return 0;


    double min_J = this->supp_min(J);
    double max_J = this->supp_max(J);

    double retval=0;
    // Sum over all transitions to eta values
    for (size_t eta = 0; eta<eii.fin.size(); eta++)
    {
        for (auto&& boundary : this->knots_between(min_J, max_J)) {
            double diff = (boundary.second - boundary.first)/2.;
            double avg = (boundary.second + boundary.first)/2.;
            double tmp=0;
            for (int j=0; j<GAUSS_ORDER_EII; j++)
            {
                double e = gaussX_EII[j]*diff + avg;
                double min_ep = max(e+eii.ionB[eta], this->supp_min(K));
                double max_ep = this->supp_max(K);
                // for (auto&& boundary2 : this->knots_between(min_ep, max_ep)) {
                double tmp2=0;
                // double diff2 = (boundary2.second - boundary2.first)/2.;
                // double avg2 = (boundary2.second + boundary2.first)/2.;
                double diff2 = (max_ep - min_ep) *0.5;
                double avg2 = (max_ep + min_ep) *0.5;
                for (int k=0; k<GAUSS_ORDER_EII; k++) {
                    double ep = gaussX_EII[k]*diff2 + avg2;
                    tmp2 += gaussW_EII[k]*(*this)(K, ep)*pow(ep,0.5)*
                        Dipole::DsigmaBEB(ep, e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
                }
                tmp += gaussW_EII[j]*raw_bspline(J, e)*tmp2 * diff2;
                // }
            }
            retval += tmp * diff;
        }
        // RHS part, 'missing' factor of 2 is incorporated into definition of sigma
        
        double min_JK = max(this->supp_min(J), this->supp_min(K));
        min_JK = max(min_JK, (double) eii.ionB[eta]);
        double max_JK = min(this->supp_max(J), this->supp_max(K));
        for (auto&& boundary : this->knots_between(min_JK, max_JK)) {
            double diff = (boundary.second - boundary.first)/2.;
            double avg = (boundary.second + boundary.first)/2.;
            double tmp=0;
            for (int k=0; k<GAUSS_ORDER_EII; k++)
            {
                double e = gaussX_EII[k]*diff + avg;
                tmp += gaussW_EII[k]*pow(e, 0.5)*raw_bspline(J, e)*(*this)(K, e)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            retval -= tmp*diff;
        }
        
    }

    retval *= 1.4142135624; 

    return retval;
} */


double SplineIntegral::calc_Q_eii( const RateData::EIIdata& eii, size_t J, size_t K) const {
    // computes sum of all Q_eii integrals away from eii.init, integrated against basis functions f_J, f_K
    // J is the 'free' index
    // sqrt(2) sum_eta \int dep sqrt(ep) f(ep) dsigma/de (e | ep) - sqrt(e) sigma * f_J(e)
    // 'missing' factor of 2 in front of RHS is not a mistake! (arises from definition of dsigma)

    double min_J = this->supp_min(J);
    double max_J = this->supp_max(J);

    double retval=0;
    // Sum over all transitions to eta values
    for (size_t eta = 0; eta<eii.fin.size(); eta++)
    {
        double tmp = 0;
        for (int j=0; j<GAUSS_ORDER_EII; j++)
        {
            double e = gaussX_EII[j]*(max_J-min_J)/2 + (max_J+min_J)/2;
            double tmp2=0;
            double min_ep = max(e+eii.ionB[eta], this->supp_min(K));
            double max_ep = this->supp_max(K);
            if (max_ep <= min_ep) continue;
            for (int k=0; k<GAUSS_ORDER_EII; k++) {
                double ep = gaussX_EII[k]*(max_ep-min_ep)*0.5 + (min_ep+max_ep)*0.5;
                tmp2 += gaussW_EII[k]*(*this)(K, ep)*pow(ep,0.5)*
                    Dipole::DsigmaBEB(ep, e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            }
            tmp2 *= (max_ep - min_ep)/2;
            // tmp += gaussW_EII[j]*(*this)(J, e)*tmp2;
            tmp += gaussW_EII[j]*raw_bspline(J, e)*tmp2;
        }
        tmp *= (max_J-min_J)/2;
        retval += tmp;
        // RHS part, 'missing' factor of 2 is incorporated into definition of sigma
        tmp=0;
        double min_JK = max(this->supp_min(J), this->supp_min(K));
        min_JK = max(min_JK, (double) eii.ionB[eta]);
        double max_JK = min(this->supp_max(J), this->supp_max(K));
        if (max_JK <= min_JK) continue;
        for (int k=0; k<GAUSS_ORDER_EII; k++)
        {
            double e = gaussX_EII[k]*(max_JK-min_JK)/2 + (min_JK+max_JK)/2;
            // tmp += gaussW_EII[k]*pow(e, 0.5)*(*this)(J, e)*(*this)(K, e)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
            tmp += gaussW_EII[k]*pow(e, 0.5)*raw_bspline(J, e)*(*this)(K, e)*Dipole::sigmaBEB(e, eii.ionB[eta], eii.kin[eta], eii.occ[eta]);
        }
        tmp *= (max_JK-min_JK)/2;
        retval -= tmp;
    }

    retval *= 1.4142135624; 

    return retval;
}

// JKL indices, but the Q_TBR coeffs are sparse btw J and L
// (overlap by at most BasisSet::BSPLINE_ORDER either way)
// Returns a vector of pairs (Q-idx, nonzzero-idx)
SplineIntegral::sparse_matrix SplineIntegral::calc_Q_tbr( const RateData::InverseEIIdata& tbr, size_t J) const {
    // J is the 'free' index, K and L label the given indices of the free distribution
    //
    double J_min = this->supp_min(J);
    double J_max = this->supp_max(J);

    size_t num_nonzero=0;

    SplineIntegral::sparse_matrix retval; // a vector of SparseTriple
    SparseTriple tup;
    for (size_t K=0; K<num_funcs; K++) {
        for (size_t L=0; L<num_funcs; L++) {
            // First part of the integral. Represents the source term corresponding to
            // the production of an energetic electron of energy ep + s - B
            double K_min = this->supp_min(K);
            double K_max = this->supp_max(K);
            double L_min = this->supp_min(L);
            double L_max = this->supp_max(L);


            double tmp=0;

            // Sum over all transitions in tbr
            for (size_t xi=0; xi<tbr.fin.size(); xi++) {
                double B =tbr.ionB[xi];
                // double aJ_eff = (aJ < B) ? B : aJ;
                // First half of integral:
                // 0.5 * \int dep \int ds B_J(e) B_K(ep) B_L(s) dsigma(ep | s + ep + B) (s+ep+B)/sqrt(s ep)
                // where e = s+B+ep
                // Integration ranges:
                // K_min < ep < K_max
                // L_min < s  < L_max
                // J_min < e  < J_max  ===> J_min - ep - B < s < J_max - ep - B
                // where the ep integral is done first
                for (int j=0; j<GAUSS_ORDER_TBR; j++) {
                    double ep = gaussX_TBR[j]*(K_max- K_min)/2 + (K_max + K_min)/2;
                
                    double a = max(L_min, J_min - ep - B); // integral lower bound
                    double b = min(L_max, J_max - ep - B); // integral upper bound
                    if (a >= b) continue;
                    double tmp2=0;
                    for (int k=0; k<GAUSS_ORDER_TBR; k++) {
                        double s = gaussX_TBR[k]*(b-a)/2 + (b + a)/2;
                        double e = s + ep + B;
                        // tmp2 += gaussW_TBR[k] * ( e ) * pow(s*ep,-0.5) * Dipole::DsigmaBEB(e, ep, B, tbr.kin[xi], tbr.occ[xi])*(*this)(J, e)*(*this)(L, s);
                        tmp2 += gaussW_TBR[k] * ( e ) * pow(s*ep,-0.5) * Dipole::DsigmaBEB(e, ep, B, tbr.kin[xi], tbr.occ[xi])*raw_bspline(J, e)*(*this)(L, s);
                    }
                    tmp += 0.5*gaussW_TBR[j]*tmp2*(*this)(K,ep)*(b-a)*(K_max - K_min)/4;
                }
            }


            // Second part of the integral
            // Represents the sink terms corresponding to the absorption of low energy electrons at this point
            double min_JK = max(J_min, K_min);
            double max_JK = min(J_max, K_max);
            if (max_JK > min_JK ) {
                for (size_t xi=0; xi<tbr.fin.size(); xi++) {
                    double B =tbr.ionB[xi];
                    for (int j=0; j<GAUSS_ORDER_TBR; j++) {
                        double e = gaussX_TBR[j]*(max_JK-min_JK)/2 + (max_JK+min_JK)/2;
                        double tmp2=0;
                        for (int k=0; k<GAUSS_ORDER_TBR; k++) {
                            double s = gaussX_TBR[k]*(L_max-L_min)/2 + (L_max+L_min)/2;
                            tmp2 += gaussW_TBR[k]*(s+e+B)*pow(e * s,-0.5)*(*this)(L,s)*Dipole::DsigmaBEB(s+e+B, e, B, tbr.kin[xi], tbr.occ[xi]);
                        }
                        // tmp -= gaussW_TBR[j]*tmp2*(*this)(J,e)*(*this)(K,e)*(L_max-L_min)*(max_JK-min_JK)/4;
                        tmp -= gaussW_TBR[j]*tmp2*raw_bspline(J,e)*(*this)(K,e)*(L_max-L_min)*(max_JK-min_JK)/4;
                    }
                }
            }

            tmp *= 2*Constant::Pi*Constant::Pi;


            if (fabs(tmp) > DBL_CUTOFF_TBR) {
                tup.K = K;
                tup.L = L;
                tup.val=tmp;
                retval.push_back(tup);
                num_nonzero++;
            }
        }
    }
    return retval;
}


inline double SplineIntegral::Q_ee_F(double e, size_t K) const {
    double K_min = this->supp_min(K);
    double K_max = this->supp_max(K);

    // 2e^-1/2  int_0^e dx x phiK(x) 
    double max_K_0E = min(K_max, e);
    double min_K_0E = K_min;
    double tmp2 = 0;
    for (auto&& boundary : this->knots_between(min_K_0E, max_K_0E)){
        if ( boundary.first < DBL_CUTOFF_QEE) {
            double integ = 0;
            for (size_t j = 0; j < GAUSS_ORDER_SQRT; j++) {
                double x = boundary.second*(gaussX_sqrt[j] + 1)/2.;
                integ += gaussW_sqrt[j]*x*raw_bspline(K, x);
                static_assert(USING_SQRTE_PREFACTOR == true);
            }
            tmp2 += integ*pow(boundary.second/2.,1.5);
        } else {
            double integ = 0;
            double diff = (boundary.second - boundary.first)/2.;
            double avg = (boundary.second + boundary.first)/2.;
            for (size_t j = 0; j < GAUSS_ORDER_EE; j++) {
                double x = gaussX_EE[j] * diff + avg;
                integ += gaussW_EE[j]*x*(*this)(K, x);
            }
            tmp2 += integ*diff;
        }
    }
    tmp2 *= 2*pow(e,-0.5);

    // 2e int_e^inf x^-1/2 phiK(x)
    double max_I2 = K_max;
    double min_I2 = max(K_min, e);
    double tmp3 = 0;
    for (auto&& boundary : this->knots_between(min_I2, max_I2)){
        double integ = 0;
        double diff = (boundary.second - boundary.first)/2.;
        double avg = (boundary.second + boundary.first)/2.;
        for (size_t j = 0; j < GAUSS_ORDER_EE; j++) {
            double x = gaussX_EE[j] * diff + avg;
            integ += gaussW_EE[j]*pow(x,-0.5)*(*this)(K, x);
        }
        tmp3 += integ*diff;
    }
    tmp3 *= 2*e;

    return tmp2 + tmp3;
}    

inline double SplineIntegral::Q_ee_G(double e, size_t K) const {   
    // Begin G integral
    double max_K_0E = min(this->supp_max(K), e);
    double min_K_0E = this->supp_min(K);

    double Gint = 0;
    for (auto&& boundary : this->knots_between(min_K_0E, max_K_0E)) {
        double tmp=0;
        if ( boundary.first < DBL_CUTOFF_QEE) {
            // At the origin, special care must be taken with the sqrt singularity
            for (size_t j = 0; j < GAUSS_ORDER_SQRT; j++) {
                double x = boundary.second*(gaussX_sqrt[j] + 1)/2.;
                tmp += gaussW_sqrt[j]*raw_bspline(K, x);
                static_assert(USING_SQRTE_PREFACTOR == true);
            }
            tmp *= pow(boundary.second/2.,1.5);
            
        } else {
            double diff = (boundary.second - boundary.first)/2.;
            double avg = (boundary.second + boundary.first)/2.;
            for (size_t j = 0; j < GAUSS_ORDER_EE; j++) {
                double x = gaussX_EE[j] * diff + avg;
                tmp += gaussW_EE[j]*(*this)(K, x);
            }
            tmp *= diff;
        }
        Gint += tmp;
    }
    Gint *= 3*pow(e,-0.5);
    // G integral done.

    return  Gint;
}

// Uses a modified form of the result quoted by Rockwood 1973, 
// "Elastic and Inelastic Cross Sections for Electron-Hg Scattering From Hg Transport Data"
// Originally attributable to Rosenbluth
// N.B. this is NOT an efficient implementation, but it only runs once -  \_( '_')_/    <<*Cries in trying to load simulation state* - S.P.>>
// J, K, L refer to spline indices.
SplineIntegral::pair_list SplineIntegral::calc_Q_ee(size_t J, size_t K) const {
    double J_min = this->supp_min(J);
    double J_max = this->supp_max(J);
    pair_list p(0);
    for (size_t L = 0; L < num_funcs; L++) {
        double total=0;
        // F component
        double min_JL = max(J_min, this->supp_min(L));
        double max_JL = min(J_max, this->supp_max(L));

        if (max_JL <= min_JL) continue;

        
        for (auto&& boundary : this->knots_between(min_JL,max_JL)){
            double tmp=0;
            if (boundary.first < DBL_CUTOFF_QEE) {
                for (size_t i = 0; i < GAUSS_ORDER_SQRT; i++) {
                    double e = boundary.second*(gaussX_sqrt[i] + 1)/2.;
                    tmp += gaussW_sqrt[i]*raw_Dbspline(J, e)*(Q_ee_F(e, K) * ( -this->raw_Dbspline(L,e) ) - Q_ee_G(e, K) * this->raw_Dbspline(L,e));
                }
                tmp *= pow(boundary.second/2, 0.5);
                static_assert(USING_SQRTE_PREFACTOR);
            } else {
                double diff = (boundary.second - boundary.first)/2.;
                double avg = (boundary.second + boundary.first)/2.;
                for (size_t i = 0; i < GAUSS_ORDER_EE; i++) {
                    double e = gaussX_EE[i] * diff + avg;
                    // tmp += gaussW_EE[i]*raw_Dbspline(J, e)*(Q_ee_F(e, K) * ( (*this)(L, e)/2/e - this->D(L, e) ) - Q_ee_G(e, K) * (*this)(L,e));
                    static_assert(USING_SQRTE_PREFACTOR);
                    tmp += gaussW_EE[i]*raw_Dbspline(J, e)*(Q_ee_F(e, K) * ( -pow(e,0.5)*this->raw_Dbspline(L,e) ) - Q_ee_G(e, K) * (*this)(L,e));
                }
                tmp *= diff;
            }
            total += tmp;
            
        }
        
        // Multiply by alpha
        total *= 2./3.*Constant::Pi*sqrt(2);
        // N.B. this value still needs to be multiplied by the Coulomb logarithm.
        if (fabs(total) > DBL_CUTOFF_QEE){
            SparsePair s;
            s.idx = L;
            s.val = total;
            p.push_back(s);
        }
    }
    return p;
}
