/** @file HybridIntegrator.hpp
 * @authors Alaric Sanders & Spencer Passmore
 * @brief Defines the Hybrid class which adds a Moulton step after each step from the inherited method of Adams_BM.
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

#ifndef HYBRID_INTEGRATE_HPP
#define HYBRID_INTEGRATE_HPP

#include "AdamsIntegrator.hpp"
#include "Constant.h"
#include "FreeDistribution.h"
#include "RateSystem.h"
#ifdef NCURSES
#include "Display.h"
#endif
#include "GridSpacing.hpp" // for checkpoint
#include <Eigen/Dense>
#include <iostream>

struct Checkpoint {
    size_t n;
    std::vector<double> knots;
    FeatureRegimes regimes;
    std::vector<state_type> last_few_states;
    std::vector<double> last_few_times;
};

namespace ode {
template<typename T>
class Hybrid : public Adams_BM<T>{
    public:
    Hybrid<T>(unsigned int order=4, double rtol=1e-5, unsigned max_iter = 50): 
        Adams_BM<T>(order), stiff_rtol(rtol), stiff_max_iter(max_iter){};
    const static unsigned MAX_BEULER_ITER = 50;  // unused?
    FeatureRegimes regimes;
    double timestep_reached = 0;       
    private:
    virtual void sys_ee(const T& q, T& qdot) =0;
    // virtual void Jacobian2(const T& q, T& qdot, double t) =0; 
    protected:

    double stiff_rtol = 1e-4;  // relative tolerance
    unsigned stiff_max_iter = 300;//200; 
    double intolerable_stiff_divergence =0;//0.5; // allowed (relative) divergence  if stiff_max_iter exceeded
    
    // stiff ode intermediate steps (i.e. steps it does without the nonstiff part)
    int mini_n;
    int old_mini_n;
    #if defined NO_EE || defined NO_MINISTEPS
    int num_stiff_ministeps = 1;
    #elif defined NO_MINISTEP_UPDATING
    int num_stiff_ministeps = 500;
    #else
    int num_stiff_ministeps = 50;
    #endif
    double mini_dt;

    // Grid/timestep dynamics    
    size_t steps_per_grid_transform;    
    Checkpoint old_checkpoint; // the checkpoint that is loaded.
    Checkpoint checkpoint; 
    int checkpoint_gap = 100; 
    size_t increasing_dt = 0;
    // TODO need to separate good_state from euler_exceeded, currently we have the case of good_state is false and 
    // euler is exceeded, which means euler exceeded. And we have the case of good state is false and euler 
    // exceeded is false, which means a nan was encountered.    
    bool good_state = true;
    bool euler_exceeded = false;  // Whether exceeded the max number of euler iterations (on prev. step)    
    // [not in use as bugged] For gaussian pulses, this stores the times that we decreased dt so that it can be increased later.
    //TODO this should be saved if want to load from a point in time, but really it's meant to replace loading.
    std::vector<double> times_to_increase_dt;  
    
    void run_steps(ofstream& _log, const double t_resume, const int steps_per_time_update);  // TODO clean up bootstrapping. -S.P.
    void iterate(ofstream& _log, double t_initial, size_t npoints_initial, const double t_resume, const int steps_per_time_update);
    void initialise_transient_y(int n); // Approximates initial stiff intermediate steps
    void initialise_transient_y_v2(int n); // Uses lagrange interpolation to approximate initial stiff intermediate steps
    void modify_ministeps(const int n,const int num_ministeps);
    #ifdef NO_MINISTEP_UPDATING
    int max_ministep_reductions = 0;
    #else
    int max_ministep_reductions = 5;
    #endif
    std::vector<T> lagrange_interpolate(const std::vector<double> x, const std::vector<T> f, const std::vector<double> new_x);
    std::vector<T> y_transient; // stores transient/intermediate steps for stiff solver
    std::vector<T> old_y_transient; // stores final transient/intermediate steps of last step. 
    // More virtual funcs defined by ElectronRateSolver:
    virtual state_type get_ground_state()=0;
    virtual void pre_ode_step(ofstream& _log, size_t& n,const int steps_per_time_update)=0;
    virtual int post_ode_step(ofstream& _log, size_t& n)=0;
    /// Unused
    void backward_Euler(unsigned n); 
    void step_stiff_part(unsigned n);
};

template<typename T>
// t_resume = the time to resume simulation from if loading a sim. -S.P.
// _log only used for cross-section recalcs atm.
/**
 * @brief 
 * 
 * @param _log 
 * @param t_initial 
 * @param npoints Initial number of time steps, defining the end time along with this->dt. Thus allows flexibility in chosen end time. 
 * @param t_resume 
 * @param steps_per_time_update 
 */
void Hybrid<T>::iterate(ofstream& _log, double t_initial, size_t npoints, const double t_resume, const int steps_per_time_update) {
    if (this->dt < 1E-16) {
        std::cerr<<"WARN: step size "<<this->dt<<"is smaller than machine precision"<<std::endl;
    } else if (this->dt < 0) {
        throw std::runtime_error("Step size is negative!");
    }
    
    bool resume_sim = (t_resume == t_initial) ? false : true;

    size_t resume_n = 0;
    if (resume_sim){
        for (size_t n=1; n<npoints; n++){
            if (this->t[n] >= t_resume){
                resume_n = n;
                break;
            }
        }
        // The time step size does not depend on previous run's time step size. i.e. step size is same as if there was no loading.
        // TODO implement assertion that density isn't empty.
        size_t resume_n_if_const_dt = (t_resume-t_initial)/this->dt ;
        npoints -= (resume_n_if_const_dt + 1);
        npoints += this->t.size(); // 
        // Try to set checkpoint to be at the starting step 
        // Check if need to set it before the last knot change.
        size_t checkpoint_n = resume_n;
        while (Distribution::most_recent_knot_change_idx(checkpoint_n) >= resume_n - this->order + 1){ 
            checkpoint_n--;
            assert(checkpoint_n > this->order);
            assert(checkpoint_n > resume_n - this->order); // may trigger if knot changes too close together.
        }        
        std::vector<state_type> check_states;
        std::vector<double> check_times;
        {
            std::vector<state_type>::const_iterator start_vect_idx = this->y.begin() - this->order + 1 + checkpoint_n;  
            std::vector<state_type>::const_iterator end_vect_idx = this->y.begin() + checkpoint_n;  
            check_states = std::vector<state_type>(start_vect_idx, end_vect_idx+1);
        }
        {
            std::vector<double>::const_iterator start_vect_idx = this->t.begin() - this->order + 1 + checkpoint_n;  
            std::vector<double>::const_iterator end_vect_idx = this->t.begin() + checkpoint_n;  
            check_times = std::vector<double>(start_vect_idx, end_vect_idx+1);   
        }
        assert(check_states.size() == this->order);
        assert(check_times.size() == this->order);
        assert(check_states.front().F.container_size() == check_states.back().F.container_size());

        checkpoint = {checkpoint_n, Distribution::get_knot_energies(),this->regimes,check_states,check_times};
    }
    old_checkpoint = checkpoint; 

    npoints = (npoints >= this->order) ? npoints : this->order;

    // Set up the containers
    this->t.resize(npoints,INFINITY);
    this->y.resize(npoints);

    // Set up the t grid       
    this->t[0] = t_initial;

    for (size_t n=1; n<npoints; n++){
        if (resume_sim && n <= resume_n){
            continue; // Don't reset already simulated states
        }
        this->t[n] = this->t[n-1] + this->dt;
    }
    this->run_steps(_log,t_resume, steps_per_time_update);
}


// Overrides the underlying Adams method, adding a more refined but computationally expensive treatment for the stiff Q^{ee} contribution to deltaf.
/**
 * @brief  Iterates through the (t_final-t_initial = simulation_timespan)/dt timesteps.
 * .s
 * @tparam T 
 */
template<typename T>
void Hybrid<T>::run_steps(ofstream& _log, const double t_resume, const int steps_per_time_update){
    assert(this->y.size() == this->t.size());
    assert(this->t.size() >= this->order);

    size_t n = 0; // Current step index.
    if (t_resume < this->t[this->order]){
        assert(t_resume == this->t[0]);
        // 4th Order Runge-Kutta 
        // initialise enough points for multistepping to get going
        while(n < this->order) {
            this->step_rk4(n);
            n++;
        }
        // Set first checkpoint. 
        std::vector<state_type> check_states;
        std::vector<double> check_times;
        {
            std::vector<state_type>::const_iterator start_vect_idx = this->y.begin() + 1;  
            std::vector<state_type>::const_iterator end_vect_idx = this->y.begin() + this->order;  
            check_states = std::vector<state_type>(start_vect_idx, end_vect_idx+1);
        }
        {
            std::vector<double>::const_iterator start_vect_idx = this->t.begin() + 1;  
            std::vector<double>::const_iterator end_vect_idx = this->t.begin() + this->order;  
            check_times = std::vector<double>(start_vect_idx, end_vect_idx+1);
        }        
        assert(check_states.size() == this->order);
        assert(check_times.size() == this->order);
        assert(check_states.front().F.container_size() == check_states.back().F.container_size());
        checkpoint = {this->order, Distribution::get_knot_energies(), this->regimes, check_states,check_times};
        old_checkpoint = checkpoint;
    }
    else{
        // Loaded simulation
        n+=this->order;
    }

    // Start with n = last step with data.
    while(this->t[n] < t_resume)
        n++;
    initialise_transient_y((int)n);

    // Set up display
    std::stringstream tol;
    tol << "[ sim ] Implicit solver uses relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<"\n\r";
    std::cout << tol.str();  // Display in regular terminal even after ncurses screen is gone.
    #ifdef NCURSES
    Display::header += tol.str(); 
    Display::display_stream = std::stringstream(Display::header, ios_base::app | ios_base::out); // text shown that updates with frequency 1/steps_per_time_update.
    Display::popup_stream = std::stringstream(std::string(), ios_base::app | ios_base::out);  // text displayed during step of special events like grid updates
    Display::create_screen(); 
    #endif

    // Run hybrid multistepping (nonstiff -> bound dynamics, stiff -> free dynamics). 
    while (n < this->t.size()-1) {
        // 1. Display general information
        // 2. Handle live plotting
        // 3. Handle dynamic updates to dt (and checkpoints, which also handle case of NaN being encountered by solver).
        this->pre_ode_step(_log, n,steps_per_time_update);
        this->step_nonstiff_part(n); 
        #ifdef DEBUG_BOUND
        for(size_t a = 0; a < this->y[n+1].atomP.size();a++)
            for(size_t i=0;i < this->y[n+1].atomP[a].size();i++){
                assert(this->y[n+1].atomP[a][i] >= 0);
            }        
        #endif    
        // this->y[n+1].from_backwards_Euler(this->dt, this->y[n], stiff_rtol, stiff_max_iter);
        this->step_stiff_part(n);

        // 1. Handle dynamic grid updates 
        // 2. Check if user wants to exit simulation early  
        if (this->post_ode_step(_log,n) == 1){
            break; // user asked to end the simulation early.
        } 
        #ifdef DEBUG_BOUND
        for(size_t a = 0; a < this->y[n+1].atomP.size();a++)
            for(size_t i=0;i < this->y[n+1].atomP[a].size();i++){
                assert(this->y[n+1].atomP[a][i] >= 0);
            }        
        #endif   
        n++;     
    }
#ifdef NCURSES
    Display::close();
#endif // NCURSES
    std::cout<<"\n";
    // std::cout <<"[ sim ] Implicit solver used relative tolerance "<<stiff_rtol<<", max iterations "<<stiff_max_iter<<"\n";
    std::cout<<"[ sim ] final t = "<<this->t.back() * Constant::fs_per_au<<" fs"<< endl;  
}
/**
 * @brief Use a true implicit method to estimate change to bad part of system based solely on its own action.
 * 
 * @param n n+1 is the index of the step to predict.
 * @return template<typename T> 
 */
template<typename T>
void Hybrid<T>::step_stiff_part(unsigned n){
    if(!good_state) return;
    
    // Loop explanation:
    // mini_n == 'ministep' index (that was *last* updated)
    // y_transient contains the ministeps needed for next step's calculation.
    // y_transient = [1,2,3,4] for cubic splines  (order elements)
    // We use the order=3 most recent elements, and replace the oldest element with the newest
    // Hence we are targeting (mini_n+1)%(order) in Adams-Moulton step
    T delta_bound;
    // 'rel idx' = relative index within transient y
    int last_rel_idx = mini_n%(this->order);
    int next_rel_idx; 

    
    delta_bound = this->y[n];
    delta_bound *= -1.;
    delta_bound += this->y[n+1];    
    // Interpolate bound contribution to y dot so that we can approximate the bound contribution at each intermediate step.
    std::vector<T> delta_bound_interpolated(num_stiff_ministeps);
    for (int i = 0; i < num_stiff_ministeps;i++){
        T tmp = delta_bound;
        tmp *= ((i+1)/num_stiff_ministeps); 
        delta_bound_interpolated[i] += tmp;
    }

    #ifdef DEBUG
    assert(this->y[n].F[3] == y_transient[last_rel_idx].F[3]);
    #endif

    int excess_count = 0;
    int under_count = 0;
    int num_ministep_reductions = 0;
    old_mini_n = mini_n;
    old_y_transient = y_transient; // Stores final ministeps of step n-1.
    while (mini_n < old_mini_n + num_stiff_ministeps){  // mini_n = n if num_stiff_ministeps = 1.
        
        // Adams-Moulton step I believe -S.P.
        T tmp;
        tmp = this->zero_y;
        // tmp acts as an aggregator
        for (int i = 1; i < this->order; i++){  // work through last N=order-1 ministeps. i.e. Order = 3 corresponds to 2 step method.
            T ydot; // ydot stores the change this loop.
            this->sys_ee(y_transient[(1-i+mini_n)%(this->order)], ydot); 
            ydot *= this->b_AM[i];
            tmp += ydot;
        }
        //this->y[n+1] = this->y[n] + tmp * dt;
        tmp *= mini_dt;
        tmp += y_transient[last_rel_idx];

        // tmp now stores the additive constant:  y_n+1 = tmp + h b0 f(y_n+1, t_n+1)

        // Picard iteration
        next_rel_idx = (mini_n+1)%(this->order); 
        T prev;
        double diff = stiff_rtol*2;
        unsigned idx=0;
        while (diff > stiff_rtol/num_stiff_ministeps && idx < stiff_max_iter){ //stiff_rtol is for the full step, so ministeps have smaller tolerance. 
            // solve by Picard iteration
            prev = y_transient[next_rel_idx];
            prev *= -1;
            T dydt;
            this->sys_ee(y_transient[next_rel_idx], dydt);
            dydt *= this->b_AM[0]*(mini_dt);
            y_transient[next_rel_idx] = tmp;
            y_transient[next_rel_idx] += dydt;
            prev += y_transient[next_rel_idx];
            diff = prev.norm()/y_transient[next_rel_idx].norm(); // Seeking convergence
            idx++;
        }
        if(idx==stiff_max_iter){
            if (num_ministep_reductions < max_ministep_reductions){
                // Try again with more ministeps
                modify_ministeps(n,min(num_stiff_ministeps*2,1000));
                mini_n = old_mini_n;
                num_ministep_reductions++;
            }
            //std::cerr<<"Max Euler iterations exceeded, err = "<<diff<<std::endl;
            else if (diff > intolerable_stiff_divergence){
                //std::cerr << "Max error ("<<intolerable_stiff_divergence<<") exceeded, ending simulation early." <<std::endl; // moved to 
                this->good_state = false;
                this->timestep_reached = this->t[n+1]*Constant::fs_per_au; // t[n+1] is equiv. to t in bound !good_state case, where error condition this is modelled off is found.
                this->euler_exceeded = true;  
                break;
            }
        }
        y_transient[next_rel_idx] += delta_bound_interpolated[mini_n-old_mini_n];  // Add interpolated bound state contribution
        #ifndef NO_MINISTEP_UPDATING
        else if (idx >= stiff_max_iter/4){
            excess_count++;
        }
        else if (idx == 1){
            under_count++;
        }
        #endif 
        mini_n++;
        last_rel_idx = mini_n%(this->order);
    }
    #ifndef NO_MINISTEP_UPDATING
    if (excess_count > under_count)
        modify_ministeps(n,min((int)(num_stiff_ministeps*1.5),1000));
    else if (excess_count*2 < under_count)
        modify_ministeps(n,max((int)(this->order-1),(int)(num_stiff_ministeps*0.9)));
    #endif
    this->y[n+1] = y_transient[next_rel_idx];
    
}

/// Initialises intermediate steps needed for stiff solver to get going.
// For order 3, transient y is indexed as: [mini_n-3,mini_n-2,mini_n-1,mini_n]
template<typename T>
void Hybrid<T>::initialise_transient_y(int n) {
    assert(this->y[n-1].F.container_size() == this->y[n].F.container_size());

    y_transient.resize(this->order);
    mini_n = this->order-1; // last index of simulation's first "loop" of transient y. 
    if ((this->order-1)/num_stiff_ministeps <= 1){
        mini_dt = this->dt/(double)num_stiff_ministeps; 
        // Approximate intermediate steps with by interpolating between n and n - 1.
        // We store y[n] and flag it as the most recent ministep. 
        T ydot;
        ydot = this->y.at(n-1);
        ydot*= -1.;
        ydot += this->y.at(n);
        ydot *= 1./(double)num_stiff_ministeps;    
        // Set transient_y  = [y[n]-3*ydot,y[n]-2*ydot,y[n]-ydot,y[n]] for order 3. Value at the first index will be given the value at the next ministep, and so on.
        for(int i = 0; i < this->order; i++){  // Technically we don't need to fill the first idx, but we may as well.
            // y[n-i] = y[n] - ydot*i   (here y and n refer to the intermediate steps)
            y_transient.at(mini_n-i) = this->y[n];
            T tmp;
            tmp = ydot; tmp*=-i;
            y_transient[mini_n-i] += tmp;
        }
    }
    // else if <can get a whole number of steps>{

    //}
    else{
        mini_dt = this->dt;
        // Not enough ministeps to remain within a single step.
        for(int i = 0; i < this->order; i++){
            y_transient[mini_n-i] = this->y[n-i];
        }
    }
    // sample to catch the error of y(n) != y_transient(mini_n).
    assert(this->y[n].F[0] == y_transient[mini_n].F[0]);
    if (this->y[n].F.container_size() > 3)
        assert(this->y[n].F[3] == y_transient[mini_n].F[3]);
    if (this->y[n].F.container_size() > 5)
        assert(this->y[n].F[5] == y_transient[mini_n].F[5]);
}


// TODO generalise with lagrange polynomial and larger y_transient?
/**
 * @brief Interpolates f(x) to points given by new_x using barycentric form of lagrange polynomial
 * @note Questionable if worth it. 
 */

template<typename T>
std::vector<T> Hybrid<T>::lagrange_interpolate(const std::vector<double> x,  const std::vector<T> y, const std::vector<double> new_x){
    // y(x) = sum(wn/(x-x_n).y_n)/D, where D = sum(w_n/(x-x_n))
    double D = 0; 
    std::vector<T> new_y(new_x.size()); 
    // Compute (inverse) weights
    std::vector<double> inv_w(new_x.size()); 
    assert(new_x.size()>1);
    for (size_t n = 0; n < new_x.size(); n++){
        inv_w[n] = 1;
        for(size_t m = 0; m < new_x.size(); m++){
            if (m == n) continue;
            inv_w[n] *= (x[n]-x[m]);
        }
    }
    for(size_t i = 0; i < new_x.size(); i++){
        bool found_node = false;
        for (size_t n = 0; n < new_x.size(); n++){
            // if seek value at a node, return value at that node.
            if (new_x[i] == x[n]){
                new_y[i] = y[n]; 
                found_node = true;
                break;
            }
        }
        if (found_node) continue;
        // Calculate W = sum w_n/(x-x_n).y_n, and M  = sum w_n/(x-x_n)
        T W = this->zero_y;
        double M = 0;
        for (size_t k = 0; k < inv_w.size();k++){
            double M_k = 1/((new_x[i]-x[k])*inv_w[k]);
            T W_k = y[k];
            W_k*=M_k;
            W += W_k;
            M += M_k;
        }
        new_y[i] = W;
        new_y[i]*= (1/M);
    }
    return new_y;
}

/**
 * @brief 
 * @details Assumes that the time steps between transient y is indeed mini_dt.
 * @tparam T 
 * @param n 
 * @param num_ministeps 
 */
template<typename T>
void Hybrid<T>::modify_ministeps(const int n,const int num_ministeps){     
    #ifdef NO_MINISTEP_UPDATING
    std::runtime_error("Unexpected call to ministep modification when ministep updating was turned off");
    #else
    assert(this->y[n].F[3] == old_y_transient[(old_mini_n)%(this->order)].F[3]);
    
    double old_mini_dt = mini_dt;

    num_stiff_ministeps = num_ministeps;
    if ((this->order-1)/num_stiff_ministeps <= 1){
        mini_dt = this->dt/(double)num_stiff_ministeps;  
        std::vector<double> old_times(old_y_transient.size());
        std::vector<double> new_times(old_y_transient.size());
        std::vector<T> old_y(old_y_transient.size());
        
        for(int i = 0; i < this->order; i++){ 
            // unwind last few ministeps of last step (in old_y_transient)
            int old_idx = (old_mini_n-i)%(this->order);
            int new_idx = (mini_n-i)%(this->order);
            old_y[old_idx] = old_y_transient[old_idx];
            // old_times correspond to same indices as old_y
            old_times[old_idx] = this->t[n] - old_mini_dt*i; 
            new_times[new_idx] = this->t[n] - mini_dt*i;
        }
        y_transient.resize(this->order);
        y_transient = lagrange_interpolate(old_times,old_y,new_times);
    }
    else{
        // Not enough ministeps to remain within a single step.
        mini_dt = this->dt;
        for(int i = 0; i < this->order; i++){
            y_transient[(mini_n-i)%(this->order)] = this->y[n-i];
        }        
    }
    assert(this->y[n].F[3] == y_transient[(mini_n)%(this->order)].F[3]); // test it
    #endif
}
template<typename T>
void Hybrid<T>::initialise_transient_y_v2(int n){ 
    mini_n = this->order-1;
    if ((this->order-1)/num_stiff_ministeps <= 1){
        mini_dt = this->dt/(double)num_stiff_ministeps; 
        std::vector<double> old_times(y_transient.size());
        std::vector<double> new_times(y_transient.size());    

        for(int i = 0; i < this->order; i++){ 
            old_times[mini_n-i] = this->t[n-i];
            new_times[mini_n-i] = this->t[n] - mini_dt*i;
        }

        y_transient.resize(this->order);
        y_transient = lagrange_interpolate();
    }
    // else if <can get a whole number of steps>{

    //}    
    else{
        // Not enough ministeps to remain within a single step.
        mini_dt = this->dt;
        for(int i = 0; i < this->order; i++){
            y_transient[mini_n-i] = this->y[n-i];
        }
    }
    assert(this->y[n].F[3] == y_transient[(mini_n)%(this->order)].F[3]);  
}



template<typename T>
void Hybrid<T>::backward_Euler(unsigned n){
    // Assumes that y_n+1 contains a guess based on sys, and estimates a solution to y_n
    unsigned idx = 0;
    //old = y[n+1] - y[n] 
    T old = this->y[n];
    old *= -1.;
    old += this->y[n+1];  // 
    //old caches the step from the regular part

    // Guess. (Naive Euler)
    sys_ee(this->y[n+1], this->y[n+1]);
    this->y[n+1] *= this->dt;
    this->y[n+1] += this->y[n]; 

    // Picard iteration TODO seek better way to do this - reduce modularisation - reduce to solving linear system?
    double diff = stiff_rtol*2;
    while (diff > stiff_rtol && idx < stiff_max_iter){ // TODO more efficient mapping to use by considering known form of equation  c.f. ELENDIF
        T tmp = this->y[n+1];
        sys_ee(tmp, this->y[n+1]);
        this->y[n+1] *= this->dt;
        this->y[n+1] += this->y[n];
        
        tmp *= -1;
        tmp += this->y[n+1];
        diff = tmp.norm()/this->y[n].norm();
        idx++;
    }
    this->y[n+1] += old;
    if(idx==stiff_max_iter) std::cerr<<"Max Euler iterations exceeded (B)"<<std::endl;
}


    





} // end of namespace 'ode'


#endif
