#include "src/HybridIntegrator.hpp"
#include <iostream>
#include <cmath>

using namespace std;
using namespace ode;

template<int N>
struct lvector {
    double* X;
    lvector() {
        X = new double[N];
    }
    lvector(const lvector<N>& u) {
        X = new double[N];
        this->operator=(u);
    }
    ~lvector() {
        delete[] X;
    }
    // Critical vector-space devices
    lvector<N>& operator+=(const lvector<N> &s) {
        for (size_t i=0; i<N; i++) {
            X[i] += s.X[i];
        }
        return *this;
    }
    lvector<N>& operator*=(const double x) {
        for (size_t i=0; i<N; i++) {
            X[i] *= x;
        }
        return *this;
    }
    lvector<N>& operator=(const lvector<N> &s) {
        // self-assignment guard
        if (this == &s)
            return *this;

        for (size_t i=0; i<N; i++) {
            X[i] = s.X[i];
        }
        return *this;
    }
    double norm() {
        double n=0;
        for (size_t i=0; i<N; i++) {
            n += abs(X[i]);
        }
        return n;
    }
    // lvector<N>& operator=(lvector<N>&& other) {
    //     // If we're not trying to move the object into itself...
    //     if (this != &other) {
    //       delete[] this->X;  // Delete the string's original data.
    //       this->X = other.X;  // Copy the other string's data into this string.
    //       other.X = nullptr;  // Finally, reset the other string's data pointer.
    //     }
    //     return *this;
    // }

};

typedef lvector<2> vec_t;

class Harmonic : ode::Hybrid<vec_t>{
    // Classic harmonic oscillator.
public:
    Harmonic(double _omega, size_t num_points, double _dt, int _order) : Hybrid<vec_t>(_order) {
        this->omega = _omega;
        this->omega2 = _omega*_omega;
        vec_t init;
        init.X[0] = 1.;
        init.X[1] = 0.;
        this->setup(init, _dt);
        this->solve_dynamics(0, num_points*_dt); // returns final time
        // this->y and this->t are now populated
        comp_analytic();
    }
    void comp_analytic() {
        analytic.resize(this->t.size());
        for (size_t i=0; i<t.size(); i++) {
            analytic[i].X[0] = cos(omega*t[i]);
            analytic[i].X[1] = -omega*sin(omega*t[i]);
        }
    }
    void print() {
        cout<<"t\ty_c\ty_a\tdiff\n"<<endl;
        for (size_t i = 0; i < y.size(); i++) {
            double delta = y[i].X[0] - analytic[i].X[0];
            cout<<t[i]<<'\t'<<y[i].X[0]<<'\t'<<analytic[i].X[0]<<'\t'<<delta<<endl;
        }
    }
    double avgdev() {
        double tmp=0;
        for (size_t i = 0; i < y.size(); i++) {
            tmp += (y[i].X[0] - analytic[i].X[0])*(y[i].X[0] - analytic[i].X[0]);
            //tmp += (y[i].X[1] - analytic[i].X[1])*(y[i].X[1] - analytic[i].X[1]);
        }
        tmp /= (y.size() -1);
        return pow(tmp, 0.5);
    }

protected:
    vector<vec_t> analytic;
    double omega=1;
    double omega2=1;
    void sys_bound(const vec_t& q, vec_t& qdot, const double t) {
        // q[0] = x
        // q[1] = v
        qdot.X[1] = -q.X[0]*omega2;
        // qdot.X[0] = q.X[1];
        // qdot.X[1] = 0;
        qdot.X[0] = 0;
    }
    void sys_ee(const vec_t& q, vec_t& qdot) {
        // q[0] = x
        // q[1] = v
        // qdot.X[1] = -q.X[0]*omega2;
        qdot.X[0] = q.X[1];
        qdot.X[1] = 0;
    }
};

int main(int argc, char const *argv[]) {
    if (argc <3) {
        cerr<<"Usage: abm_verif [omega] [step] [order]"<<endl;
        return 1;
    }
    double w = atof(argv[1]);
    double step = atof(argv[2]);
    int ord = atoi(argv[3]);
    cerr<<"[abm test] "<<ord<<"th order method, h="<<step<<", omega="<<w<<endl;
    Harmonic h(w, 100000, step, ord);
    h.print();
    cerr<<"Average deviation: "<<h.avgdev()<<endl;
    return 0;
}
