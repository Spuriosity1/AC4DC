#include "../src/Plasma.h"
#include <cstdio>

void print_est(elec_state_t &es) {
    std::printf("{N=%lf, E=%lf, Np=%lf, Ep=%lf}\n",es.N, es.E, es.Np, es.Ep);
}

int main(int argc, char const *argv[]) {
    elec_state_t a, b, c;
    a = 0;
    printf("a=0 ==> ");
    print_est(a);

    printf("\nb = {12, 30, 4, 3} ==> ");
    b.N = 12;
    b.E = 30;
    b.Np = 4;
    b.Ep = 3;
    print_est(b);

    printf("\nc = {-2, 0, -4, -10} ==> ");
    c.N = -2;
    c.E = 0;
    c.Np = -4;
    c.Ep = -12;
    print_est(c);

    printf("\nb = b*2  ==> ");
    b = b*2;
    print_est(b);

    printf("\n c *= 0.4 ==> ");
    c *= 0.4;
    print_est(c);

    printf(" a= b * c ==> ");
    a = b*c;
    print_est(a);

    return 0;
}
