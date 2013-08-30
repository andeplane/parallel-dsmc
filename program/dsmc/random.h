#pragma once

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

class Random {
private:
    inline double bessel_i0(double x);
public:
    long     iy;
    long     iv[NTAB];
    long     idum[1];
    double alpha_n;
    double alpha_t;
    double sqrt_one_minus_alpha;
    double sqrt_one_minus_alpha_over_alpha;
    double sqrt_alpha_over_two;
    long   cercignani_lampis_normal_component_trials;

    Random(long seed, double alpha_n_, double alpha_t_);
    double next_double();
    double next_gauss();
    double next_cercignani_lampis_normal_component(double v_norm_in, double factor, double max_v_out);
};
