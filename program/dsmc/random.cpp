#include <random.h>
#include <math.h>
#include <iostream>
#include <math.h>

using namespace std;

Random::Random(long seed, double alpha_n_, double alpha_t_) {
   *idum = seed;
   alpha_n = alpha_n_;
   alpha_t = alpha_t_;
   sqrt_one_minus_alpha = sqrt(1-alpha_n);
   sqrt_one_minus_alpha_over_alpha = sqrt(1-alpha_n)/alpha_n;
   sqrt_alpha_over_two = sqrt(alpha_n/2);
   iy = 0;
   cercignani_lampis_normal_component_trials = 0;
}

double Random::next_gauss() {

    return sqrt( -2.0*log(1.0 - next_double()) )
              * cos( 6.283185307 * next_double() );
}

inline double Random::bessel_i0(double x) {
   double ax, ans;
   double y;

   if ((ax=fabs(x)) < 3.75) {
      y = x/3.75;
      y *=y;
      ans = 1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492+y*(0.2659732+y*(0.360768e-1+y*0.45813e-2)))));
   } else {
      y = 3.75/ax;
      ans = (exp(ax)/sqrt(ax))*(0.39894228+y*(0.1328592e-1+y*(0.225319e-2+y*(-0.157565e-2+y*(0.916281e-2+y*(-0.2057706e-1+y*(0.2635537e-1+y*(-0.1647633e-1+y*0.392377e-2))))))));
   }

   return ans;
}

double Random::next_double()
{
   int             j;
   long            k;
   double          temp;

   if (*idum <= 0 || !iy) {
      if (-(*idum) < 1) *idum=1;
      else *idum = -(*idum);
      for(j = NTAB + 7; j >= 0; j--) {
         k     = (*idum)/IQ;
         *idum = IA*(*idum - k*IQ) - IR*k;
         if(*idum < 0) *idum += IM;
         if(j < NTAB) iv[j] = *idum;
      }
      iy = iv[0];
   }
   k     = (*idum)/IQ;
   *idum = IA*(*idum - k*IQ) - IR*k;
   if(*idum < 0) *idum += IM;
   j     = iy/NDIV;
   iy    = iv[j];
   iv[j] = *idum;
   if((temp=AM*iy) > RNMX) return RNMX;
   else return temp;
}

double Random::next_cercignani_lampis_normal_component(double v_norm_in, double factor, double max_v_out) {
    max_v_out /= factor;
    v_norm_in /= factor;
    double std_dev = sqrt_alpha_over_two;
    double mu = sqrt_one_minus_alpha*v_norm_in;
    while(true) {
      cercignani_lampis_normal_component_trials++;

      double v_out = std_dev*next_gauss() + mu;
      double random_uniform = next_double();
      double bessel_argument = 2*sqrt_one_minus_alpha_over_alpha*v_norm_in*v_out;
      double b_max = max_v_out;

      double b = v_out*exp(-bessel_argument)*bessel_i0(bessel_argument)/b_max;
      if(b>=random_uniform) return v_out*factor;
    }
}
