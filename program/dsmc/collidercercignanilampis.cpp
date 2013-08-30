#include "collidercercignanilampis.h"
#include <random.h>
#include <math.h>
#include <system.h>
#include <unitconverter.h>

ColliderCercignaniLampis::ColliderCercignaniLampis(double sqrt_wall_temp_over_mass_, System *system_)
{
    sqrt_wall_temp_over_mass = sqrt_wall_temp_over_mass_;
    system = system_;
}

void ColliderCercignaniLampis::collide(Random *rnd, double *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2)
{
    double k = 1.3806488e-23;
    double T = 100;
    double m = 6.63352088e-26;
    double factor = sqrt(2*k*T/m);
    factor = sqrt(2)*sqrt_wall_temp_over_mass;

    double v_norm_in = v[0]*normal_vector[0] + v[1]*normal_vector[1] + v[2]*normal_vector[2];
    double v_t1_in   = v[0]*tangent_vector_1[0] + v[1]*tangent_vector_1[1] + v[2]*tangent_vector_1[2];
    double v_t2_in   = v[0]*tangent_vector_2[0] + v[1]*tangent_vector_2[1] + v[2]*tangent_vector_2[2];
    // double v_normal   = system->unit_converter->velocity_from_SI(rnd->next_cercignani_lampis_normal_component(system->unit_converter->velocity_to_SI(v_norm_in), factor, 1500));
    double v_max = system->unit_converter->velocity_from_SI(1500);
    double v_normal   = rnd->next_cercignani_lampis_normal_component(v_norm_in, factor, v_max);
    double sigma = 1 - sqrt(1-rnd->alpha_t);
    double v_tangent1 = rnd->next_gauss()*sqrt(rnd->alpha_t)*sqrt_wall_temp_over_mass + v_t1_in*(1 - sigma);
    double v_tangent2 = rnd->next_gauss()*sqrt(rnd->alpha_t)*sqrt_wall_temp_over_mass + v_t2_in*(1 - sigma);

    // Normal vector
    float n_x = normal_vector[0];
    float n_y = normal_vector[1];
    float n_z = normal_vector[2];

    // Tangent vector 1
    float t1_x = tangent_vector_1[0];
    float t1_y = tangent_vector_1[1];
    float t1_z = tangent_vector_1[2];

    // Tangent vector 2
    float t2_x = tangent_vector_2[0];
    float t2_y = tangent_vector_2[1];
    float t2_z = tangent_vector_2[2];

    v[0] = v_normal*n_x + v_tangent1*t1_x + v_tangent2*t2_x;
    v[1] = v_normal*n_y + v_tangent1*t1_y + v_tangent2*t2_y;
    v[2] = v_normal*n_z + v_tangent1*t1_z + v_tangent2*t2_z;
}
