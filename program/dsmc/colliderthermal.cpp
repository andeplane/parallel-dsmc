#include "colliderthermal.h"
#include <math.h>
#include <random.h>


ColliderThermal::ColliderThermal(double sqrt_wall_temp_over_mass_)
{
    num_collisions = 0;
    sqrt_wall_temp_over_mass = sqrt_wall_temp_over_mass_;
}

void ColliderThermal::collide(Random *rnd, data_type *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2)
{
    data_type v_normal   = sqrt_wall_temp_over_mass*sqrt(-2*log(rnd->next_double()));
    data_type v_tangent1 = sqrt_wall_temp_over_mass*rnd->next_gauss();
    data_type v_tangent2 = sqrt_wall_temp_over_mass*rnd->next_gauss();

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
    num_collisions++;
}
