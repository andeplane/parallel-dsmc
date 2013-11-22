#include "colliderspecular.h"
#include <iostream>

ColliderSpecular::ColliderSpecular()
{
}

void ColliderSpecular::collide(Random *rnd, data_type *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2, bool print_details)
{
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

    data_type v_normal   = n_x*v[0] + n_y*v[1] + n_z*v[2];
    data_type v_tangent1 = t1_x*v[0] + t1_y*v[1] + t1_z*v[2];
    data_type v_tangent2 = t2_x*v[0] + t2_y*v[1] + t2_z*v[2];

    v[0] = -v_normal*n_x + v_tangent1*t1_x + v_tangent2*t2_x;
    v[1] = -v_normal*n_y + v_tangent1*t1_y + v_tangent2*t2_y;
    v[2] = -v_normal*n_z + v_tangent1*t1_z + v_tangent2*t2_z;
}
