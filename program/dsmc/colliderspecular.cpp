#include "colliderspecular.h"
#include <iostream>

ColliderSpecular::ColliderSpecular()
{
}

void ColliderSpecular::collide(Random *rnd, double *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2)
{
    // Vector along normal vector
    float n_x = normal_vector[0]*v[0];
    float n_y = normal_vector[1]*v[1];
    float n_z = normal_vector[2]*v[2];

    // Vector along tangent vector 1
    float t1_x = tangent_vector_1[0]*v[0];
    float t1_y = tangent_vector_1[1]*v[1];
    float t1_z = tangent_vector_1[2]*v[2];

    // Vector along tangent vector 2
    float t2_x = tangent_vector_2[0]*v[0];
    float t2_y = tangent_vector_2[1]*v[1];
    float t2_z = tangent_vector_2[2]*v[2];

    v[0] = -n_x + t1_x + t2_x;
    v[1] = -n_y + t1_y + t2_y;
    v[2] = -n_z + t1_z + t2_z;
}
