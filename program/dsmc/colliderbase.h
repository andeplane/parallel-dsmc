#pragma once
#include <system.h>

class Random;

class ColliderBase
{
protected:
    double sqrt_wall_temp_over_mass;
public:
    unsigned long num_collisions;
    ColliderBase();
    virtual void collide(Random *rnd, data_type *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2);
};
