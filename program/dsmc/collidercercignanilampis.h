#pragma once
#include <colliderbase.h>

class Random;
class ColliderCercignaniLampis : public ColliderBase
{
public:
    ColliderCercignaniLampis(double sqrt_wall_temp_over_mass_);
    virtual void collide(Random *rnd, double *v_in, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2);
};
