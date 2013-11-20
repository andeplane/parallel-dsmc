#pragma once
#include <colliderbase.h>

class Random;
class System;

class ColliderCercignaniLampis : public ColliderBase
{
private:
    System *system;
public:
    ColliderCercignaniLampis(double sqrt_wall_temp_over_mass_, System *system_);
    virtual void collide(Random *rnd, data_type *v_in, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2);
};
