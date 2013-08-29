#pragma once
#include <colliderbase.h>
class Random;

class ColliderSpecular : public ColliderBase
{
public:
    ColliderSpecular();
    virtual void collide(Random *rnd, double *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2);
};
