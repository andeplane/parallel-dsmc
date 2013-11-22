#pragma once
#include <colliderbase.h>

class ColliderSpecular : public ColliderBase
{
public:
    ColliderSpecular();
    virtual void collide(Random *rnd, data_type *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2, bool print_details=false);
};
