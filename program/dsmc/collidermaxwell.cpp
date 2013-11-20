#include "collidermaxwell.h"
#include <math.h>
#include <random.h>
#include <iostream>
using namespace std;

ColliderMaxwell::ColliderMaxwell(double sqrt_wall_temp_over_mass_)
{
    sqrt_wall_temp_over_mass = sqrt_wall_temp_over_mass_;
}

void ColliderMaxwell::collide(Random *rnd, data_type *v, float *normal_vector, float *tangent_vector_1, float *tangent_vector_2)
{
    cout << "Maxwell collision model is not implemented yet." << endl;
}
