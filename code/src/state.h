#pragma once
#include "vec3.h"

class State
{
private:
    unsigned int numberOfParticlesAllocatedMemory;
    unsigned int numberOfParticles;
    double *x;
    double *y;
    double *z;
    double *vx;
    double *vy;
    double *vz;
    double *mass;

public:
    State();
    ~State();
    void resize(unsigned int maxNumberOfParticles);
    vec3 position(unsigned int particleIndex);
    vec3 velocity(unsigned int particleIndex);
};
