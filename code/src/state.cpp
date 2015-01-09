#include "state.h"

State::State() :
    numberOfParticlesAllocatedMemory(0),
    numberOfParticles(0),
    x(NULL), y(NULL), z(NULL),
    vx(NULL), vy(NULL), vz(NULL),
    mass(NULL)
{
    resize(1e5); // Initial array size
}

State::~State()
{

}

void State::resize(unsigned int maxNumberOfParticles)
{
    // Allocate new memory
    double *newX    = new double[maxNumberOfParticles];
    double *newY    = new double[maxNumberOfParticles];
    double *newZ    = new double[maxNumberOfParticles];
    double *newVx   = new double[maxNumberOfParticles];
    double *newVy   = new double[maxNumberOfParticles];
    double *newVz   = new double[maxNumberOfParticles];
    double *newMass = new double[maxNumberOfParticles];

    // Copy current data into new location
    memcpy(newX, x, numberOfParticles*sizeof(double));
    memcpy(newY, y, numberOfParticles*sizeof(double));
    memcpy(newZ, z, numberOfParticles*sizeof(double));
    memcpy(newVx, vx, numberOfParticles*sizeof(double));
    memcpy(newVy, vy, numberOfParticles*sizeof(double));
    memcpy(newVz, vz, numberOfParticles*sizeof(double));
    memcpy(newMass, mass, numberOfParticles*sizeof(double));

    // Free memory for previous allocations
    delete[] x; delete[] y; delete[] z;
    delete[] vx; delete[] vy; delete[] vz;
    delete[] mass;

    // Point data to these new pointers
    x = newX; y = newY; z = newZ;
    vx = newVx; vy = newVy; vz = newVz;
    mass = newMass;
    numberOfParticlesAllocatedMemory = maxNumberOfParticles;
}

vec3 State::position(unsigned int particleIndex)
{
    return vec3(x[particleIndex], y[particleIndex], z[particleIndex]);
}

vec3 State::velocity(unsigned int particleIndex)
{
    return vec3(vx[particleIndex], vy[particleIndex], vz[particleIndex]);
}

