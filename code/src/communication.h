#pragma once
#include "vec3.h"

class Communication
{
private:
    vec3 m_origin;
    unsigned int m_numberOfCPUs;
    unsigned int m_myid;
public:
    static Communication *createDefault();
    Communication();
    ~Communication();
};
