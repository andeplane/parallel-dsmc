#pragma once
#include <vector>

class System {
private:
    Grid *m_grid;
    Settings *m_settings;
    Random *m_random;
    ParticleMover *m_particleMover;
    Geometry *m_geometry;
    Communication *m_communication;
    State *m_state;
public:
    System();
    ~System();
};
