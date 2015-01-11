#pragma once
#include <vector>
class Grid; class Settings; class ParticleMover; class Geometry; class Communication; class State;
class System {
private:
    bool m_isInitialized;
    Grid *m_grid;
    Settings *m_settings;
    ParticleMover *m_particleMover;
    Geometry *m_geometry;
    Communication *m_communication;
    State *m_state;
public:
    System();
    ~System();
    void initialize();
};
