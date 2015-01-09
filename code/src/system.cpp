#include "system.h"

System::System() {

}

System::~System() {
    delete m_grid;
    delete m_settings;
    delete m_random;
    delete m_particleMover;
    delete m_geometry;
    delete m_communication;
    delete m_state;
}
