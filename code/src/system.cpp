#include "system.h"
#include "grid.h"
#include "settings.h"
#include "particlemover.h"
#include "geometry.h"
#include "communication.h"
#include "state.h"

System::System() :
    m_isInitialized(false),
    m_grid(NULL),
    m_settings(NULL),
    m_particleMover(NULL),
    m_geometry(NULL),
    m_communication(NULL),
    m_state(NULL)
{

}

System::~System() {
    delete m_grid;
    delete m_settings;
    delete m_particleMover;
    delete m_geometry;
    delete m_communication;
    delete m_state;
}

void System::initialize()
{
    if(!m_settings) m_settings = Settings::createDefault();
    if(!m_grid) m_grid = Grid::createDefault();
    if(!m_particleMover) m_particleMover = ParticleMover::createDefault();
    if(!m_geometry) m_geometry = Geometry::createDefault();
    if(!m_communication) m_communication = Communication::createDefault();
    if(!m_state) m_state = State::createDefault();

    m_isInitialized = true;
}
