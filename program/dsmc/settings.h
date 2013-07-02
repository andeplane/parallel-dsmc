#pragma once

#include <cinifile.h>

class Settings
{
public:
    Settings(string filename);
    CIniFile ini_file;
    bool load_previous_state;
    bool create_movie;
    bool maintain_pressure;
    int statistics_interval;
    int movie_every_n_frame;
    int atoms_per_molecule;
    int movie_molecules;
    int timesteps;
    int cells_per_node_x, cells_per_node_y, cells_per_node_z;
    int gravity_direction;
    int nodes_x, nodes_y, nodes_z;
    double viscosity;
    double mass;
    double temperature;
    double wall_temperature;
    double gravity;
    double density;
    double diam;
    double Lx, Ly, Lz;
    double reservoir_fraction;
    double pressure_A, pressure_B;
    double dt;
};
