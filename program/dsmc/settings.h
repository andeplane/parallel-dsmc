#pragma once

#include <cinifile.h>

typedef enum {
    velocity_profile_box = 0,
    velocity_profile_cylinder = 1,
    velocity_profile_other = 2
} velocity_profile;

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
    int cells_x;
    int cells_y;
    int cells_z;
    int flow_direction;
    int velocity_bins;
    int nx, ny, nz;
    int seed;

    double viscosity;
    double mass;
    double temperature;
    double wall_temperature;
    double gravity;
    double density;
    double diam;
    double Lx;
    double Ly;
    double Lz;
    double pressure_A;
    double pressure_B;
    double dt;
    string velocity_profile_type;
    string surface_interaction_model;
    double alpha_n, alpha_t;

};
