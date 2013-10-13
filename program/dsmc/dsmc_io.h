#pragma once
#include <stdio.h>
#include <fstream>

using namespace std;

class Settings;
class System;
class Grid;

class DSMC_IO
{
public:
    Settings *settings;
    System *system;

    DSMC_IO(System *system_);
    void save_state_to_file_binary();
    void load_state_from_file_binary();
    void save_state_to_movie_file();
    void finalize();
    void read_grid_matrix(string filename, Grid *grid);
    bool movie_file_open;
    int  movie_frames;
    double *data;

    FILE *energy_file;
    FILE *velocity_file;
    FILE *flux_file;
    FILE *permeability_file;
    FILE *density_file;
    FILE *linear_density_file;
    FILE *num_molecules_file;
    FILE *pressure_file;
    FILE *linear_pressure_file;
    ofstream *time_statistics_file;
    ofstream *movie_file;
};
