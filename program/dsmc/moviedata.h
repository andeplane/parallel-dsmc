#pragma once
#include <vector>
#include <string>
#include <defines.h>

using std::vector;
using std::string;
class CVector;
class CTexture;
class Solver;
class Timestep
{
public:
    Timestep();
    Timestep *next;
    Timestep *previous;
    Solver *solver;
    vector<double> positions;
    vector<double> velocities;
    void add_molecule_data(vector<float> new_positions);
    void render_points();
    void render_billboards(CTexture *texture, float scale);
};

class MovieData
{
public:
    int cpus;
    int timesteps;
    Timestep *first_timestep;
    MovieData(int cpus_, int timesteps_);
    void load_movie_files(string state_folder);
};
