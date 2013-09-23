#pragma once

class Cell;
class Sorter;
class Grid;
class DSMC_IO;
class Random;
class Settings;
class UnitConverter;
class DSMCTimer;
class MoleculeMover;

#include <iostream>
#include <fstream>
#include <vector>
#include <cinifile.h>

#define MAX_MOLECULE_NUM 10000000

using namespace std;

class System {
private:
    void init_positions();
    void init_velocities();
    void init_molecules();
    void init_cells();
	void move();
    void init_randoms();
    void collide();
	void accelerate();
    void maintain_pressure_A();
    void maintain_pressure_B();
    void maintain_pressure();
    bool remove_molecule_in_pressure_reservoir(bool remove_from_source);
    void find_position_in_reservoirs(double *r, bool find_position_in_source);
    void add_molecule_in_pressure_reservoirs(bool add_in_source);
    inline void find_position(double *r);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    void update_cell_volume();
    void setup_molecules();
    void setup_cells();
    void calculate_porosity();
    void update_molecule_cells();
    void count_reservoir_particles();
public:
    int cell_index_from_position(double *r);

    DSMC_IO *io;
    DSMCTimer *timer;
    Grid *world_grid;
    Settings *settings;
    UnitConverter * unit_converter;
    Random *rnd;
    MoleculeMover *mover;

    vector<Cell*> active_cells;
    vector<Cell*> all_cells;
    vector<Cell*> reservoir_A_cells;
    vector<Cell*> reservoir_B_cells;

    double *mpi_receive_buffer;
    double *r;
    double *v;
    double *r0;
    unsigned long *molecule_index_in_cell;
    unsigned long *molecule_cell_index;

    int num_molecules;

    double reservoir_size;
    double grid_origo_x, grid_origo_y, grid_origo_z;
    double length[3];
    double half_length[3];
    double atoms_per_molecule;
	double mpv; 	// Most probable velocity
	double mfp; 	// Mean free path
	double dt;
	double t;
    double temperature;
    double diam, density;
    double wall_temperature;
    double cell_length_x;
    double cell_length_y;
    double cell_length_z;
    double porosity;
    double volume;

    unsigned long collisions;
	int steps;
    int myid;
    int cells_x, cells_y, cells_z;
    int num_cells_vector[3];
    long reservoir_b_particle_count;
    long flux_count;

    void initialize(Settings *settings_, int myid_);
	void step();
    System() { }
};
