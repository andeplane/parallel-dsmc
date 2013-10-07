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
class Topology;

#include <iostream>
#include <fstream>
#include <vector>
#include <cinifile.h>
#include <mpi.h>

#define CYLINDER_RADIUS_SQUARED 0.0064
#define BOX_FRACTION 0.2
#define MAX_MOLECULE_NUM 20000000

using namespace std;

class System {
private:
    void init_positions();
    void init_velocities();
    void init_molecules();
    void init_cells();
	void move();
    void mpi_move();
    void init_randoms();
    void collide();
	void accelerate();
    void maintain_pressure_A();
    void maintain_pressure_B();
    void maintain_pressure();
    void find_position_in_reservoirs(double *r, bool find_position_in_source);
    void add_molecules_in_inlet_reservoir(Cell *cell, const double &velocity_std_dev, const int &delta_num_molecules);
    void remove_molecules_in_inlet_reservoir(Cell *cell, const int &delta_num_molecules);
    void add_molecules_in_outlet_reservoir(Cell *cell, const double &velocity_std_dev, const int &delta_num_molecules);
    void remove_molecules_in_outlet_reservoir(Cell *cell, const int &delta_num_molecules);
    void add_molecule_to_cell(Cell *cell, const int &molecule_index);
    void add_molecules_from_mpi(vector<double> &data, const int &num_new_molecules);
    void remove_molecule_from_system(const long &molecule_index);
    inline void find_position(double *r);
    inline void find_position_in_cell(Cell *cell, double *r);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    void update_cell_volume();
    void setup_molecules();
    void setup_cells();
    void calculate_porosity();
    void calculate_global_porosity();
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
    Topology *topology;

    vector<Cell*> active_cells;
    vector<Cell*> all_cells;
    vector<Cell*> reservoir_A_cells;
    vector<Cell*> reservoir_B_cells;

    vector<double> mpi_receive_buffer;
    vector<double> mpi_send_buffer;
    vector<double> r;
    vector<double> v;
    vector<double> r0;
    vector<unsigned long> molecule_index_in_cell;
    vector<unsigned long> molecule_cell_index;
    vector<int> node_num_new_molecules;
    vector<vector<double> > node_molecule_data;

    long num_molecules_local;
    long num_molecules_global;

    double grid_origo_x, grid_origo_y, grid_origo_z;
    double length[3];
    double one_over_length[3];
    double half_length[3];
    double atoms_per_molecule;
    double most_probable_velocity; 	// Most probable velocity
	double dt;
	double t;
    double temperature;
    double diam, density;
    double wall_temperature;
    double cell_length_x;
    double cell_length_y;
    double cell_length_z;
    double porosity;
    double porosity_global;
    double volume;
    double volume_global;

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
