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
    void set_topology();
    void sync_mpi_initialize();
    void set_initial_positions_and_mark_as_not_moved();
    bool molecule_is_on_this_node(double *r);
    int neighbor_index_of_molecule(double *r);
    void mpi_move();
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
    double *mpi_send_buffer;
    double *r;
    double *v;
    double *r0;
    unsigned long *molecule_index_in_cell;
    unsigned long *molecule_cell_index;

    double *tmp_r;
    double *tmp_v;
    bool   *tmp_molecule_moved;
    bool   *molecule_moved;

    unsigned long num_molecules_local;
    unsigned long num_molecules_global;

    double reservoir_size;
    double grid_origo_x, grid_origo_y, grid_origo_z;
    double origo[3];
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
    double porosity, porosity_global;
    double volume;
    double cell_length[3];
    int node_index[3];
    int num_processors[3];
    int neighbor_nodes[6];
    int num_cells[3];
    int num_cells_per_node_total;
    int num_cells_total;
    int num_cells_per_node[3];
    int num_moved_molecules_indices[6];
    int *moved_molecules_indices[6];
    double node_length[3];
    short my_parity[3];

    unsigned long collisions;
    unsigned long num_active_cells;
	int steps;
    int myid;
    int num_nodes;

    void initialize(Settings *settings_, int myid_);
	void step();
    System() { }
};
