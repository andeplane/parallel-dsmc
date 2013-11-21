#pragma once
#include <mpi.h>
#include <defines.h>

class Cell;
class Sorter;
class Grid;
class DSMC_IO;
class Random;
class Settings;
class DSMCTimer;
class MoleculeMover;
class Topology;

#include <iostream>
#include <fstream>
#include <vector>
#include <cinifile.h>
#include <mpi.h>
#include <unitconverter.h>
#include <stdexcept>      // std::out_of_range

#define CYLINDER_RADIUS_SQUARED 0.0064
#define BOX_FRACTION 0.2
#define MAX_MPI_DATA 100000

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
    void add_molecule_to_cell(Cell *cell, const int &molecule_index);
    void add_molecules_from_mpi(vector<data_type> &data, const int &num_new_molecules);
    void remove_molecule_from_system(const long &molecule_index);
    bool validate_number_of_cells();
    inline void find_position(const int &index);
    inline int cell_index_from_ijk(const int &i, const int &j, const int &k);
    void setup_molecules();
    void setup_cells();
    void update_molecule_cells();
    void count_reservoir_particles();

    inline Cell *cell_that_should_contain_molecule(const int &molecule_index) {
        int global_cell_index = cell_index_from_position(molecule_index);
        int local_cell_index = cell_index_map.at(global_cell_index);
        return active_cells.at(local_cell_index);
    }

    inline Cell *cell_currently_containing_molecule(const int &molecule_index) {
        long global_cell_index;
        try {
            global_cell_index = molecule_cell_index.at(molecule_index);
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
            cout << "Molecule index " << molecule_index << " was larger than molecule_cell_index accepts on node "  << myid << endl;
            exit(1);
        }
        long local_cell_index;
        try {
            local_cell_index = cell_index_map.at(global_cell_index);
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
            cout << "Global cell index " << global_cell_index << " was larger than cell_index_map on node " << myid << endl;
            exit(1);
        }
        try {
            return active_cells.at(local_cell_index);
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
            cout << "Local cell index " << local_cell_index << " was larger than active_cells on node " << myid << endl;
            exit(1);
        }
    }

public:
    int cell_index_from_position(const int &index);

    DSMC_IO *io;
    DSMCTimer *timer;
    Grid *world_grid;
    Settings *settings;
    UnitConverter * unit_converter;
    Random *rnd;
    MoleculeMover *mover;
    Topology *topology;

    vector<Cell*> active_cells;

    vector<data_type> mpi_receive_buffer;
    vector<data_type> r;
    vector<data_type> v;
    vector<int> molecule_index_in_cell;
    vector<int> molecule_cell_index;
    vector<int> node_num_new_molecules;
    vector<int> cell_index_map;
    vector<vector<data_type> > node_molecule_data;
    // data_type **node_molecule_data;
    vector<int> steps_since_collision;

    int MAX_MOLECULE_NUM;
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
    double t0;
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

    double t_in_nano_seconds() { return unit_converter->time_to_SI(t)*1e9; }
    void initialize(Settings *settings_, int myid_);
	void step();
    System() { }
};
