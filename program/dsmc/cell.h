#pragma once
#include <vector>
#include <defines.h>

class Random;
class System;
class Cell;

using namespace std;

class Cell {
public:
    double volume;
    int pixels; // Used to calculate volume in a cell
    int total_pixels;
    int index;
    int collision_pairs;
    double collision_rest;
    bool is_reservoir;

    double x0, y0, z0;
    double Lx, Ly, Lz;
    data_type vr_max;
    double collision_coefficient;
    double origin[3];
    int index_vector[3];
    vector<data_type> average_velocity;

    System *system;
    int *molecules;
    int num_molecules_allocated_memory;
    int num_molecules;

    Cell(System *system);
    ~Cell() {
        average_velocity.clear();
        delete molecules;
    }
	void reset();
    unsigned long prepare();
    void resize(int n);
    int collide(Random *rnd);
    inline void collide_molecules(data_type *v0, data_type *v1, const data_type &v_rel, Random *rnd);
    void update_volume();
    void add_molecule(const int &molecule_index, unsigned long *index_in_cell, unsigned long *cell_index);
    void remove_molecule(const int &molecule_index, unsigned long *index_in_cell);
    vector<data_type> &update_average_velocity();

    double calculate_kinetic_energy();
    static bool cmp(Cell *c1, Cell *c2);
};


