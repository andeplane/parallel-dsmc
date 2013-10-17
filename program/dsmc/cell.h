#pragma once
#include <vector>

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
	double vr_max;
    double collision_coefficient;
    double origin[3];
    int index_vector[3];
    vector<double> average_velocity;

    System *system;
    int *molecules;
    int num_molecules_allocated_memory;
    int num_molecules;

    Cell(System *system);
	void reset();
    unsigned long prepare();
    void resize(int n);
    int collide(Random *rnd);
    inline void collide_molecules(double *v0, double *v1, const double &v_rel, Random *rnd);
    void update_volume();
    void add_molecule(const int &molecule_index, vector<unsigned long> &index_in_cell, vector<unsigned long> &cell_index);
    void remove_molecule(const int &molecule_index, vector<unsigned long> &index_in_cell);
    vector<double> &update_average_velocity();

    double calculate_kinetic_energy();
    static bool cmp(Cell *c1, Cell *c2);
};


