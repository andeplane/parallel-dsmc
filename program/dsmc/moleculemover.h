#pragma once
class Cell;
class System;
class Random;
class Grid;

class MoleculeMover
{
public:
    System *system;
    unsigned char *voxels;
    Grid *grid;
    long count_periodic[3];

    MoleculeMover();
    void initialize(System *system_);
    void move_molecules(double dt, Random *rnd);
    void move_molecule(const int &idx);
    inline void do_move(double *r, double *v, double *r0, const double &dt);
    void move_molecule(unsigned long &idx, double dt, Random *rnd, int depth);
};
