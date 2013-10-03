#pragma once
#include <vector>
#include <string>
class System;
class DSMC_IO;
class CVector;

using namespace std;

typedef enum {
    voxel_type_empty = 0,
    voxel_type_wall = 1,
    voxel_type_boundary = 2
} voxel_type;

class Grid
{
public:
    int Nx;
    int Ny;
    int Nz;
    int points;

    vector<double> voxel_size;
    vector<CVector> unit_normal_vectors;
    vector<CVector> point_list;

    System *system;
    unsigned char *voxels;
    float *normals;
    float *tangents1;
    float *tangents2;

    Grid(string filename, System *system_);
    unsigned char *get_voxel(const int &i, const int &j, const int &k);
    unsigned char *get_voxel(const double &x, const double &y, const double &z);
    unsigned char *get_voxel(double *r);
    void get_index_vector_from_index(const int &index, int &i, int &j, int &k);
    int get_index_of_voxel(double *r);
    double get_time_until_collision(double *r, double *v, const int &voxel_index);
    void read_matrix(string filename, DSMC_IO *io);
};
