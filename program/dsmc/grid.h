#pragma once
#include <vector>
#include <string>
#include <cvector.h>
class System;
class DSMC_IO;

using namespace std;

typedef enum {
    voxel_type_empty = 0,
    voxel_type_wall = 1,
    voxel_type_boundary = 2
} voxel_type;

class Grid
{
public:
    unsigned int nx; // num voxels in this matrix, includes a copy for each of the 26 neighboring nodes
    unsigned int ny;
    unsigned int nz;
    unsigned int global_nx; // total number of voxels in the global matrix
    unsigned int global_ny;
    unsigned int global_nz;
    unsigned int nx_per_cpu; // num unique voxels per cpu
    unsigned int ny_per_cpu;
    unsigned int nz_per_cpu;
    float nx_divided_by_global_nx;
    float ny_divided_by_global_ny;
    float nz_divided_by_global_nz;
    unsigned int num_voxels;

    vector<double> voxel_size;
    vector<CVector> unit_normal_vectors;
    vector<CVector> point_list;

    System *system;
    unsigned char *voxels;
    float *normals;
    float *tangents1;
    float *tangents2;
    CVector voxel_origin;

    Grid(string filename, System *system_);
    unsigned char *get_voxel(const CVector &voxel_indices);
    unsigned char *get_voxel(const int &i, const int &j, const int &k);
    unsigned char *get_voxel(const double &x, const double &y, const double &z);
    unsigned char *get_voxel(double *r);
    void get_index_vector_from_index(const int &index, int &i, int &j, int &k);
    int get_index_of_voxel(double *r);
    double get_time_until_collision(double *r, double *v, const int &voxel_index);
    void read_matrix(string filename, DSMC_IO *io);
};
