#include <grid.h>

#include <fstream>
#include <system.h>
#include <dsmc_io.h>
#include <cvector.h>
#include <topology.h>

#define min(a,b)                      (((a) < (b)) ? (a) : (b))
#define max(a,b)                      (((a) > (b)) ? (a) : (b))
#define clamp(value, lb, ub)          max( lb, min( ub, value ))

void Grid::read_matrix(string filename, DSMC_IO *io) {
    io->read_grid_matrix(filename, this);
}

Grid::Grid(string foldername, System *system_)
{
    char filename[1000];
    system = system_;
    sprintf(filename, "%s/%04d.bin",foldername.c_str(), system->myid);

    read_matrix(filename,system->io);
    voxel_size.resize(3);
    voxel_size[0] = system->length[0] / global_nx;
    voxel_size[1] = system->length[1] / global_ny;
    voxel_size[2] = system->length[2] / global_nz;

    // The unit normal vectors are 6 vectors pointing in the following directions
    // x-, x+, y-, y+, z-, z+
    unit_normal_vectors.resize(6,CVector(0,0,0));

    unit_normal_vectors[0].x = -1;
    unit_normal_vectors[1].x = 1;
    unit_normal_vectors[2].y = -1;
    unit_normal_vectors[3].y = 1;
    unit_normal_vectors[4].z = -1;
    unit_normal_vectors[5].z = 1;

    point_list.resize(8, CVector(0,0,0));
}

unsigned char *Grid::get_voxel(const CVector &voxel_indices) {
    return get_voxel(int(voxel_indices.x), int(voxel_indices.y), int(voxel_indices.z));
}

unsigned char *Grid::get_voxel(const int &i, const int &j, const int &k) {
    return &voxels[i*ny*nz + j*nz + k];
}

unsigned char *Grid::get_voxel(const double &x, const double &y, const double &z) {
    int i =  (x-system->topology->origin[0])*system->one_over_length[0]*global_nx + voxel_origin.x;
    int j =  (y-system->topology->origin[1])*system->one_over_length[1]*global_ny + voxel_origin.y;
    int k =  (z-system->topology->origin[2])*system->one_over_length[2]*global_nz + voxel_origin.z;

    return get_voxel(i,j,k);
}

unsigned char *Grid::get_voxel(vector<data_type> &r, int index) {
    int i =  (r.at(3*index + 0)-system->topology->origin[0])*system->one_over_length[0]*global_nx + voxel_origin.x;
    int j =  (r.at(3*index + 1)-system->topology->origin[1])*system->one_over_length[1]*global_ny + voxel_origin.y;
    int k =  (r.at(3*index + 2)-system->topology->origin[2])*system->one_over_length[2]*global_nz + voxel_origin.z;

    return get_voxel(i,j,k);
}

void Grid::get_index_vector_from_index(const int &index, int &i, int &j, int &k) {
    i = index/(ny*nz);
    j = (index / nz) % ny;
    k = index % nz;
}


