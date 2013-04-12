#include <grid.h>

#include <fstream>
#include <system.h>
#include <dsmc_io.h>

#define min(a,b)                      (((a) < (b)) ? (a) : (b))
#define max(a,b)                      (((a) > (b)) ? (a) : (b))
#define clamp(value, lb, ub)          max( lb, min( ub, value ))

void Grid::read_matrix(string filename, DSMC_IO *io) {
    io->read_grid_matrix(filename, this);
}

Grid::Grid(string filename, System *system_)
{
    system = system_;
    read_matrix(filename,system->io);
}

unsigned char *Grid::get_voxel(const int &i, const int &j, const int &k) {
    if(i < 0 || i >= Nx || j < 0 || j >= Ny || k < 0 || k >= Nz) {
        return &voxels[((i+Nx)%Nx) + ((j+Ny)%Ny)*Nx+ ((k+Nz)%Nz)*Nx*Ny];
    }

    return &voxels[i + j*Nx + k*Nx*Ny];
}

unsigned char *Grid::get_voxel(const double &x, const double &y, const double &z) {
    int i =  x/system->length[0]*Nx;
    int j =  y/system->length[1]*Ny;
    int k =  z/system->length[2]*Nz;

    return get_voxel(i,j,k);
}

unsigned char *Grid::get_voxel(double *r) {
    int i =  r[0]/system->length[0]*Nx;
    int j =  r[1]/system->length[1]*Ny;
    int k =  r[2]/system->length[2]*Nz;

    return get_voxel(i,j,k);
}

int Grid::get_index_of_voxel(double *r) {
    int i =  r[0]/system->length[0]*Nx;
    int j =  r[1]/system->length[1]*Ny;
    int k =  r[2]/system->length[2]*Nz;

    return i + j*Nx + k*Nx*Ny;
}
