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
    voxel_size.resize(3);
    voxel_size[0] = system->length[0] / Nx;
    voxel_size[1] = system->length[1] / Ny;
    voxel_size[2] = system->length[2] / Nz;
}

unsigned char *Grid::get_voxel(const int &i, const int &j, const int &k) {
    if(i < 0 || i >= Nx || j < 0 || j >= Ny || k < 0 || k >= Nz) {
        return &voxels[((i+Nx)%Nx) + ((j+Ny)%Ny)*Nx+ ((k+Nz)%Nz)*Nx*Ny];
    }

    return &voxels[i + j*Nx + k*Nx*Ny];
}

unsigned char *Grid::get_voxel(const double &x, const double &y, const double &z) {
    int i =  x*system->one_over_length[0]*Nx;
    int j =  y*system->one_over_length[1]*Ny;
    int k =  z*system->one_over_length[2]*Nz;

    return get_voxel(i,j,k);
}

unsigned char *Grid::get_voxel(double *r) {
    int i =  r[0]*system->one_over_length[0]*Nx;
    int j =  r[1]*system->one_over_length[1]*Ny;
    int k =  r[2]*system->one_over_length[2]*Nz;

    return get_voxel(i,j,k);
}

int Grid::get_index_of_voxel(double *r) {
    int i =  r[0]*system->one_over_length[0]*Nx;
    int j =  r[1]*system->one_over_length[1]*Ny;
    int k =  r[2]*system->one_over_length[2]*Nz;

    return i + j*Nx + k*Nx*Ny;
}

inline double time_until_intersection_with_plane(vector<double> &point, vector<double> &point_in_plane, vector<double> &normal) {

}

double Grid::get_time_until_collision(double *r, double *v) {
    /*
     * Strategy:
     *          First calculate time until intersection with every facet of the voxel, choose the face with lowest TOC.
     *          Then check if that point is inside the facet. If yes, return the time until collision with that facet.
     */

    double time_until_collision = 1e9;
    // Voxel index vector
    int i =  r[0]*system->one_over_length[0]*Nx;
    int j =  r[1]*system->one_over_length[1]*Ny;
    int k =  r[2]*system->one_over_length[2]*Nz;
    vector<vector<double> > cube;

    vector<double> p1(3,0); vector<double> p2(3,0); vector<double> p3(3,0); vector<double> p4(3,0);
    vector<double> p5(3,0); vector<double> p6(3,0); vector<double> p7(3,0); vector<double> p8(3,0);

    p1[0] = i*voxel_size[0]; p1[1] = j*voxel_size[1]; p1[2] = k*voxel_size[2];
    p2[0] = (i+1)*voxel_size[0]; p2[1] = j*voxel_size[1]; p2[2] = k*voxel_size[2];
    p3[0] = i*voxel_size[0]; p3[1] = j*voxel_size[1]; p3[2] = (k+1)*voxel_size[2];
    p4[0] = (i+1)*voxel_size[0]; p4[1] = j*voxel_size[1]; p4[2] = (k+1)*voxel_size[2];
    p5[0] = i*voxel_size[0]; p5[1] = (j+1)*voxel_size[1]; p5[2] = k*voxel_size[2];
    p6[0] = (i+1)*voxel_size[0]; p6[1] = (j+1)*voxel_size[1]; p6[2] = k*voxel_size[2];
    p7[0] = i*voxel_size[0]; p7[1] = (j+1)*voxel_size[1]; p7[2] = (k+1)*voxel_size[2];
    p8[0] = (i+1)*voxel_size[0]; p8[1] = (j+1)*voxel_size[1]; p8[2] = (k+1)*voxel_size[2];




}
