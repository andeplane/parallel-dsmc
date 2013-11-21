#pragma once
#include <system.h>
#include <vector>
#include <string>
#include <cvector.h>
#include <topology.h>
class DSMC_IO;

using namespace std;

typedef enum {
    voxel_type_empty = 0,
    voxel_type_wall = 1,
    voxel_type_boundary = 2
} voxel_type;

inline double time_until_collision_with_plane(CVector &r, CVector &v, CVector &point_in_plane, CVector &normal) {
    return (point_in_plane - r).dot(normal) / v.dot(normal);
}

inline bool same_side(CVector &p1, CVector &p2, CVector &a, CVector &b) {
    CVector cp1 = (b-a).cross(p1-a);
    CVector cp2 = (b-a).cross(p2-a);
    if(cp1.dot(cp2) >= 0) return true;
    else return false;
}

inline bool is_point_within_square(CVector &p1, CVector &p2, CVector &p3, CVector &p4, CVector point) {
    CVector center = (p1 + p2 + p3 + p4)*0.25;

    if (same_side(point, center, p1, p2) && same_side(point, center, p2, p3) &&
        same_side(point, center, p3, p4) && same_side(point, center, p4, p1)) return true;
    else return false;
}

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
    unsigned char *get_voxel(vector<data_type> &r, int index);
    void get_index_vector_from_index(const int &index, int &i, int &j, int &k);
    int get_index_of_voxel(const data_type *r);
    double get_time_until_collision(data_type *r, data_type *v, double dt, const int &voxel_index);

    void read_matrix(string filename, DSMC_IO *io);
};

inline double Grid::get_time_until_collision(data_type *r, data_type *v, double dt, const int &voxel_index) {
    /*
     * Strategy:
     *          First calculate time until intersection with every facet of the voxel.
     *          Then, for each facet, check if collision point is inside the facet. Among those points that are inside a facet, choose the one with lowest collision time.
     */

    CVector r_vec(r[0]-system->topology->origin[0], r[1]-system->topology->origin[1], r[2]-system->topology->origin[2]);
    CVector v_vec(v[0], v[1], v[2]);

    double time_until_collision = 1e9;
    // Voxel index vector
    int i,j,k;
    get_index_vector_from_index(voxel_index, i, j, k);
    i -= nx_per_cpu;
    j -= ny_per_cpu;
    k -= nz_per_cpu;

    point_list[0].x = i*voxel_size[0]; point_list[0].y = j*voxel_size[1]; point_list[0].z = k*voxel_size[2];
    point_list[1].x = (i+1)*voxel_size[0]; point_list[1].y = j*voxel_size[1]; point_list[1].z = k*voxel_size[2];
    point_list[2].x = i*voxel_size[0]; point_list[2].y = j*voxel_size[1]; point_list[2].z = (k+1)*voxel_size[2];
    point_list[3].x = (i+1)*voxel_size[0]; point_list[3].y = j*voxel_size[1]; point_list[3].z = (k+1)*voxel_size[2];
    point_list[4].x = i*voxel_size[0]; point_list[4].y = (j+1)*voxel_size[1]; point_list[4].z = k*voxel_size[2];
    point_list[5].x = (i+1)*voxel_size[0]; point_list[5].y = (j+1)*voxel_size[1]; point_list[5].z = k*voxel_size[2];
    point_list[6].x = i*voxel_size[0]; point_list[6].y = (j+1)*voxel_size[1]; point_list[6].z = (k+1)*voxel_size[2];
    point_list[7].x = (i+1)*voxel_size[0]; point_list[7].y = (j+1)*voxel_size[1]; point_list[7].z = (k+1)*voxel_size[2];

    double time_facet_1 = time_until_collision_with_plane(r_vec, v_vec, point_list[2], unit_normal_vectors[0]);
    double time_facet_2 = time_until_collision_with_plane(r_vec, v_vec, point_list[1], unit_normal_vectors[1]);
    double time_facet_3 = time_until_collision_with_plane(r_vec, v_vec, point_list[0], unit_normal_vectors[2]);
    double time_facet_4 = time_until_collision_with_plane(r_vec, v_vec, point_list[4], unit_normal_vectors[3]);
    double time_facet_5 = time_until_collision_with_plane(r_vec, v_vec, point_list[0], unit_normal_vectors[4]);
    double time_facet_6 = time_until_collision_with_plane(r_vec, v_vec, point_list[2], unit_normal_vectors[5]);

    bool will_hit_facet_1 = is_point_within_square(point_list[4], point_list[0], point_list[2], point_list[6], r_vec+v_vec*time_facet_1);
    bool will_hit_facet_2 = is_point_within_square(point_list[1], point_list[5], point_list[7], point_list[3], r_vec+v_vec*time_facet_2);
    bool will_hit_facet_3 = is_point_within_square(point_list[0], point_list[1], point_list[3], point_list[2], r_vec+v_vec*time_facet_3);
    bool will_hit_facet_4 = is_point_within_square(point_list[5], point_list[4], point_list[6], point_list[7], r_vec+v_vec*time_facet_4);
    bool will_hit_facet_5 = is_point_within_square(point_list[4], point_list[5], point_list[1], point_list[0], r_vec+v_vec*time_facet_5);
    bool will_hit_facet_6 = is_point_within_square(point_list[2], point_list[3], point_list[7], point_list[6], r_vec+v_vec*time_facet_6);

    if( time_facet_1 > 0 && time_facet_1 < time_until_collision && will_hit_facet_1 && !std::isnan(time_facet_1)) time_until_collision = time_facet_1;
    if( time_facet_2 > 0 && time_facet_2 < time_until_collision && will_hit_facet_2 && !std::isnan(time_facet_2)) time_until_collision = time_facet_2;
    if( time_facet_3 > 0 && time_facet_3 < time_until_collision && will_hit_facet_3 && !std::isnan(time_facet_3)) time_until_collision = time_facet_3;
    if( time_facet_4 > 0 && time_facet_4 < time_until_collision && will_hit_facet_4 && !std::isnan(time_facet_4)) time_until_collision = time_facet_4;
    if( time_facet_5 > 0 && time_facet_5 < time_until_collision && will_hit_facet_5 && !std::isnan(time_facet_5)) time_until_collision = time_facet_5;
    if( time_facet_6 > 0 && time_facet_6 < time_until_collision && will_hit_facet_6 && !std::isnan(time_facet_6)) time_until_collision = time_facet_6;

    if(time_until_collision > 1000) {
        cout << system->myid << " didn't collide with anything :/" << endl;
        cout << "Time 1: " << time_facet_1 << endl;
        cout << "Time 2: " << time_facet_2 << endl;
        cout << "Time 3: " << time_facet_3 << endl;
        cout << "Time 4: " << time_facet_4 << endl;
        cout << "Time 5: " << time_facet_5 << endl;
        cout << "Time 6: " << time_facet_6 << endl;
        cout << "Particle at:" << endl;
        cout << "r=[" << r[0] << " " << r[1] << " " << r[2] << "]" << endl;
        cout << "v=[" << v[0] << " " << v[1] << " " << v[2] << "]" << endl;
        cout << "with local i,j,k = " << i << ", " << j << ", " << k << ", " << endl;
        cout << "which is local voxel index " << voxel_index << endl;

        exit(1);
    }

    if(time_until_collision < 1e-7) return 0;
    else return time_until_collision - 1e-7; // Subtract a small number to avoid being exactly at the boundary of a voxel
}

inline int Grid::get_index_of_voxel(const data_type *r) {
    const int i =  (r[0]-system->topology->origin[0])*system->one_over_length[0]*global_nx + voxel_origin.x;
    const int j =  (r[1]-system->topology->origin[1])*system->one_over_length[1]*global_ny + voxel_origin.y;
    const int k =  (r[2]-system->topology->origin[2])*system->one_over_length[2]*global_nz + voxel_origin.z;

    return i*ny*nz + j*nz + k;
}
