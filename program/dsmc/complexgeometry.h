#pragma once
#include <string>
#include <vector>
#include <cvector.h>

using std::string;
using std::vector;

typedef enum {
    voxel_type_empty = 0,
    voxel_type_wall = 1,
    voxel_type_boundary = 2
} voxel_type;

class ComplexGeometry
{
public:
    float *vertices;
    float *normals;
    float *tangents1;
    float *tangents2;
    unsigned int nx, ny, nz, num_vertices;

    // Used to save binary matrix for dsmc
    unsigned char *vertices_unsigned_char;
    float global_porosity;
    CVector num_voxels;
    bool has_normals_tangents_and_boundary;

    float max_value;

    ComplexGeometry();
    void allocate(int nx_, int ny_, int nz_);
    void load_from_binary_file_without_normals_and_tangents(string filename, bool do_calculate_normals_tangents_and_inner_points, int number_of_neighbor_averages);
    void load_text_files(string base_filename, CVector matrix_size, double threshold);
    void create_perlin_geometry(int nx, int ny, int nz, int octave, int frequency, int amplitude , int seed, float threshold, bool do_calculate_normals_tangents_and_boundary, int number_of_neighbor_averages);
    void create_sphere(int nx_, int ny_, int nz_, float radius, bool inverse, bool do_calculate_normals_tangents_and_inner_points, int number_of_neighbor_averages);
    void create_empty_space(int nx_, int ny_, int nz_, bool do_calculate_normals_tangents_and_inner_points, int number_of_neighbor_averages);
    void calculate_normals_tangents_and_inner_points(int number_of_neighbor_averages);
    void calculate_global_porosity();
    void find_boundary_points();
    void calculate_tangents();
    void calculate_normals(int number_of_neighbor_average);
    // void save_to_file(string filename);
    void save_to_file(string foldername, CVector num_processors_vector);
    void save_vtk(string filename);
    void load_vtk(string filename);
};
