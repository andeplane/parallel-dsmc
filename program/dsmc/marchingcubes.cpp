#include <marchingcubes.h>
#include <progressbar.h>
#include <cmath>
#include <complexgeometry.h>
#include <random.h>
#include <cisosurface.h>

void MarchingCubes::save_to_file(string filename) {
//    ofstream file (filename.c_str(), ios::out | ios::binary);
//    file.write (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
//    file.write (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
//    file.write (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
//    file.write (reinterpret_cast<char*>(&vertices), 3*num_vertices*sizeof(float));
//    file.write (reinterpret_cast<char*>(&normals),   3*num_vertices*sizeof(float));
//    file.close();
}

void MarchingCubes::load_from_file(string filename) {
//    ifstream file (filename.c_str(), ios::in | ios::binary);
//    file.read (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
//    file.read (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
//    file.read (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
//    num_vertices = nx*ny*nz;

//    file.read (reinterpret_cast<char*>(&vertices), 3*num_vertices*sizeof(float));
//    file.read (reinterpret_cast<char*>(&normals),   3*num_vertices*sizeof(float));
//    file.close();
}
template <typename T>
void MarchingCubes::create_marching_cubes_from_array(const T* scalar_field, CVector mesh_dimensions, CVector box_length, double threshold, bool larger_than) {
    int nx = mesh_dimensions.x;
    int ny = mesh_dimensions.y;
    int nz = mesh_dimensions.z;
    initialize(nx*ny*nz);
    num_vertices = 0;
    num_triangles = 0;

    Random *rnd = new Random(-time(NULL), 0, 0);


    float r = 89/255.0;
    float g = 193/255.0;
    float b = 235/255.0;

//    float r = rnd->next_double();
//    float g = rnd->next_double();
//    float b = rnd->next_double();

    CIsoSurface<T> surf;
    surf.GenerateSurface(scalar_field,threshold,nx-1,ny-1,nz-1,box_length.x, box_length.y, box_length.z, larger_than);

    for(int triangle=0; triangle<surf.m_nTriangles; triangle++) {
        unsigned int p1_index = surf.m_piTriangleIndices[3*triangle+0];
        unsigned int p2_index = surf.m_piTriangleIndices[3*triangle+1];
        unsigned int p3_index = surf.m_piTriangleIndices[3*triangle+2];

        CVector p1(surf.m_ppt3dVertices[p1_index][0], surf.m_ppt3dVertices[p1_index][1], surf.m_ppt3dVertices[p1_index][2]);
        CVector p2(surf.m_ppt3dVertices[p2_index][0], surf.m_ppt3dVertices[p2_index][1], surf.m_ppt3dVertices[p2_index][2]);
        CVector p3(surf.m_ppt3dVertices[p3_index][0], surf.m_ppt3dVertices[p3_index][1], surf.m_ppt3dVertices[p3_index][2]);
        VECTOR3D &norm1 = surf.m_pvec3dNormals[p1_index];
        VECTOR3D &norm2 = surf.m_pvec3dNormals[p2_index];
        VECTOR3D &norm3 = surf.m_pvec3dNormals[p3_index];

        CVector n1(norm1[0], norm1[1], norm1[2]);
        CVector n2(norm2[0], norm2[1], norm2[2]);
        CVector n3(norm3[0], norm3[1], norm3[2]);
//        int i = p1.x;
//        int j = p1.z;
//        int k = p1.y;
//        int index = i + j*cg.nx + k*cg.nx*cg.ny;
//        float relative_value = 0.5 + 0.5*cg.vertices_float[index]/cg.max_value;

        // CVector color(r+0.7*relative_value, g+0.8*relative_value, b+0.9*relative_value);
        CVector color(r,g,b);

        add_vertex(p1);
        add_normal(n1);
        add_color(color, 1.0);

        add_vertex(p2);
        add_normal(n2);
        add_color(color, 1.0);

        add_vertex(p3);
        add_normal(n3);
        add_color(color, 1.0);
    }

    cout << "Marching cubes created with " << surf.m_nTriangles << " triangles." << endl;
}

void MarchingCubes::create_marching_cubes_from_complex_geometry(ComplexGeometry &cg, CVector system_length, double threshold, bool larger_than) {
    CVector box_length = system_length;
    box_length.x /= cg.nx;
    box_length.y /= cg.ny;
    box_length.z /= cg.nz;
//    float *field = new float[cg.num_vertices];
//    for(int n=0; n<cg.num_vertices; n++) {
//        field[n] = cg.vertices_unsigned_char[n];
//    }
    CVector mesh_dimensions(cg.nx, cg.ny, cg.nz);
    int num_vertices = cg.nx*cg.ny*cg.nz;
    float *vertices = new float[num_vertices];
    for(int i=0; i<cg.nx; i++) {
        for(int j=0; j<cg.ny; j++) {
            for(int k=0; k<cg.nz; k++) {
                int idx1 = i*cg.nz*cg.ny + j*cg.nz + k;
                int idx2 = i + j*cg.nx + k*cg.nx*cg.ny;
                vertices[idx2] = cg.vertices[idx1];
            }
        }
    }

    create_marching_cubes_from_array(vertices, mesh_dimensions, box_length, threshold, larger_than);
}

void MarchingCubes::create_marching_cubes_from_positions(vector<float> &positions, CVector num_points, CVector system_length, double threshold) {
    int total_num_vertices = num_points.x*num_points.y*num_points.z;
    float *number_of_particles_list = new float[total_num_vertices];
    int num_particles = positions.size() / 3;
    int nx = num_points.x; int ny = num_points.y; int nz = num_points.z;
    float max_x = 0;

    for(int n=0; n<num_particles; n++) {
        float x = positions[3*n+0];
        float y = positions[3*n+1];
        float z = positions[3*n+2];
        max_x = max(max_x,x);

        int i = x / system_length.x*num_points.x;
        int j = y / system_length.y*num_points.y;
        int k = z / system_length.z*num_points.z;

        int delta = 1;
        for(int di=-delta; di<=delta; di++) {
            for(int dj=-delta; dj<=delta; dj++) {
                for(int dk=-delta; dk<=delta; dk++) {
                    int index = ((i+di+nx)%nx) + ((j+dj+ny)%ny)*nx+ ((k+dk+nz)%nz)*nx*ny;
                    number_of_particles_list[index] += 1;
                }
            }
        }
    }

    CVector box_length = system_length;
    box_length = box_length/num_points;

    create_marching_cubes_from_array(number_of_particles_list, num_points, box_length, threshold, true);
}
