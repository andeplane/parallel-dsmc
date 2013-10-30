#include <complexgeometry.h>
#include <perlin.h>
#include <progressbar.h>
#include <fstream>
#include <cmath>
#include <random.h>
#include <fstream>
#include <cvector.h>
#include <cstring>
#include <cinifile.h>
#include <diamondSquare.h>

using std::ifstream;
using std::ofstream;
using std::ios;
using std::fstream;

ComplexGeometry::ComplexGeometry()
{
    vertices = NULL;
    normals = NULL;
    tangents1 = NULL;
    tangents2 = NULL;
    nx = 0; ny = 0; nz = 0; num_vertices = 0;
    num_voxels = CVector(0,0,0);
    has_normals_tangents_and_boundary = false;
    max_value = 0;
    global_porosity = -1;
}

void ComplexGeometry::allocate(int nx_, int ny_, int nz_) {
    nx = nx_; ny = ny_; nz = nz_;
    num_vertices = nx*ny*nz;
    if(vertices != NULL) {
        delete vertices;
        delete normals;
        delete tangents1;
        delete tangents2;
        delete vertices_unsigned_char;
    }
    num_voxels = CVector(nx,ny,nz);
    vertices_unsigned_char = new unsigned char[num_vertices];
    vertices = new float[num_vertices];
    normals = new float[3*num_vertices];
    tangents1 = new float[3*num_vertices];
    tangents2 = new float[3*num_vertices];
}

void ComplexGeometry::load_text_files(string base_filename, CVector matrix_size, double threshold) {
    allocate(matrix_size.x, matrix_size.y, matrix_size.z);
    int file_start = 1;
    char filename[1000];
    int num_vertices_so_far = 0;
    for(int file_num = file_start; file_num<file_start+nz; file_num++) {
        sprintf(filename, "%s%02d.txt",base_filename.c_str(), file_num);
        ifstream infile(filename);
        for(int line_num=0; line_num<nx; line_num++) {
            for(int i=0; i<ny; i++) {
                float value;
                infile >> value;
                max_value = max(max_value,value);
                vertices[num_vertices_so_far] = value;
                vertices_unsigned_char[num_vertices_so_far++] = value >= threshold;
            }
        }
        infile.close();
    }
}

void ComplexGeometry::load_vtk(string filename) {
    cout << "VTK loader is not implemented yet" << endl;
    exit(1);

//    string line;
//    ifstream file (filename.c_str());

//    string re1="(DIMENSIONS)";	// Word 1
//    string re2=".*?";	// Non-greedy match on filler
//    string re3="\\d+";	// Uninteresting: int
//    string re4=".*?";	// Non-greedy match on filler
//    string re5="\\d+";	// Uninteresting: int
//    string re6=".*?";	// Non-greedy match on filler
//    string re7="(\\d+)";	// Integer Number 1

//    PME re(re1+re2+re3+re4+re5+re6+re7,"gims");
//    int n;

//    if(file.is_open()) {
//        while(getline(file,line)) {
//            if ((n=re.match(txt))>0)
//            {
//                string word1=re[1].c_str();
//                string int1=re[2].c_str();
//                cout << "("<<word1<<")"<<"("<<int1<<")"<< std::endl;
//            }
//        }
//    }

////    ofile << "# vtk DataFile Version 2.0" << endl;
////    ofile << "structured point" << endl;
////    ofile << "ASCII" << endl;
////    ofile << endl;
////    ofile << "DATASET STRUCTURED_POINTS" << endl;
////    ofile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
////    ofile << "ORIGIN 0.0 0.0 0.0" << endl;
////    // ofile << "SPACING 1 1 1" << endl;
////    ofile << "SPACING " << 1.0/double(nx) << " " << 1.0/double(ny) << " " << 1.0/double(nz) << endl;
////    ofile << "POINT_DATA " << N << endl;
////    ofile << "SCALARS atomdist double" << endl;
////    ofile << "LOOKUP_TABLE default" << endl;
////    ofile << endl;

////    // column-major ordering...
////    for (int k = 0; k < nz; k++) {
////        for (int j = 0; j < ny; j++) {
////            for (int i = 0; i < nx; i++) {
////                int index = i + j*nx + k*nx*ny;
////                ofile << vertices[index] << endl;
////            }
////        }
////    }

//    file.close();
}

void ComplexGeometry::save_vtk(string filename) {
    ofstream ofile(filename.c_str());

    int N = nx*ny*nz;

    ofile << "# vtk DataFile Version 2.0" << endl;
    ofile << "structured point" << endl;
    ofile << "ASCII" << endl;
    ofile << endl;
    ofile << "DATASET STRUCTURED_POINTS" << endl;
    ofile << "DIMENSIONS " << nx << " " << ny << " " << nz << endl;
    ofile << "ORIGIN 0.0 0.0 0.0" << endl;
    // ofile << "SPACING 1 1 1" << endl;
    ofile << "SPACING " << 1.0/double(nx) << " " << 1.0/double(ny) << " " << 1.0/double(nz) << endl;
    ofile << "POINT_DATA " << N << endl;
    ofile << "SCALARS atomdist double" << endl;
    ofile << "LOOKUP_TABLE default" << endl;
    ofile << endl;

    // column-major ordering...
    for (int k = 0; k < nz; k++) {
        for (int j = 0; j < ny; j++) {
            for (int i = 0; i < nx; i++) {
                int index = i + j*nx + k*nx*ny;
                ofile << vertices[index] << endl;
            }
        }
    }

    ofile.close();
}

void ComplexGeometry::calculate_global_porosity() {
    // Porosity voxels are those with type voxel_type_empty
    int sum_wall_voxels = 0;
    for(int i=0; i<num_vertices; i++) sum_wall_voxels += vertices_unsigned_char[i] == voxel_type_empty;

    global_porosity = (float)sum_wall_voxels / (nx*ny*nz);
}

void ComplexGeometry::load_from_binary_file_without_normals_and_tangents(string filename, bool do_calculate_normals_tangents_and_inner_points, int number_of_neighbor_averages) {
    ifstream file (filename.c_str(), ios::in | ios::binary);
    file.read (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
    file.read (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
    file.read (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
    allocate(nx, ny, nz);

    file.read (reinterpret_cast<char*>(vertices_unsigned_char), num_vertices*sizeof(unsigned char));
    file.close();

    for(int i=0; i<num_vertices; i++) vertices[i] = vertices_unsigned_char[i];

    if(do_calculate_normals_tangents_and_inner_points) calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::calculate_normals_tangents_and_inner_points(int number_of_neighbor_averages) {
    calculate_normals(number_of_neighbor_averages);
    calculate_tangents();
    find_boundary_points();
    has_normals_tangents_and_boundary = true;
}

CVector index_vector_from_index(const int &index, CVector num_processors) {
    CVector index_vector(0,0,0);
    int num_proc_x = num_processors.x;
    int num_proc_y = num_processors.y;
    int num_proc_z = num_processors.z;

    index_vector.x = index/(num_proc_y*num_proc_z);   // Index in x-direction
    index_vector.y = (index/num_proc_z)%num_proc_y;   // Index in y-direction
    index_vector.z = index%num_proc_z;              // Index in z-direction

    return index_vector;
}

int index_from_ijk(CVector voxel_index_vector, CVector num_voxels) {
    return voxel_index_vector.x*num_voxels.y*num_voxels.z + voxel_index_vector.y*num_voxels.z + voxel_index_vector.z;
}

void ComplexGeometry::save_to_file(CIniFile &ini) {
    string foldername = ini.getstring("binary_output_folder");
    CVector num_processors_vector(ini.getint("num_processors_x"), ini.getint("num_processors_y"), ini.getint("num_processors_z"));

    int num_processors_total = num_processors_vector.x*num_processors_vector.y*num_processors_vector.z;
    CVector num_voxels_vector_per_node = CVector(nx,ny,nz)/num_processors_vector;
    int num_voxels_total_per_node = 3*3*3*num_voxels_vector_per_node.x*num_voxels_vector_per_node.y*num_voxels_vector_per_node.z;

    calculate_global_porosity();
    char filename[1000];
    unsigned char *local_voxels = new unsigned char[num_voxels_total_per_node];
    float *local_normals = new float[3*num_voxels_total_per_node];
    float *local_tangents1 = new float[3*num_voxels_total_per_node];
    float *local_tangents2 = new float[3*num_voxels_total_per_node];

    cout << "Saving geometry to " << foldername << " on cpus=(" << num_processors_vector.x << ", " << num_processors_vector.y << ", " << num_processors_vector.z << ")" << endl;
    ProgressBar p(num_processors_total,"Saving world files");
    for(int index = 0; index<num_processors_total; index++) {
        sprintf(filename,"%s/%04d.bin",foldername.c_str(), index);
        ofstream file (filename, ios::out | ios::binary);
        if(!file.is_open()) {
            cout << "Could not open " << filename << ". Check your permissions and that the folder exists." << endl;
            exit(1);
        }

        CVector index_vector = index_vector_from_index(index, num_processors_vector);
        CVector voxel_origo = index_vector*num_voxels_vector_per_node;

        int output_data_array_index = 0;
        int num_empty_voxels = 0;
        int num_local_voxels = 0;
        for(int di = -num_voxels_vector_per_node.x; di < 2*num_voxels_vector_per_node.x; di++) {
            for(int dj = -num_voxels_vector_per_node.y; dj < 2*num_voxels_vector_per_node.y; dj++) {
                for(int dk = -num_voxels_vector_per_node.z; dk < 2*num_voxels_vector_per_node.z; dk++) {
                    int i = int(voxel_origo.x + di + 10*nx) % nx;
                    int j = int(voxel_origo.y + dj + 10*ny) % ny;
                    int k = int(voxel_origo.z + dk + 10*nz) % nz;

                    CVector global_voxel_index_vector(i,j,k);

                    int global_voxel_index = index_from_ijk(global_voxel_index_vector, num_voxels);

                    local_voxels[output_data_array_index]  = vertices_unsigned_char[global_voxel_index];

                    for(int a=0; a<3; a++) {
                        local_normals[3*output_data_array_index + a] = normals[3*global_voxel_index + a];
                        local_tangents1[3*output_data_array_index + a] = tangents1[3*global_voxel_index + a];
                        local_tangents2[3*output_data_array_index + a] = tangents2[3*global_voxel_index + a];
                    }

                    // Calculate porosity of this node
                    bool is_local_node = (di>=0 && di < num_voxels_vector_per_node.x) && (dj>=0 && dj < num_voxels_vector_per_node.y) && (dk>=0 && dk < num_voxels_vector_per_node.z);

                    if(is_local_node) {
                        num_empty_voxels += vertices_unsigned_char[global_voxel_index]==voxel_type_empty;
                        num_local_voxels++;
                    }

                    output_data_array_index++;
                }
            }
        }

        float porosity = 0;
        if(num_local_voxels > 0) porosity = (float)num_empty_voxels / num_local_voxels;
        // Each dimension will have the neighbor node's voxel as well
        unsigned int local_nx = 3*num_voxels_vector_per_node.x;
        unsigned int local_ny = 3*num_voxels_vector_per_node.y;
        unsigned int local_nz = 3*num_voxels_vector_per_node.z;
        file.write (reinterpret_cast<char*>(&global_porosity), sizeof(float));
        file.write (reinterpret_cast<char*>(&porosity), sizeof(float));
        file.write (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
        file.write (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
        file.write (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
        file.write (reinterpret_cast<char*>(&local_nx), sizeof(unsigned int));
        file.write (reinterpret_cast<char*>(&local_ny), sizeof(unsigned int));
        file.write (reinterpret_cast<char*>(&local_nz), sizeof(unsigned int));

        file.write (reinterpret_cast<char*>(local_voxels), output_data_array_index*sizeof(unsigned char));
        file.write (reinterpret_cast<char*>(local_normals),   3*output_data_array_index*sizeof(float));
        file.write (reinterpret_cast<char*>(local_tangents1), 3*output_data_array_index*sizeof(float));
        file.write (reinterpret_cast<char*>(local_tangents2), 3*output_data_array_index*sizeof(float));

        file.close();

        p.update(index);
    }
}

void ComplexGeometry::create_sphere(CIniFile &ini) {
    int nx_ = ini.getint("num_voxels_x");
    int ny_ = ini.getint("num_voxels_y");
    int nz_ = ini.getint("num_voxels_z");
    float radius = ini.getdouble("sphere_radius");
    bool inverse = ini.getbool("sphere_inverse");
    int number_of_neighbor_averages = ini.getint("number_of_neighbor_averages");
    allocate(nx_, ny_, nz_);

    cout << "Creating sphere with radius=" << radius<< ", inverse=" << inverse << " on num_voxels=(" << nx << ", " << ny << ", " << nz << ")." << endl;
    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                double x = 2*(i-nx/2.0)/(double)nx;
                double y = 2*(j-ny/2.0)/(double)ny;
                double z = 2*(k-nz/2.0)/(double)nz;
                double r2 = x*x + y*y + z*z;
                int index = i + j*nx + k*nx*ny;

                if( (!inverse && r2 > radius*radius) || (inverse && r2 < radius*radius)) {
                    vertices_unsigned_char[index] = 1;
                    vertices[index] = 1;
                } else {
                    vertices_unsigned_char[index] = 0;
                    vertices[index] = 0;
                }
            }
        }
    }
    calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::create_diamond_square(CIniFile &ini) {
    int nx_ = ini.getint("num_voxels_x");
    int ny_ = ini.getint("num_voxels_y");
    int nz_ = ini.getint("num_voxels_z");
    float hurst_exponent = ini.getdouble("hurst_exponent");
    long seed = ini.getdouble("perlin_seed");
    int number_of_neighbor_averages = ini.getint("number_of_neighbor_averages");
    float h0 = ini.getdouble("diamond_square_h0");
    allocate(nx_, ny_, nz_);

    if(nx != nz) {
        cout << "Warning, nx != ny, this will cause problems with the diamond square algorithm" << endl;
    }

    float power2_x = log2(nx);
    float power2_z = log2(nz);

    if(pow(2.0f,power2_x) != nx) cout << "Warning, nx is not a power of 2" << endl;
    if(pow(2.0f,power2_z) != nz) cout << "Warning, nx is not a power of 2" << endl;

    cout << "Creating diamond square with hurst exponent " << hurst_exponent << " on num_voxels=(" << nx << ", " << ny << ", " << nz << ")." << endl;

    float sigma = 0.3;
    long seed_1 = -(seed + 1);
    long seed_2 = -(seed + 2);
    bool addition = true;
    bool periodic_boundary_conditions = true;
    int rng = 2;
    vector<double> corners(4,h0);

    DiamondSquare generator;
    vector<vector<double> > height_map_1 = generator.generate(power2_x, hurst_exponent, corners, seed_1, sigma, addition, periodic_boundary_conditions, rng);
    vector<vector<double> > height_map_2 = generator.generate(power2_x, hurst_exponent, corners, seed_2, sigma, addition, periodic_boundary_conditions, rng);

    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int index = i + j*nx + k*nx*ny;

                float h1 = 5*float(j) / ny;
                float h2 = 5*float(ny - j) / ny;

                if( h1 < height_map_1[k][i] || h2 < height_map_2[k][i] ) {
                    vertices_unsigned_char[index] = 1;
                    vertices[index] = 1;
                } else {
                    vertices_unsigned_char[index] = 0;
                    vertices[index] = 0;
                }

                if(j == 0 || j == ny-1) {
                    vertices_unsigned_char[index] = 1;
                    vertices[index] = 1;
                }
            }
        }
    }
    calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::create_empty_space(CIniFile &ini) {
    int nx_ = ini.getint("num_voxels_x");
    int ny_ = ini.getint("num_voxels_y");
    int nz_ = ini.getint("num_voxels_z");
    int number_of_neighbor_averages = ini.getint("number_of_neighbor_averages");
    allocate(nx_, ny_, nz_);
    cout << "Creating empty space on num_voxels=(" << nx << ", " << ny << ", " << nz << ")." << endl;
    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int index = i + j*nx + k*nx*ny;
                vertices_unsigned_char[index] = 0;
                vertices[index] = 0;
            }
        }
    }

    calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::create_box(CIniFile &ini) {
    int nx_ = ini.getint("num_voxels_x");
    int ny_ = ini.getint("num_voxels_y");
    int nz_ = ini.getint("num_voxels_z");
    int number_of_neighbor_averages = ini.getint("number_of_neighbor_averages");
    float porosity = ini.getdouble("box_porosity");
    allocate(nx_, ny_, nz_);
    cout << "Creating box with porosity=" << porosity << " on num_voxels=(" << nx << ", " << ny << ", " << nz << ")." << endl;
    int mid_y = ny/2;
    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int index = i + j*nx + k*nx*ny;
                int dy = abs(j-mid_y);
                if(dy > mid_y*porosity) {
                    vertices_unsigned_char[index] = 1;
                    vertices[index] = 1;
                } else {
                    vertices_unsigned_char[index] = 0;
                    vertices[index] = 0;
                }
            }
        }
    }

    calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::create_poiseuille(CIniFile &ini) {
    int nx_ = ini.getint("num_voxels_x");
    int ny_ = ini.getint("num_voxels_y");
    int nz_ = ini.getint("num_voxels_z");
    int number_of_neighbor_averages = ini.getint("number_of_neighbor_averages");

    allocate(nx_, ny_, nz_);
    cout << "Creating Poiseuille plates on num_voxels=(" << nx << ", " << ny << ", " << nz << ")." << endl;
    int mid_y = ny/2;
    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int index = i + j*nx + k*nx*ny;

                if(j == 0 || j == ny-1) {
                    vertices_unsigned_char[index] = 1;
                    vertices[index] = 1;
                } else {
                    vertices_unsigned_char[index] = 0;
                    vertices[index] = 0;
                }
            }
        }
    }

    calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::create_perlin_geometry(CIniFile &ini) {
    int nx_ = ini.getint("num_voxels_x");
    int ny_ = ini.getint("num_voxels_y");
    int nz_ = ini.getint("num_voxels_z");
    int number_of_neighbor_averages = ini.getint("number_of_neighbor_averages");
    int octave = ini.getint("perlin_octave");
    int frequency = ini.getint("perlin_frequency");
    int amplitude = ini.getint("perlin_amplitude");
    int seed = ini.getint("perlin_seed");
    int num_scales = ini.getint("perlin_num_scales");
    int constant = ini.getdouble("perlin_constant");
    float scale_factor = ini.getdouble("perlin_scale_factor");
    float threshold = ini.getdouble("perlin_threshold");

    Perlin p(octave, frequency, amplitude, seed);
    allocate(nx_, ny_, nz_);

    cout << "Creating perlin noise with " << endl;
    cout << "octave=" << octave << endl;
    cout << "frequency=" << frequency << endl;
    cout << "amplitude=" << amplitude << endl;
    cout << "seed=" << seed << endl;
    cout << "num_scales=" << num_scales << endl;
    cout << "constant=" << constant << endl;
    cout << "scale_factor=" << scale_factor << endl;
    cout << "threshold=" << threshold << endl;
    cout << " on num_voxels=(" << nx << ", " << ny << ", " << nz << ")." << endl;

    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                double x = (i-nx/2.0)/(double)nx;
                double y = (j-ny/2.0)/(double)ny;
                double z = (k-nz/2.0)/(double)nz;
                float s = 1.0;

                int index = i*ny*nz + j*nz + k; //new

                double val = 0;
                for (int a=0; a<num_scales  ; a++) {
                    s = scale_factor*a + constant;

                    val += p.Get(x*s, y*s, z*s);
                }

                vertices[index] = val;
                if(val >= threshold) vertices_unsigned_char[index] = 1;
                else vertices_unsigned_char[index] = 0;
            }
        }
    }

    calculate_normals_tangents_and_inner_points(number_of_neighbor_averages);
}

void ComplexGeometry::calculate_normals(int number_of_neighbor_average) {
    bool at_least_one_wall_neighbor;
    bool all_neighbors_are_walls;
    memset(normals, 0, 3*num_vertices*sizeof(float));
    memset(tangents1, 0, 3*num_vertices*sizeof(float));
    memset(tangents2, 0, 3*num_vertices*sizeof(float));

    ProgressBar progress_bar(nx, "Creating normals");
    for(int i=0;i<nx;i++) {
        progress_bar.update(i);

        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                at_least_one_wall_neighbor = false;
                all_neighbors_are_walls = true;

                // int idx = i + j*nx + k*nx*ny;
                int idx = i*ny*nz + j*nz + k; //new
                for(int di=-1;di<=1;di++) {
                    for(int dj=-1;dj<=1;dj++) {
                        for(int dk=-1;dk<=1;dk++) {
                            int idx2 = ((i+di+nx)%nx)*ny*nz + ((j+dj+ny)%ny)*nz + ((k+dk+nz) % nz); //new

                            if(vertices_unsigned_char[idx2]>0) {
                                // If at least one wall neighbor, this is not a single wall voxel
                                at_least_one_wall_neighbor = true;
                            }
                            if(vertices_unsigned_char[idx2]==0) {
                                // If not all neighbors are walls and we have norm=1, this is a
                                // single plane that has no defined normal vector.
                                all_neighbors_are_walls = false;
                            }

                            normals[3*idx+0] -= vertices_unsigned_char[idx2]*di;
                            normals[3*idx+1] -= vertices_unsigned_char[idx2]*dj;
                            normals[3*idx+2] -= vertices_unsigned_char[idx2]*dk;
                        }
                    }
                }

                double norm_squared = normals[3*idx+0]*normals[3*idx+0] + normals[3*idx+1]*normals[3*idx+1] + normals[3*idx+2]*normals[3*idx+2];
                if(norm_squared > 0) {
                    normals[3*idx+0] = 0;
                    normals[3*idx+1] = 0;
                    normals[3*idx+2] = 0;
                    for(int di=-number_of_neighbor_average; di<=number_of_neighbor_average; di++) {
                        for(int dj=-number_of_neighbor_average; dj<=number_of_neighbor_average; dj++) {
                            for(int dk=-number_of_neighbor_average; dk<=number_of_neighbor_average; dk++) {
                                int idx2 = ((i+di+nx)%nx)*ny*nz + ((j+dj+ny)%ny)*nz + ((k+dk+nz) % nz); //new

                                normals[3*idx+0] -= vertices_unsigned_char[idx2]*di;
                                normals[3*idx+1] -= vertices_unsigned_char[idx2]*dj;
                                normals[3*idx+2] -= vertices_unsigned_char[idx2]*dk;
                            }
                        }
                    }
                }


                double norm = sqrt(normals[3*idx+0]*normals[3*idx+0] + normals[3*idx+1]*normals[3*idx+1] + normals[3*idx+2]*normals[3*idx+2]);

                if(norm > 0) {
                    normals[3*idx+0] /= norm;
                    normals[3*idx+1] /= norm;
                    normals[3*idx+2] /= norm;
                } else if(!at_least_one_wall_neighbor || !all_neighbors_are_walls) {
                    // Single point or single pixel-plane, should not be a wall
                    vertices_unsigned_char[idx] = 0;
                }
            }
        }
    }
}

void ComplexGeometry::calculate_tangents() {
    Random *rnd = new Random(-1,0,0);
    ProgressBar progress_bar(nx, "Creating tangent vectors");
    for(int i=0;i<nx;i++) {
        progress_bar.update(i);
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int idx = i*ny*nz + j*nz + k; //new

                tangents1[3*idx+0] = rnd->next_double();
                tangents1[3*idx+1] = rnd->next_double();
                tangents1[3*idx+2] = rnd->next_double();

                double dot_product = normals[3*idx+0]*tangents1[3*idx+0] + normals[3*idx+1]*tangents1[3*idx+1] + normals[3*idx+2]*tangents1[3*idx+2];

                // Perform gram-schmidt
                tangents1[3*idx+0] -= normals[3*idx+0]*dot_product;
                tangents1[3*idx+1] -= normals[3*idx+1]*dot_product;
                tangents1[3*idx+2] -= normals[3*idx+2]*dot_product;

                // Normalize
                double norm = sqrt(tangents1[3*idx+0]*tangents1[3*idx+0] + tangents1[3*idx+1]*tangents1[3*idx+1] + tangents1[3*idx+2]*tangents1[3*idx+2]);

                if(norm>0) {
                    tangents1[3*idx+0] /= norm;
                    tangents1[3*idx+1] /= norm;
                    tangents1[3*idx+2] /= norm;
                }

                // t2 = n x t1
                tangents2[3*idx+0] = tangents1[3*idx+1]*normals[3*idx+2] - tangents1[3*idx+2]*normals[3*idx+1];
                tangents2[3*idx+1] = tangents1[3*idx+2]*normals[3*idx+0] - tangents1[3*idx+0]*normals[3*idx+2];
                tangents2[3*idx+2] = tangents1[3*idx+0]*normals[3*idx+1] - tangents1[3*idx+1]*normals[3*idx+0];

                // Normalize
                norm = sqrt(tangents2[3*idx+0]*tangents2[3*idx+0] + tangents2[3*idx+1]*tangents2[3*idx+1] + tangents2[3*idx+2]*tangents2[3*idx+2]);

                if(norm>0) {
                    tangents2[3*idx+0] /= norm;
                    tangents2[3*idx+1] /= norm;
                    tangents2[3*idx+2] /= norm;
                }
            }
        }
    }
}

void ComplexGeometry::find_boundary_points() {
    ProgressBar progress_bar(nx, "Finding boundary points");
    for(int i=0;i<nx;i++) {
        progress_bar.update(i);
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                int idx = i + j*nx + k*nx*ny;
                double normal_norm = normals[3*idx+0]*normals[3*idx+0] + normals[3*idx+1]*normals[3*idx+1] + normals[3*idx+2]*normals[3*idx+2];

                if(vertices_unsigned_char[idx] > 0 && normal_norm>0) {
                    vertices_unsigned_char[idx] = 2;
                }
            }
        }
    }
}
