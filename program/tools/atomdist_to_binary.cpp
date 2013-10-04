/*! 
    \file atomdist_to_binary.cpp  <!-- This is needed for doxygen to include documentation of a file without classes, structs etc. -->

    \brief Loads binary files created by distance_to_atom, and puts all the data into one large binary matrix.

    #### Strategy
        We loop over all the nodes used to create the binary files, calculate the linear node index using our library,
        and load the corresponding linear matrix from each file. For each matrix we calculate the origin of the node 
        in the global voxel matrix. From here we loop over the (n_voxel_x, n_voxel_y, n_voxel_z)-cube inside the global 
        matrix and copy the data from the loaded matrix to the global matrix, while keeping in mind that both matrices 
        are linear.
        
    #### Output
        A geometry_m.bin file that can be used as input file in create_world to calculate normal and tangent vectors.

    \param[in] nodes_x
    \param[in] nodes_y
    \param[in] nodes_z
            The number of nodes used to create the binary input files.
    \param[in] threshold
            The distance that defines the voxel boolean value.
    \param[in] binary_data_directory
            The directory of the binary files you want to analyze.
    \param[in] binary_data_filename
            The filename of the binary files you want to analyze (binary_data_directory/binary_data_filename00000)
    \param[out] output_file
            Output file name
*/

#include <cstdlib> // exit(1), atof() etc.
#include <iostream>
#include <iomanip> // setw(), setfill()
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "lib.h" // convert linear index to 3d index and vice versa

using namespace std;

vector<float> load_binary_data(string binary_data_directory, string binary_data_filename, int nodes_x, int nodes_y, int nodes_z, int &n_voxel_x, int &n_voxel_y, int &n_voxel_z);
void print_binary(vector<float> &linear_matrix, string output_file, unsigned int nx, unsigned int ny, unsigned int nz, double threshold, bool binary_output=0);

int main(int args, char *argv[]) {
    if (args != 8) {
        cout << "Run program with " << endl;
        cout << "./atomdist_to_binary  nx  ny  nz  threshold  binary_data_directory  binary_data_filename_base  output_file" << endl;
        exit(1);
    }    
    int nodes_x = atoi(argv[1]);
    int nodes_y = atoi(argv[2]);
    int nodes_z = atoi(argv[3]);
    double threshold = atof(argv[4]);
    string binary_data_directory = argv[5];
    string binary_data_filename = argv[6];
    string output_file = argv[7];
    
    char filename[1000];
    sprintf(filename, "%s", output_file.c_str());
    ofstream file (filename);
    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << " for writing, aborting!" << endl;
        exit(1);
    }
    file.close();

    int n_voxel_x, n_voxel_y, n_voxel_z;
    vector<float> distance_to_nearest_atom_linear_matrix = load_binary_data(binary_data_directory, binary_data_filename, nodes_x, nodes_y, nodes_z, n_voxel_x, n_voxel_y, n_voxel_z);

    int n_points_x = nodes_x*n_voxel_x;
    int n_points_y = nodes_y*n_voxel_y;
    int n_points_z = nodes_z*n_voxel_z;
    print_binary(distance_to_nearest_atom_linear_matrix, output_file, n_points_x, n_points_y, n_points_z, threshold, false);
}

/*!
    Loads linear matrices from the binary files of each node. For each node we calculate the position of the node
    in the system. From this position we loop over a (n_voxel_x, n_voxel_y, n_voxel_z) - cube and copy the data
    from the local linear matrix to the global linear matrix.
*/
vector<float> load_binary_data(string binary_data_directory, string binary_data_filename, int nodes_x, int nodes_y, int nodes_z, int &n_voxel_x, int &n_voxel_y, int &n_voxel_z) {

    ostringstream filename;

    // loading file 0 to find max_distance_angstrom, and the number of voxels (this is the same in all files)
    filename << binary_data_directory << "/" << binary_data_filename << setfill('0') << setw(6) << 0;
    ifstream file(filename.str().c_str(), ios::binary | ios::in);

    float max_distance_angstrom;
    file.read(reinterpret_cast<char*>(&max_distance_angstrom), sizeof(float));
    file.read(reinterpret_cast<char*>(&n_voxel_x), sizeof(int));
    file.read(reinterpret_cast<char*>(&n_voxel_y), sizeof(int));
    file.read(reinterpret_cast<char*>(&n_voxel_z), sizeof(int));
    file.close();

    cout << "Read the following data from file 0:" << endl;
    cout << "max_distance_angstrom = " << max_distance_angstrom << endl;
    cout << "n_voxel_x = " << n_voxel_x << endl;
    cout << "n_voxel_y = " << n_voxel_y << endl;
    cout << "n_voxel_z = " << n_voxel_z << endl;

    int local_linear_matrix_length = n_voxel_x*n_voxel_y*n_voxel_z;
    int total_num_voxels = local_linear_matrix_length*nodes_x*nodes_y*nodes_z;
    int numprocs = nodes_x*nodes_y*nodes_z;
    int global_linear_matrix_length = local_linear_matrix_length*numprocs;;
    int global_matrix_dimensions[3] = {nodes_x*n_voxel_x, nodes_y*n_voxel_y, nodes_z*n_voxel_z};

    vector<float> local_linear_matrix(local_linear_matrix_length);
    vector<float> global_linear_matrix(global_linear_matrix_length);

    // loop over nodes
	for (int i = 0; i < nodes_x; i++) {
		for (int j = 0; j < nodes_y; j++) {
			for (int k = 0; k < nodes_z; k++) {
				int node_id = convert_3d_indices_to_linear_index(i, j, k, nodes_x, nodes_y, nodes_z);

				filename.str(""); // reset stringstream
				filename << binary_data_directory << "/" << binary_data_filename << setfill('0') << setw(6) << node_id;

                cout << "Loading binary file " << filename.str() << endl;

				ifstream file(filename.str().c_str(), ios::binary | ios::in);

                if (!file.is_open()) {
                    cout << "Error: could not open binary file " << filename.str() << ", aborting!" << endl;
                    exit(1);
                }

                float garbage;
                int garbageint;
                file.read(reinterpret_cast<char*>(&garbage), sizeof(float));
                file.read(reinterpret_cast<char*>(&garbageint), sizeof(int));
                file.read(reinterpret_cast<char*>(&garbageint), sizeof(int));
                file.read(reinterpret_cast<char*>(&garbageint), sizeof(int));
                file.read(reinterpret_cast<char*>(&local_linear_matrix[0]), local_linear_matrix_length*sizeof(float));
                file.close();

                // putting the local matrix in the correct place in the global linear matrix
                int local_matrix_origin_in_global_matrix[3];
                convert_linear_index_to_3d_indices(node_id, nodes_x, nodes_y, nodes_z, &local_matrix_origin_in_global_matrix[0]);
                local_matrix_origin_in_global_matrix[0] *= n_voxel_x;
                local_matrix_origin_in_global_matrix[1] *= n_voxel_y;
                local_matrix_origin_in_global_matrix[2] *= n_voxel_z;

                int index_in_global_linear_matrix;
                int index_in_local_linear_matrix;
                int position_in_global_matrix[3] = {
                    local_matrix_origin_in_global_matrix[0],
                    local_matrix_origin_in_global_matrix[1],
                    local_matrix_origin_in_global_matrix[2]
                };

                // looping over the loaded local linear matrix
                for (int i = 0; i < n_voxel_x; i++) {
                    position_in_global_matrix[0] = local_matrix_origin_in_global_matrix[0] + i;
                    for (int j = 0; j < n_voxel_y; j++) {
                        position_in_global_matrix[1] = local_matrix_origin_in_global_matrix[1] + j;

                        index_in_local_linear_matrix = convert_3d_indices_to_linear_index(i, j, 0, n_voxel_x, n_voxel_y, n_voxel_z);
                        index_in_global_linear_matrix = convert_3d_indices_to_linear_index(position_in_global_matrix, global_matrix_dimensions);
                        for (int k = 0; k < n_voxel_z; k++) {
                            global_linear_matrix[index_in_global_linear_matrix + k] = local_linear_matrix[index_in_local_linear_matrix + k];
                        }
                    }
                }
			}
		}
	}
    return global_linear_matrix;
}

/*! 
    Prints a linear matrix to binary format to be used as input in create_world.
*/
void print_binary(vector<float> &linear_matrix, string output_file, unsigned int nx, unsigned int ny, unsigned int nz, double threshold, bool binary_output) {

    stringstream filename;

    int N = nx*ny*nz;

    cout << "Printing to binary file" << endl;

    filename << output_file;
    ofstream *file = new ofstream();
    file->open(filename.str().c_str(), ios::out | ios::binary);
    file->write(reinterpret_cast<char*>(&nx), sizeof(unsigned int));
    file->write(reinterpret_cast<char*>(&ny), sizeof(unsigned int));
    file->write(reinterpret_cast<char*>(&nz), sizeof(unsigned int));

    // column-major ordering...
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                int index = convert_3d_indices_to_linear_index(i, j, k, nx, ny, nz);
                unsigned char voxel_value = linear_matrix[index] < threshold;
                file->write(reinterpret_cast<char*>(&voxel_value), sizeof(unsigned char));
            }
        }
    }

    file->close();
}

