#include <topology.h>
#include <system.h>

Topology::Topology(int myid_, int nx, int ny, int nz, System *system_)
{
    // Fill in node_id to facet_id map
    system = system_;
    myid = myid_;
    num_processors = nx*ny*nz;
    num_processors_vector.resize(3,0);
    length.resize(3,0);
    origin.resize(3,0);
    num_processors_vector[0] = nx; num_processors_vector[1] = ny; num_processors_vector[2] = nz;
    length[0] = system->length[0]/num_processors_vector[0];
    length[1] = system->length[1]/num_processors_vector[1];
    length[2] = system->length[2]/num_processors_vector[2];
    index_vector = index_vector_from_index(myid);
    origin[0] = length[0]*index_vector[0];
    origin[1] = length[1]*index_vector[1];
    origin[2] = length[2]*index_vector[2];
    node_id_to_facet_id_list.resize(num_processors,-1);
    for(int i=0; i<num_processors_vector[0]; i++) {
        for(int j=0; j<num_processors_vector[1]; j++) {
            for(int k=0; k<num_processors_vector[2]; k++) {
                int index = index_from_ijk(i,j,k);
                int dx = i - index_vector[0];
                int dy = j - index_vector[1];
                int dz = k - index_vector[2];

                if(dx == -1 || i == num_processors_vector[0]-1 && index_vector[0] == 0) {
                    // The node is to our left in x-dir
                    node_id_to_facet_id_list[index] = 0;
                    continue;
                }
                if(dx == 1 || i == 0 && index_vector[0] == num_processors_vector[0]-1) {
                    // The node is to our right in x-dir
                    node_id_to_facet_id_list[index] = 1;
                    continue;
                }
                if(dy == -1 || j == num_processors_vector[1]-1 && index_vector[1] == 0) {
                    // The node is to our left in y-dir
                    node_id_to_facet_id_list[index] = 2;
                    continue;
                }
                if(dy == 1 || j == 0 && index_vector[1] == num_processors_vector[1]-1) {
                    // The node is to our right in y-dir
                    node_id_to_facet_id_list[index] = 3;
                    continue;
                }
                if(dz == -1 || k == num_processors_vector[2]-1 && index_vector[2] == 0) {
                    // The node is to our left in z-dir
                    node_id_to_facet_id_list[index] = 4;
                    continue;
                }
                if(dz == 1 || k == 0 && index_vector[2] == num_processors_vector[2]-1) {
                    // The node is to our right in z-dir
                    node_id_to_facet_id_list[index] = 5;
                    continue;
                }
            }
        }
    }
    facet_id_to_node_id_list.resize(6,0);
    int a = (index_vector[0] - 1 + num_processors_vector[0]) % num_processors_vector[0]; // i - 1
    facet_id_to_node_id_list[0] = index_from_ijk(a, index_vector[1], index_vector[2]);
    a = (index_vector[0] + 1 + num_processors_vector[0]) % num_processors_vector[0];  // i + 1
    facet_id_to_node_id_list[1] = index_from_ijk(a, index_vector[1], index_vector[2]);

    a = (index_vector[1] - 1 + num_processors_vector[1]) % num_processors_vector[1];  // j - 1
    facet_id_to_node_id_list[2] = index_from_ijk(index_vector[0], a, index_vector[2]);
    a = (index_vector[1] + 1 + num_processors_vector[1]) % num_processors_vector[1];  // j + 1
    facet_id_to_node_id_list[3] = index_from_ijk(index_vector[0], a, index_vector[2]);
    a = (index_vector[2] - 1 + num_processors_vector[2]) % num_processors_vector[2];  // k - 1
    facet_id_to_node_id_list[4] = index_from_ijk(index_vector[0], index_vector[1], a);
    a = (index_vector[2] + 1 + num_processors_vector[2]) % num_processors_vector[2];  // k + 1
    facet_id_to_node_id_list[5] = index_from_ijk(index_vector[0], index_vector[1], a);
}

vector<int> Topology::index_vector_from_index(const int &index) {
    vector<int> index_vector(3,0);

    index_vector[0] = index/(num_processors_vector[1]*num_processors_vector[2]);   // Index in x-direction
    index_vector[1] = (index/num_processors_vector[2])%num_processors_vector[1];   // Index in y-direction
    index_vector[2] = index%num_processors_vector[2];                              // Index in z-direction

    return index_vector;
}

int Topology::index_from_ijk(int i, int j, int k) {
    return i*num_processors_vector[1]*num_processors_vector[2] + j*num_processors_vector[2] + k;
}

int Topology::index_from_position(const double *r) {
    int i = r[0] / system->length[0] * num_processors_vector[0];
    int j = r[1] / system->length[1] * num_processors_vector[1];
    int k = r[2] / system->length[2] * num_processors_vector[2];

    return index_from_ijk(i,j,k);
}

vector<int> Topology::index_vector_from_position(const double *r) {
    vector<int> index_vector(3,0);

    index_vector[0] = r[0] / system->length[0] * num_processors_vector[0];
    index_vector[1] = r[1] / system->length[1] * num_processors_vector[1];
    index_vector[2] = r[2] / system->length[2] * num_processors_vector[2];

    return index_vector;
}

bool Topology::is_position_inside(double *r) {
    return (r[0] >= origin[0] && r[0] < origin[0]+length[0] &&
            r[1] >= origin[1] && r[1] < origin[1]+length[1] &&
            r[2] >= origin[2] && r[2] < origin[2]+length[2]);
}
