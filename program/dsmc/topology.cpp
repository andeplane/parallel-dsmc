#include "topology.h"
#include <system.h>

vector<int> Topology::index_vector_from_index(const int &index) {
    vector<int> index_vector(3,0);

    index_vector[0] = index/(Topology::num_processors_vector[1]*Topology::num_processors_vector[2]);   // Index in x-direction
    index_vector[1] = (index/Topology::num_processors_vector[2])%Topology::num_processors_vector[1];   // Index in y-direction
    index_vector[2] = index%Topology::num_processors_vector[2];                              // Index in z-direction

    return index_vector;
}

Topology::Topology(int myid_, int nx, int ny, int nz, System *system)
{
    myid = myid_;
    num_processors = nx*ny*nz;
    num_processors_vector.resize(3,0);
    length.resize(3,0);
    origin.resize(3,0);
    num_processors_vector[0] = nx; num_processors_vector[1] = ny; num_processors_vector[2] = nz;
    length[0] = system->length[0]/num_processors_vector[0];
    length[1] = system->length[1]/num_processors_vector[1];
    length[2] = system->length[2]/num_processors_vector[2];
    index_vector = Topology::index_vector_from_index(myid);
    origin[0] = length[0]*index_vector[0];
    origin[1] = length[1]*index_vector[1];
    origin[2] = length[2]*index_vector[2];
}

bool Topology::is_position_inside(double *r) {

    return (r[0] >= origin[0] && r[0] < origin[0]+length[0] &&
            r[1] >= origin[1] && r[1] < origin[1]+length[1] &&
            r[2] >= origin[2] && r[2] < origin[2]+length[2]);
}
