#include <topology.h>
#include <system.h>

Topology::Topology(int myid_, int nx, int ny, int nz, System *system_)
{
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
}

vector<int> Topology::index_vector_from_index(const int &index) {
    vector<int> index_vector(3,0);

    index_vector[0] = index/(num_processors_vector[1]*num_processors_vector[2]);   // Index in x-direction
    index_vector[1] = (index/num_processors_vector[2])%num_processors_vector[1];   // Index in y-direction
    index_vector[2] = index%num_processors_vector[2];                              // Index in z-direction

    return index_vector;
}

int Topology::index_from_position(const double *r) {
    int i = r[0] / system->length[0] * num_processors_vector[0];
    int j = r[1] / system->length[1] * num_processors_vector[1];
    int k = r[2] / system->length[2] * num_processors_vector[2];

    return i*num_processors_vector[1]*num_processors_vector[2] + j*num_processors_vector[2] + k;
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
