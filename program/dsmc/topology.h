#pragma once
#include <vector>
#include <defines.h>
using std::vector;

class System;


class Topology
{
public:
    vector<double> length;
    vector<int>    num_processors_vector;
    vector<int>    node_id_to_facet_id_list;
    vector<int>    facet_id_to_node_id_list;
    vector<int>    my_parity;
    int            num_processors;

    vector<double> origin;
    vector<int>    index_vector;
    int            myid;
    System        *system;

    Topology(int myid_, int nx, int ny, int nz, System *system);
    vector<int> index_vector_from_index(const int &index);
    vector<int> index_vector_from_position(const vector<double> &r, const int &index);
    bool is_position_inside(const vector<double> &r, const int &index);
    int index_from_position(data_type *r);
    int index_from_ijk(int i, int j, int k);
    int index_from_molecule_index(const int &index);

};
