#pragma once
#include <vector>
using std::vector;

class System;


class Topology
{
public:
    vector<double> length;
    vector<int>    num_processors_vector;
    int            num_processors;

    vector<double> origin;
    vector<int>    index_vector;
    int            myid;
    System        *system;

    Topology(int myid_, int nx, int ny, int nz, System *system);
    vector<int> index_vector_from_index(const int &index);
    vector<int> index_vector_from_position(const double *r);
    bool is_position_inside(double *r);
    int index_from_position(const double *r);

};
