#ifndef LIB_H
#define LIB_H

int convert_3d_indices_to_linear_index(              int* index_vector, int *nx_ny_nz);
int convert_3d_indices_to_linear_index(              int* index_vector, int nx, int ny, int nz);
int convert_3d_indices_to_linear_index(              int index_x, int index_y, int index_z, int *nx_ny_nz);
int convert_3d_indices_to_linear_index(              int index_x, int index_y, int index_z, int nx, int ny, int nz);
inline int convert_3d_indices_to_linear_index_inline(int index_x, int index_y, int index_z, int nx, int ny, int nz);

void convert_linear_index_to_3d_indices(              int linear_index, int *nx_ny_nz,          int *index_vector);
void convert_linear_index_to_3d_indices(              int linear_index, int nx, int ny, int nz, int *index_vector);
void convert_linear_index_to_3d_indices(              int linear_index, int *nx_ny_nz,          int &index_x, int &index_y, int &index_z);
inline void convert_linear_index_to_3d_indices_inline(int linear_index, int nx, int ny, int nz, int &index_x, int &index_y, int &index_z);

int convert_3d_indices_to_linear_index(int* index_vector, int *nx_ny_nz) {
    return convert_3d_indices_to_linear_index_inline(index_vector[0], index_vector[1], index_vector[2], nx_ny_nz[0], nx_ny_nz[1], nx_ny_nz[2]);
}

int convert_3d_indices_to_linear_index(int* index_vector, int nx, int ny, int nz) {
    return convert_3d_indices_to_linear_index_inline(index_vector[0], index_vector[1], index_vector[2], nx, ny, nz);
}

int convert_3d_indices_to_linear_index(int index_x, int index_y, int index_z, int *nx_ny_nz) {
    return convert_3d_indices_to_linear_index_inline(index_x, index_y, index_z, nx_ny_nz[0], nx_ny_nz[1], nx_ny_nz[2]);
}

int convert_3d_indices_to_linear_index(int index_x, int index_y, int index_z, int nx, int ny, int nz) {
    return convert_3d_indices_to_linear_index_inline(index_x, index_y, index_z, nx, ny, nz);
}

void convert_linear_index_to_3d_indices(int linear_index, int *nx_ny_nz, int *index_vector) {
    convert_linear_index_to_3d_indices_inline(linear_index, nx_ny_nz[0], nx_ny_nz[1], nx_ny_nz[2], index_vector[0], index_vector[1], index_vector[2]);
}

void convert_linear_index_to_3d_indices(int linear_index, int nx, int ny, int nz, int *index_vector) {
    convert_linear_index_to_3d_indices_inline(linear_index, nx, ny, nz, index_vector[0], index_vector[1], index_vector[2]);
}

void convert_linear_index_to_3d_indices(int linear_index, int *nx_ny_nz, int &index_x, int &index_y, int &index_z) {
    convert_linear_index_to_3d_indices_inline(linear_index, nx_ny_nz[0], nx_ny_nz[1], nx_ny_nz[2], index_x, index_y, index_z);
}

void convert_linear_index_to_3d_indices(int linear_index, int nx, int ny, int nz, int &index_x, int &index_y, int &index_z) {
    convert_linear_index_to_3d_indices_inline(linear_index, nx, ny, nz, index_x, index_y, index_z);
}


/*
 Below are "private" inline functions that are used by the functions above.
 This is so we get the functions inlined in the functions above, but only have
 to edit one function
*/

inline int convert_3d_indices_to_linear_index_inline(int index_x, int index_y, int index_z, int nx, int ny, int nz) {
    /* Converts 3d indices to linear index */
    return index_x*ny*nz + index_y*nz + index_z;
}

inline void convert_linear_index_to_3d_indices_inline(int linear_index, int nx, int ny, int nz, int &index_x, int &index_y, int &index_z) {
    /* Converts linear index to 3d indices */
    index_x = linear_index/(ny*nz);   // Index in x-direction
    index_y = (linear_index/nz)%ny;   // Index in y-direction
    index_z = linear_index%nz;        // Index in z-direction
}

#endif // LIB_H
