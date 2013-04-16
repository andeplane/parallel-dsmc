// array::fill example
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

class Random {
public:
	long     iy;
	long     iv[NTAB];
	long     idum[1];

	Random(long seed) {
		*idum = seed;
   		iy = 0;
	}

	double nextDouble() {
		int             j;
		long            k;
		double          temp;

		if (*idum <= 0 || !iy) {
		  if (-(*idum) < 1) *idum=1;
		  else *idum = -(*idum);
		  for(j = NTAB + 7; j >= 0; j--) {
		     k     = (*idum)/IQ;
		     *idum = IA*(*idum - k*IQ) - IR*k;
		     if(*idum < 0) *idum += IM;
		     if(j < NTAB) iv[j] = *idum;
		  }
		  iy = iv[0];
		}
		k     = (*idum)/IQ;
		*idum = IA*(*idum - k*IQ) - IR*k;
		if(*idum < 0) *idum += IM;
		j     = iy/NDIV;
		iy    = iv[j];
		iv[j] = *idum;
		if((temp=AM*iy) > RNMX) return RNMX;
		else return temp;
	}
};

void calculate_normals(int Nx, int Ny, int Nz, float *normal, unsigned char *M) {
	int idx = 0;
	int idx2 = 0;
    float norm;
    bool at_least_one_wall_neighbor;
    bool all_neighbors_are_walls;
    for(int i=0;i<Nx;i++) {
        for(int j=0;j<Ny;j++) {
            for(int k=0;k<Nz;k++) {
            	at_least_one_wall_neighbor = false;
            	all_neighbors_are_walls = true;

            	idx = i + j*Nx + k*Nx*Ny;
                for(int di=-1;di<=1;di++) {
                    for(int dj=-1;dj<=1;dj++) {
                        for(int dk=-1;dk<=1;dk++) {
                        	idx2 = ((i+di+Nx)%Nx) + ((j+dj+Ny)%Ny)*Nx+ ((k+dk+Nz)%Nz)*Nx*Ny;

                        	if(M[idx2]>0) {
                        		// If at least one wall neighbor, this is not a single wall voxel
                        		at_least_one_wall_neighbor = true;
                        	}
                        	if(M[idx2]==0) {
                        		// If not all neighbors are walls and we have norm=1, this is a
                        		// single plane that has no defined normal vector.
                        		all_neighbors_are_walls = false;
                        	}

                        	normal[3*idx+0] -= M[idx2]*di;
                        	normal[3*idx+1] -= M[idx2]*dj;
                        	normal[3*idx+2] -= M[idx2]*dk;
                        }
                    }
                }

                norm = sqrt(normal[3*idx+0]*normal[3*idx+0] + normal[3*idx+1]*normal[3*idx+1] + normal[3*idx+2]*normal[3*idx+2]);
                
                if(norm > 0) {
                	normal[3*idx+0] /= norm;
                	normal[3*idx+1] /= norm;
                	normal[3*idx+2] /= norm;
                } else if(!at_least_one_wall_neighbor || !all_neighbors_are_walls) {
                	// Single point or single pixel-plane, should not be a wall
                	M[idx] = 0;
                }

                // if(norm>0 && k==10) {
                //     cout << "n(" << i << "," << j << ",0)=(" << normal[3*idx+0] << "," << normal[3*idx+1] << "," << normal[3*idx+2] << ")." << "     M=" << int(M[idx]) << endl;
                // }
            }
        }
    }
}

void calculate_tangents(int Nx, int Ny, int Nz, float *tangent1, float *tangent2, float *normal, unsigned char *M) {
	int idx = 0;
    float norm, dot_product;
    Random *rnd = new Random(-1);

	for(int i=0;i<Nx;i++) {
        for(int j=0;j<Ny;j++) {
        	for(int k=0;k<Nz;k++) {
            	idx = i + j*Nx + k*Nx*Ny;
                tangent1[3*idx+0] = rnd->nextDouble();
                tangent1[3*idx+1] = rnd->nextDouble();
                tangent1[3*idx+2] = rnd->nextDouble();
                
                dot_product = normal[3*idx+0]*tangent1[3*idx+0] + normal[3*idx+1]*tangent1[3*idx+1] + normal[3*idx+2]*tangent1[3*idx+2];
                
                // Perform gram-schmidt
                tangent1[3*idx+0] -= normal[3*idx+0]*dot_product;
                tangent1[3*idx+1] -= normal[3*idx+1]*dot_product;
                tangent1[3*idx+2] -= normal[3*idx+2]*dot_product;

                // Normalize
                norm = sqrt(tangent1[3*idx+0]*tangent1[3*idx+0] + tangent1[3*idx+1]*tangent1[3*idx+1] + tangent1[3*idx+2]*tangent1[3*idx+2]);

                if(norm>0) {
                	tangent1[3*idx+0] /= norm;
	                tangent1[3*idx+1] /= norm;
	                tangent1[3*idx+2] /= norm;
                }
                
                // t2 = n x t1
                tangent2[3*idx+0] = tangent1[3*idx+1]*normal[3*idx+2] - tangent1[3*idx+2]*normal[3*idx+1];
                tangent2[3*idx+1] = tangent1[3*idx+2]*normal[3*idx+0] - tangent1[3*idx+0]*normal[3*idx+2];
                tangent2[3*idx+2] = tangent1[3*idx+0]*normal[3*idx+1] - tangent1[3*idx+1]*normal[3*idx+0];

                // Normalize
                norm = sqrt(tangent2[3*idx+0]*tangent2[3*idx+0] + tangent2[3*idx+1]*tangent2[3*idx+1] + tangent2[3*idx+2]*tangent2[3*idx+2]);

                if(norm>0) {
                	tangent2[3*idx+0] /= norm;
	                tangent2[3*idx+1] /= norm;
	                tangent2[3*idx+2] /= norm;
                }
            }
        }
    }
}

void calculate_inner_points(int Nx, int Ny, int Nz, float *normal, unsigned char *M) {
	int idx = 0;
    float norm;
    for(int i=0;i<Nx;i++) {
    	for(int j=0;j<Ny;j++) {
			for(int k=0;k<Nz;k++) {
                idx = i + j*Nx + k*Nx*Ny;
                norm = normal[3*idx+0]*normal[3*idx+0] + normal[3*idx+1]*normal[3*idx+1] + normal[3*idx+2]*normal[3*idx+2];

                if(M[idx] > 0 && norm>0) {
                	M[idx] = 2;
                }

                // printf("M(%d,%d,%d) = %d\n",i,j,k,M[idx]);
                // printf("N(%d,%d,%d) = (%f,%f,%f)\n",i,j,k,normal[3*idx+0],normal[3*idx+1],normal[3*idx+2]);
                // cout << norm << endl << endl;
            }
        }
    }
}

void save_to_file(char *outfile, int Nx, int Ny, int Nz, float *normal, float *tangent1, float *tangent2, unsigned char *M) {
	int points = Nx*Ny*Nz;
	ofstream file (outfile, ios::out | ios::binary);
    file.write (reinterpret_cast<char*>(&Nx), sizeof(int));
    file.write (reinterpret_cast<char*>(&Ny), sizeof(int));
    file.write (reinterpret_cast<char*>(&Nz), sizeof(int));
    file.write (reinterpret_cast<char*>(M), points*sizeof(unsigned char));
    file.write (reinterpret_cast<char*>(normal),   3*points*sizeof(float));
    file.write (reinterpret_cast<char*>(tangent1), 3*points*sizeof(float));
    file.write (reinterpret_cast<char*>(tangent2), 3*points*sizeof(float));
    file.close();
}

int main (int args, char *argv[]) {
	if(args < 3) {
		cout << "Run program with:" << endl << "./create_world infile outfile reservoir_dimension reservoir_size (fraction of total size)" << endl;
		return 0;
	}

    unsigned char N0[3];
    unsigned char N[3];

	char *infile = argv[1];
	char *outfile = argv[2];
    int reservoir_dimension = atoi(argv[3]);
    double reservoir_size_fraction = atof(argv[4]);
    int reservoir_size = 0;

	ifstream file (infile, ios::in | ios::binary);
	file.read (reinterpret_cast<char*>(&N0[0]), sizeof(unsigned char));
	file.read (reinterpret_cast<char*>(&N0[1]), sizeof(unsigned char));
	file.read (reinterpret_cast<char*>(&N0[2]), sizeof(unsigned char));
	cout << "Matrix is " << int(N0[0]) << "x" << int(N0[1]) << "x" << int(N0[2]) << endl;
	int points = N0[0]*N0[1]*N0[2];

    N[0] = N0[0];
    N[1] = N0[1];
    N[2] = N0[2];

	unsigned char *M = new unsigned char[points];
    file.read (reinterpret_cast<char*>(M), points*sizeof(unsigned char));
    file.close();

    // Create a slightly larger matrix so that we have space for the reservoirs
    reservoir_size = N[reservoir_dimension]*reservoir_size_fraction/(2-2*reservoir_size_fraction);
    N[reservoir_dimension] += 2*reservoir_size;
    points = N[0]*N[1]*N[2];
    unsigned char *M2 = new unsigned char[points];

    memset(M2,0,points*sizeof(unsigned char)); // set all values to zero

    // Copy all values from the original matrix into the extended matrix
    for(int i=0;i<N0[0];i++) {
        for(int j=0;j<N0[1];j++) {
            for(int k=0;k<N0[2];k++) {
                int index_M = i + j*N0[0] + k*N0[0]*N0[1];
                int i_M2 = i+reservoir_size*(reservoir_dimension==0);
                int j_M2 = j+reservoir_size*(reservoir_dimension==1);
                int k_M2 = k+reservoir_size*(reservoir_dimension==2);
                int index_M2 = i_M2 + j_M2*N[0] + k_M2*N[0]*N[1];
                M2[index_M2] = M[index_M];
            }
        }
    }

    int idx[3];
    if(reservoir_dimension>=0) {
        for(idx[0]=0;idx[0]<N[0];idx[0]++) {
            for(idx[1]=0;idx[1]<N[1];idx[1]++) {
                for(idx[2]=0;idx[2]<N[2];idx[2]++) {

                    int index_M2 = idx[0] + idx[1]*N[0] + idx[2]*N[0]*N[1];
                    if(idx[reservoir_dimension] == 0 || idx[reservoir_dimension] == N[reservoir_dimension]-1) {
                        // Create hard walls at the ends
                        M2[index_M2] = 1;
                    }
                }
            }
        }
    }

    float *normal = new float[3*points];
    float *tangent1 = new float[3*points];
    float *tangent2 = new float[3*points];
    // memset(normal,0,3*points*sizeof(float));
    // memset(tangent1,0,3*points*sizeof(float));
    // memset(tangent2,0,3*points*sizeof(float));

    for(int i=0;i<3*points;i++) {
    	normal[i] = 0;
    	tangent1[i] = 0;
		tangent2[i] = 0;
    }

    cout << "Creating normals..." << endl;
    calculate_normals(N[0],N[1],N[2],normal,M2);
    cout << "Creating tangents..." << endl;
    calculate_tangents(N[0],N[1],N[2],tangent1,tangent2,normal,M2);
    cout << "Creating boundary..." << endl;
    calculate_inner_points(N[0],N[1],N[2],normal,M2);
    cout << "Creating done!" << endl;
    save_to_file(outfile,N[0],N[1],N[2], normal,tangent1,tangent2,M2);
	
	return 0;
}