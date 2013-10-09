// array::fill example
#include <iostream>
#include <fstream>
#include <math.h>
#include <perlin.cpp>
#include <vector>

using namespace std;

int main (int args, char *argv[]) {
	if(args < 3) {
		cout << "Run program with:" << endl << "./create_world_perlin nx ny nz" << endl;
		return 0;
	}

    unsigned int nx = atoi(argv[1]);
    unsigned int ny = atoi(argv[2]);
    unsigned int nz = atoi(argv[3]);
    long points = nx*ny*nz;

    vector<unsigned char> M;
    M.resize(points,0);
	
    // oktav, frekvens, amplitude , seed
    Perlin p(1, 1, 1, 3);

    // Copy all values from the original matrix into the extended matrix
    for(int i=0;i<nx;i++) {
        for(int j=0;j<ny;j++) {
            for(int k=0;k<nz;k++) {
                double x = (i-nx/2.0)/(double)nx;
                double y = (j-ny/2.0)/(double)ny;
                double z = (k-nz/2.0)/(double)nz;
                float s = 1.0;

                int index = i + j*nx + k*nx*ny;
                double val = 0;
                for (int a=0; a<5  ; a++) {
                    s = 3.0*a + 2.2513531;
                    
                    val += p.Get(x*s, y*s, z*s);
                }
                // val = pow(val,4.0)*cos(val);;
                if(val < 0.2) M[index] = 0;
                else M[index] = 1;
                // M[index] = 1;
            }
        }
    }

    ofstream file ("perlin_m.bin", ios::out | ios::binary);
    file.write (reinterpret_cast<char*>(&nx), sizeof(unsigned int));
    file.write (reinterpret_cast<char*>(&ny), sizeof(unsigned int));
    file.write (reinterpret_cast<char*>(&nz), sizeof(unsigned int));
    file.write (reinterpret_cast<char*>(&M[0]), points*sizeof(unsigned char));
    file.close();
	
	return 0;
}