#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

double L0 = 1e4;

int main(int args, char *argv[]) {
	if(args < 2) {
		cout << "Please specify project folder, number of CPUs and timesteps." << endl;
		return 0;
	}

	char *base_name = argv[1];
	int cpus = atoi(argv[2]);
	int timesteps = atoi(argv[3]);
	
	float *positions = new double[3*1000000];
	ofstream file ("movie.xyz", ios::out);
	
	ifstream **movie_files = new ifstream*[cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"%s/movie_files/movie%04d.bin",base_name,cpu);
		cout << "File: " << filename << endl;
		movie_files[cpu] = new ifstream(filename,ios::in | ios::binary);
	}
	cout << cpus << " state files opened." << endl;

	for(int timestep=0;timestep<timesteps;timestep++) {
		int num_particles = 0;
		for(int cpu=0;cpu<cpus;cpu++) { 
			int N;
			movie_files[cpu]->read(reinterpret_cast<char*>(&N),sizeof(int));
			movie_files[cpu]->read(reinterpret_cast<char*>(&positions[3*num_particles]),3*N*sizeof(float));
			num_particles += N;
		}
		

		file << num_particles << endl;
		file << "sup" << endl;
		for(int n=0;n<num_particles;n++) {
        	// We return height - r(1) because system is inverted
        	file << "Ar " << L0*positions[3*n+0] << " " << L0*positions[3*n+1] << " " << L0*positions[3*n+2] << endl;
	    }

	    cout << "Wrote timestep " << timestep << endl;
	}

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		movie_files[cpu]->close();
	}
	cout << "Movie created." << endl;

	return 0;
}