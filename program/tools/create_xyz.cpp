#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

double L0 = 1e4;

int main(int args, char *argv[]) {
	if(args < 2) {
		cout << "Please specify the number of cpus." << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	
	double *molecule_data = new double[9*1000000];
	ofstream file ("state.xyz", ios::out);
	
	ifstream **state_files = new ifstream*[cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"state_files/state%04d.bin",cpu);
		state_files[cpu] = new ifstream(filename,ios::in | ios::binary);
	}
	cout << cpus << " state files opened." << endl;
	
	int num_molecules = 0;
	for(int cpu=0;cpu<cpus;cpu++) { 
		int N;
		state_files[cpu]->read(reinterpret_cast<char*>(&N),sizeof(int));
		state_files[cpu]->read(reinterpret_cast<char*>(&molecule_data[9*num_molecules]),9*N*sizeof(double));
		num_molecules += N;
	}
	

	file << num_molecules << endl;
	file << "sup" << endl;
	for(int n=0;n<num_molecules;n++) {
    	// We return height - r(1) because system is inverted
    	file << "Si " << L0*molecule_data[9*n+0] << " " << L0*molecule_data[9*n+1] << " " << L0*molecule_data[9*n+2] << " " << n << endl;
    }

    cout << "Created XYZ-file with " << num_molecules << " molecules." << endl;

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		state_files[cpu]->close();
	}

	return 0;
}