#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>

using namespace std;

double L0 = 1e4;

int main(int args, char *argv[]) {
	if(args < 4) {
		cout << "Please specify the number of cpus, state folder and output filename." << endl;
		return 0;
	}
	int cpus = atoi(argv[1]);
	string state_folder = argv[2];
	string output_filename = argv[3];

	double *molecule_data = new double[9*1000000];
	ofstream file (output_filename.c_str(), ios::out);

	ifstream **state_files = new ifstream*[cpus];
	for(int cpu=0;cpu<cpus;cpu++) {
		char *filename = new char[100];
		sprintf(filename,"%s/state%04d.bin",state_folder.c_str(), cpu);
		state_files[cpu] = new ifstream(filename,ios::in | ios::binary);
		if (state_files[cpu]->is_open()){
	        std::cout << "Reading file " << filename << "\n";
	    }
	    else {    
	        std::cout << "Failed to open file " << filename << ". Aborting." << std::endl;
	        exit(1);
	    }
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
    	file << "H " << L0*molecule_data[9*n+0] << " " << L0*molecule_data[9*n+1] << " " << L0*molecule_data[9*n+2] << endl;
    }

    cout << "Created XYZ-file with " << num_molecules << " molecules." << endl;

	file.close();
	for(int cpu=0;cpu<cpus;cpu++) {
		state_files[cpu]->close();
	}

	return 0;
}