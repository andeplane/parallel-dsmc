#include <dsmc_io.h>

#include <system.h>
#include <cell.h>
#include <mpi.h>
#include <grid.h>
#include <settings.h>
#include <dsmctimer.h>
#include <topology.h>

DSMC_IO::DSMC_IO(System *system_) {
    system = system_;
    settings = system->settings;
    movie_frames = 0;
    movie_file_open = false;

    char *time_statistics_filename = new char[100];
    sprintf(time_statistics_filename,"statistics/time/time%04d.txt",system->myid);
    time_statistics_file = new ofstream(time_statistics_filename,ios::out | ios::binary);
    delete time_statistics_filename;

    if(system->myid==0) {
        energy_file = fopen("statistics/energy.txt","w");
        velocity_file = fopen("statistics/velocity.txt","w");
        flux_file = fopen("statistics/flux.txt","w");
        permeability_file = fopen("statistics/permeability.txt","w");
        density_file = fopen("statistics/density.txt","w");
        linear_density_file = fopen("statistics/linear_density.txt","w");
        num_molecules_file = fopen("statistics/num_molecules_file.txt","w");
        pressure_file = fopen("statistics/num_molecules_file.txt","w");
        linear_pressure_file = fopen("statistics/linear_pressure.txt","w");
    }
}

void DSMC_IO::save_state_to_movie_file() {
    system->timer->start_io();
    if(settings->create_movie && !(system->steps % settings->movie_every_n_frame)) {
        if(!movie_file_open) {
            char *filename = new char[100];
            sprintf(filename,"movie_files/movie%04d.bin",system->myid);
            movie_file = new ofstream(filename,ios::out | ios::binary);
            movie_file_open = true;
            data = new double[3*MAX_MOLECULE_NUM];
            memset(data,0,3*MAX_MOLECULE_NUM*sizeof(float));
            delete filename;
        }

        int count = 0;
        for(unsigned int n=0;n<system->num_molecules_local;n++) {
            data[count++] = system->r[3*n+0];
            data[count++] = system->r[3*n+1];
            data[count++] = system->r[3*n+2];
        }

        count /= 4; // This should represent the number of particles, 4 doubles per particle
        int actual_movie_molecules = settings->movie_molecules / system->topology->num_processors;

        movie_file->write (reinterpret_cast<char*>(&actual_movie_molecules), sizeof(int));
        movie_file->write (reinterpret_cast<char*>(data), 3*actual_movie_molecules*sizeof(double));
    }
    system->timer->end_io();
}

void DSMC_IO::save_state_to_file_binary() {
    system->timer->start_io();
    if(system->myid==0) cout << "Saving state to file..." << endl;

    int N = system->num_molecules_local;

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);

    ofstream file (filename, ios::out | ios::binary);

    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }

    double *tmp_data = new double[9*N];

    int count = 0;

    for(unsigned int n=0;n<system->num_molecules_local;n++) {
        tmp_data[count++] = system->r[3*n+0];
        tmp_data[count++] = system->r[3*n+1];
        tmp_data[count++] = system->r[3*n+2];

        tmp_data[count++] = system->v[3*n+0];
        tmp_data[count++] = system->v[3*n+1];
        tmp_data[count++] = system->v[3*n+2];

        tmp_data[count++] = system->r0[3*n+0];
        tmp_data[count++] = system->r0[3*n+1];
        tmp_data[count++] = system->r0[3*n+2];
    }

    file.write (reinterpret_cast<char*>(&N), sizeof(int));
    file.write (reinterpret_cast<char*>(tmp_data), 9*N*sizeof(double));

    file.close();
    delete tmp_data;
    delete filename;
    system->timer->end_io();
}

void DSMC_IO::load_state_from_file_binary() {
    system->timer->start_io();
    if(system->myid==0) cout << "Loading state from file..." << endl;
    int N = 0;

    char *filename = new char[100];
    sprintf(filename,"state_files/state%04d.bin",system->myid);
    ifstream file (filename, ios::in | ios::binary);
    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }

    file.read(reinterpret_cast<char*>(&N),sizeof(int));
    double *tmp_data = new double[9*N];

    file.read(reinterpret_cast<char*>(tmp_data), 9*N*sizeof(double));
    file.close();

    for(int n=0;n<N;n++) {
        double *r = &system->r[3*n];
        double *v = &system->v[3*n];
        double *r0 = &system->r0[3*n];

        r[0] = tmp_data[9*n+0];
        r[1] = tmp_data[9*n+1];
        r[2] = tmp_data[9*n+2];
        v[0] = tmp_data[9*n+3];
        v[1] = tmp_data[9*n+4];
        v[2] = tmp_data[9*n+5];
        r0[0] = tmp_data[9*n+6];
        r0[1] = tmp_data[9*n+7];
        r0[2] = tmp_data[9*n+8];

        Cell *cell = system->all_cells[system->cell_index_from_position(r)];
        cell->add_molecule(n,system->molecule_index_in_cell,system->molecule_cell_index);
    }

    system->num_molecules_local = N;

    delete filename;
    delete tmp_data;
    system->timer->end_io();
}

void DSMC_IO::finalize() {
    if(movie_file_open) {
        movie_file->close();
    }
    time_statistics_file->close();

    if(system->myid != 0) return;
    fclose(energy_file);
    fclose(velocity_file);
    fclose(flux_file);
    fclose(permeability_file);
    fclose(num_molecules_file);
    fclose(linear_density_file);
    fclose(density_file);
    fclose(pressure_file);
    fclose(linear_pressure_file);
}

void DSMC_IO::read_grid_matrix(string filename, Grid *grid) {
    ifstream file (filename.c_str(), ios::in | ios::binary);
    if(!file.is_open()) {
        cout << "Error, could not open file " << filename << endl;
        exit(1);
    }
    int Nx, Ny, Nz, points;

    file.read (reinterpret_cast<char*>(&Nx), sizeof(int));
    file.read (reinterpret_cast<char*>(&Ny), sizeof(int));
    file.read (reinterpret_cast<char*>(&Nz), sizeof(int));
    points = Nx*Ny*Nz;
    grid->Nx = Nx; grid->Ny = Ny; grid->Nz = Nz; grid->points = points;

    grid->voxels = new unsigned char[points];
    grid->normals   = new float[3*points];
    grid->tangents1 = new float[3*points];
    grid->tangents2 = new float[3*points];

    file.read (reinterpret_cast<char*>(grid->voxels), points*sizeof(unsigned char));
    file.read (reinterpret_cast<char*>(grid->normals), 3*points*sizeof(float));
    file.read (reinterpret_cast<char*>(grid->tangents1), 3*points*sizeof(float));
    file.read (reinterpret_cast<char*>(grid->tangents2), 3*points*sizeof(float));
    file.close();
}
