#include <moviedata.h>
#include <fstream>
#include <iostream>
#include <cvector.h>
#include <system.h>
#include <solver.h>
#include <GL/glfw.h>      // Include OpenGL Framework library
#include <ctexture.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;

Timestep::Timestep() {
    next = NULL;
    solver = NULL;
    previous = NULL;
}

void Timestep::render_billboards(CTexture *texture, float scale) {
    int num_molecules = positions.size()/3;
    vector<int> steps_since_last_collision(num_molecules,100);
    texture->render_billboards(positions, velocities, steps_since_last_collision, num_molecules, scale);
}

void Timestep::add_molecule_data(vector<float> data) {
    int num_molecules_in_this_timestep = data.size() / 6;

    vector<double> new_positions_double;
    vector<double> new_velocities_double;
    new_positions_double.resize(num_molecules_in_this_timestep*3);
    new_velocities_double.resize(num_molecules_in_this_timestep*3);

    for(int i=0; i<num_molecules_in_this_timestep; i++) {
        for(int a=0; a<3; a++) {
            new_positions_double[3*i+a] = data[6*i+a];
            new_velocities_double[3*i+a] = data[6*i + 3 + a];
        }
    }
    positions.insert(positions.end(), new_positions_double.begin(), new_positions_double.end());
    velocities.insert(velocities.end(), new_velocities_double.begin(), new_velocities_double.end());

    new_positions_double.clear();
    new_velocities_double.clear();
}

void Timestep::render_points() {
    int num_molecules = positions.size() / 3;
    glColor4f(1.0, 1.0, 1.0, 1.0);
    glPointSize(3.0);
    glBegin(GL_POINTS);
    for(int i=0; i<num_molecules; i++) {
        glVertex3f( positions[3*i+0], positions[3*i+1], positions[3*i+2]);
    }
    glEnd();
}

MovieData::MovieData(int cpus_, int timesteps_)
{
    cpus = cpus_;
    timesteps = timesteps_;
}

void MovieData::load_movie_files(string state_folder) {
    cout << "Will load " << timesteps << " timesteps from " << cpus << " cpus at " << state_folder << endl;
    char filename[1000];
    int actual_movie_molecules = 0;
    vector<float> data;
    Timestep *previous_timestep = NULL;
    for(int timestep=0; timestep < timesteps; timestep++) {
        Timestep *timestep_object = new Timestep();
        if(previous_timestep == NULL) first_timestep = timestep_object;
        else {
            timestep_object->previous = previous_timestep;
            previous_timestep->next = timestep_object;
        }

        previous_timestep = timestep_object;
    }

    // Add periodic time
    first_timestep->previous = previous_timestep;
    previous_timestep->next = first_timestep;

    for(int cpu=0; cpu<cpus; cpu++) {
        sprintf(filename,"%s/movie_files/movie%04d.bin",state_folder.c_str(),cpu);
        ifstream file (filename, ios::in | ios::binary);
        Timestep *timestep_object = first_timestep;
        for(int timestep=0; timestep < timesteps; timestep++) {
            file.read (reinterpret_cast<char*>(&actual_movie_molecules), sizeof(int));
            if(data.size() < 3*actual_movie_molecules) data.resize(6*actual_movie_molecules);
            file.read (reinterpret_cast<char*>(&data[0]), 6*actual_movie_molecules*sizeof(float));
            timestep_object->add_molecule_data(data);
            timestep_object = timestep_object->next;
        }
        file.close();
    }

    Timestep *timestep_object = first_timestep;
    for(int timestep=0; timestep < timesteps; timestep++) {
        timestep_object = timestep_object->next;
    }

    data.clear();
}
