#include <cstddef>
#include <cell.h>
#include <math.h>
#include <time.h>
#include <system.h>
#include <random.h>
#include <settings.h>
#include <cstring>

using namespace std;

Cell::Cell(System *_system) {
    system = _system;
    vr_max = 0;
    pixels = 0;
    num_molecules = 0;
    total_pixels = 0;
    collision_rest = 0;
    average_velocity.resize(3,0);

    num_molecules_allocated_memory = 0;
    molecules = NULL;
}

bool Cell::cmp(Cell *c1, Cell *c2) {
    return c1->collision_pairs < c2->collision_pairs;
}

void Cell::update_volume() {
    if(pixels==0) {
        // cout << "No pixels in here?" << endl;
        volume = 0;
        collision_coefficient = 0;
        return;
    }

    num_molecules_allocated_memory = 100;
    molecules = new int[num_molecules_allocated_memory];
    // Update the effective cell volume. A cell may contain 50% of solid material
    volume = system->volume/(system->cells_x*system->cells_y*system->cells_z)*(float)pixels/total_pixels;
    collision_coefficient = 0.5*system->atoms_per_molecule*M_PI*system->diam*system->diam*system->dt/volume;
}

unsigned long Cell::prepare() {
    //* Determine number of candidate collision pairs to be selected in this cell
    //* Add the rest from previous timestep
    double select = collision_coefficient*num_molecules*(num_molecules-1)*vr_max + collision_rest;
    
    collision_pairs = round(select);      // Number of pairs to be selected
    collision_rest = select - collision_pairs;
    // collision_pairs = 0;

    return collision_pairs;
}

void Cell::collide_molecules(double *v0, double *v1, const double &v_rel, Random *rnd) {
    double vcmx  = 0.5*(v0[0] + v1[0]);
    double vcmy  = 0.5*(v0[1] + v1[1]);
    double vcmz  = 0.5*(v0[2] + v1[2]);

    double cos_th = 1.0 - 2.0*rnd->next_double();      // Cosine and sine of
    double sin_th = sqrt(1.0 - cos_th*cos_th);        // collision angle theta
    double phi = 2*M_PI*rnd->next_double();

    double vrelx = v_rel*cos_th;                   // Compute post-collision relative velocity
    double vrely = v_rel*sin_th*cos(phi);
    double vrelz = v_rel*sin_th*sin(phi);

    v0[0] = vcmx + 0.5*vrelx;
    v0[1] = vcmy + 0.5*vrely;
    v0[2] = vcmz + 0.5*vrelz;

    v1[0] = vcmx - 0.5*vrelx;
    v1[1] = vcmy - 0.5*vrely;
    v1[2] = vcmz - 0.5*vrelz;
}

int Cell::collide(Random *rnd) {

    //* Skip cells with only one particle

    if( num_molecules < 2 ) return 0;  // Skip to the next cell

    double crm = vr_max;     // Current maximum relative speed

	//* Loop over total number of candidate collision pairs
    int collisions = 0;
    for(int isel=0; isel<collision_pairs; isel++ ) {
		//* Pick two particles at random out of this cell
        int ip0 = (int)(rnd->next_double()*num_molecules);
        int ip1 = ((int)(ip0+1+rnd->next_double()*(num_molecules-1))) % num_molecules;

        ip0 = molecules[ip0];
        ip1 = molecules[ip1];

		//* Calculate pair's relative speed
        double *v0 = &system->v[3*ip0];
        double *v1 = &system->v[3*ip1];
        double v_rel = sqrt(pow(v0[0] - v1[0],2) + pow(v0[1] - v1[1],2) + pow(v0[2] - v1[2],2));

        if( v_rel > crm ) {         // If relative speed larger than crm,
            crm = v_rel;            // then reset crm to larger value
        }

		//* Accept or reject candidate pair according to relative speed
        if( v_rel > rnd->next_double()*vr_max ) {
			//* If pair accepted, select post-collision velocities

			collisions++;
            collide_molecules(v0, v1, v_rel, rnd);
		} // Loop over pairs
	}
	
    vr_max = crm;
	return collisions;
}

void Cell::add_molecule(const int &molecule_index, unsigned long *index_in_cell, unsigned long *cell_index) {
    if(num_molecules+1>num_molecules_allocated_memory) {
        // We need to reallocate
        num_molecules_allocated_memory *= 2;
        int *tmp = new int[num_molecules_allocated_memory];
        memcpy(tmp,molecules,num_molecules*sizeof(int));
        delete molecules;
        molecules = tmp;
    }

    molecules[num_molecules] = molecule_index;
    index_in_cell[molecule_index] = num_molecules;
    cell_index[molecule_index] = index;
    num_molecules++;
}

void Cell::remove_molecule(const int &molecule_index, unsigned long *index_in_cell) {
    if(num_molecules>1) {
        // Move the last molecule over here
        molecules[ index_in_cell[molecule_index] ] = molecules[num_molecules-1];
        index_in_cell[ molecules[num_molecules-1] ] = index_in_cell[molecule_index];
    }

    num_molecules--;
}

double Cell::calculate_kinetic_energy() {
    double kinetic_energy = 0;
    for(int n=0;n<num_molecules;n++) {
        int index = molecules[n];
        kinetic_energy += system->v[3*index+0]*system->v[3*index+0] + system->v[3*index+1]*system->v[3*index+1] + system->v[3*index+2]*system->v[3*index+2];
    }
    kinetic_energy *= 0.5*system->settings->mass*system->atoms_per_molecule;

    return kinetic_energy;
}

vector<double>& Cell::update_average_velocity() {
    vector<double> new_average_velocity(3,0);

    for(int n=0; n<num_molecules; n++) {
        int index = molecules[n];
        new_average_velocity[0] += system->v[3*index+0];
        new_average_velocity[1] += system->v[3*index+1];
        new_average_velocity[2] += system->v[3*index+2];
    }

    if(num_molecules > 0) {
        new_average_velocity[0] /= num_molecules;
        new_average_velocity[1] /= num_molecules;
        new_average_velocity[2] /= num_molecules;
    }

    average_velocity[0] = 0.8*average_velocity[0] + 0.2*new_average_velocity[0];
    average_velocity[1] = 0.8*average_velocity[1] + 0.2*new_average_velocity[1];
    average_velocity[2] = 0.8*average_velocity[2] + 0.2*new_average_velocity[2];

    return average_velocity;
}
