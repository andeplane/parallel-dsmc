#include <system.h>

#include <math.h>
#include <fstream>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <unitconverter.h>
#include <dsmc_io.h>

#include <time.h>
#include <system.inc.cpp>
#include <dsmctimer.h>

void System::step() {
    steps += 1;
    t += dt;
    accelerate();
    move();
    collide();
    if(settings->maintain_pressure) maintain_pressure();
}

void System::move() {
    timer->start_moving();
    for(int n=0;n<num_molecules;n++) {
        double v0 = v[3*n+2];
        // mover->move_molecule_cylinder(n,dt,rnd,0);
        mover->move_molecule(n,dt,rnd,0);
    }
    timer->end_moving();

    timer->start_moving();
    if(myid==0) update_molecule_cells();
    timer->end_moving();
}

void System::collide() {
    timer->start_colliding();

    for(int i=0;i<active_cells.size();i++) {
        Cell *cell = active_cells[i];

        cell->prepare();
        collisions += cell->collide(rnd);
    }

    timer->end_colliding();
}

void System::accelerate() {
    if(settings->gravity_direction < 0) return;
    timer->start_accelerate();

    int gravity_dir = settings->gravity_direction;
    double gravity = settings->gravity*dt;

    for(int n=0;n<num_molecules;n++) {
        v[3*n+gravity_dir] += gravity;
    }
    timer->end_accelerate();
}

void System::find_position_in_reservoirs(double *r, bool find_position_in_A) {
    bool did_collide = true;
    while(did_collide) {
        r[0] = length[0]*rnd->next_double();
        r[1] = length[1]*rnd->next_double();
        r[2] = length[2]*rnd->next_double();
        if(find_position_in_A)  r[settings->gravity_direction] = reservoir_size*rnd->next_double();
        else r[settings->gravity_direction] = (length[settings->gravity_direction] - reservoir_size) + reservoir_size*rnd->next_double();

        did_collide = *world_grid->get_voxel(r)>=voxel_type_wall;
    }
}

void System::add_molecule_in_pressure_reservoirs(bool add_in_A) {
    int n = num_molecules;

    v[3*n+0] = rnd->next_gauss()*sqrt(temperature/settings->mass);
    v[3*n+1] = rnd->next_gauss()*sqrt(temperature/settings->mass);
    v[3*n+2] = rnd->next_gauss()*sqrt(temperature/settings->mass);

    find_position_in_reservoirs(&r[3*n],add_in_A);
    r0[3*n+0] = r[3*n+0];
    r0[3*n+1] = r[3*n+1];
    r0[3*n+2] = r[3*n+2];

    Cell *cell = all_cells[cell_index_from_position(&r[3*n])];
    cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
    num_molecules++;
}

bool System::remove_molecule_in_pressure_reservoir(bool remove_from_A) {
    Cell *cell = NULL;

    if(remove_from_A) cell = reservoir_A_cells[ reservoir_A_cells.size()*rnd->next_double() ];
    else cell = reservoir_B_cells[ reservoir_B_cells.size()*rnd->next_double() ];

    if(cell->num_molecules>0) {
        // Remove this random molecule
        int this_molecule_index_in_cell = cell->num_molecules*rnd->next_double();
        int molecule_index = cell->molecules[this_molecule_index_in_cell];
        cell->remove_molecule(molecule_index,molecule_index_in_cell);

        // Move the last molecule into that memory location
        int last_molecule_index = num_molecules-1;

        while(last_molecule_index==molecule_index) {
            last_molecule_index--;
        }

        int last_molecule_cell_index = molecule_cell_index[last_molecule_index];
        int last_molecule_index_in_cell = molecule_index_in_cell[last_molecule_index];
        memcpy(&r[3*molecule_index],&r[3*last_molecule_index],3*sizeof(double));
        memcpy(&v[3*molecule_index],&v[3*last_molecule_index],3*sizeof(double));
        memcpy(&r0[3*molecule_index],&r0[3*last_molecule_index],3*sizeof(double));

        cell = all_cells[last_molecule_cell_index];
        molecule_cell_index[molecule_index] = last_molecule_cell_index;
        molecule_index_in_cell[molecule_index] = last_molecule_index_in_cell;
        cell->molecules[last_molecule_index_in_cell] = molecule_index;
        num_molecules--;
        return true;
    } else return false;
}

void System::update_molecule_cells() {
    for(int n=0;n<num_molecules;n++) {
        int cell_index_new = cell_index_from_position(&r[3*n]);
        int cell_index_old = molecule_cell_index[n];
//        cout << "New cell index: " << cell_index_new << endl;
//        cout << "Old cell index: " << cell_index_old << endl;
//        cout << "Position: " << r[3*n + 0] << " " << r[3*n + 1] << " " << r[3*n + 2] << endl;

        Cell *new_cell = all_cells[cell_index_new];
        Cell *old_cell = all_cells[cell_index_old];

        if(cell_index_new != cell_index_old) {
            // We changed cell
            old_cell->remove_molecule(n,molecule_index_in_cell);
            new_cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
        }
    }
}

void System::maintain_pressure() {
    timer->start_pressure();
    maintain_pressure_A();
    maintain_pressure_B();
    timer->end_pressure();
}

void System::maintain_pressure_A() {
    long num_molecules_in_reservoir = 0;
    double volume_in_reservoir = 0;
    double pressure_in_reservoir = 0;
    double kinetic_energy_in_reservoir = 0;

    for(int i=0;i<reservoir_A_cells.size();i++) {
        Cell *cell = reservoir_A_cells[i];
        num_molecules_in_reservoir += cell->num_molecules;
        volume_in_reservoir += cell->volume;
        kinetic_energy_in_reservoir += cell->calculate_kinetic_energy();
    }

    double temperature_in_reservoir = 2.0/3.0*kinetic_energy_in_reservoir/(num_molecules_in_reservoir*atoms_per_molecule);
    double temp_over_volume = 0;
    if(volume_in_reservoir>0) {
        temp_over_volume = temperature_in_reservoir/volume_in_reservoir;
        pressure_in_reservoir = atoms_per_molecule*num_molecules_in_reservoir*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_A);
        long wanted_num_molecules = wanted_pressure*volume_in_reservoir/temperature_in_reservoir/atoms_per_molecule;

        long delta = wanted_num_molecules-num_molecules_in_reservoir;

        if(pressure_in_reservoir<wanted_pressure) {
            int num_add = abs(delta)*0.1;
            for(int i=0;i<num_add;i++) {
                add_molecule_in_pressure_reservoirs(true);
            }
        } else {
            int num_remove = abs(delta)*0.1;
            for(int i=0;i<num_remove;i++) {
                if(remove_molecule_in_pressure_reservoir(true)) { }
                else i--;
            }
        }
    }
}

void System::maintain_pressure_B() {
    long num_molecules_in_reservoir = 0;
    double volume_in_reservoir = 0;
    double pressure_in_reservoir = 0;
    double kinetic_energy_in_reservoir = 0;

    for(int i=0;i<reservoir_B_cells.size();i++) {
        Cell *cell = reservoir_B_cells[i];
        num_molecules_in_reservoir += cell->num_molecules;
        volume_in_reservoir += cell->volume;
        kinetic_energy_in_reservoir += cell->calculate_kinetic_energy();
    }
    double temperature_in_reservoir = 2.0/3.0*kinetic_energy_in_reservoir/(num_molecules_in_reservoir*atoms_per_molecule);

    double temp_over_volume = 0;
    if(volume_in_reservoir>0) {
        temp_over_volume = temperature_in_reservoir/volume_in_reservoir;
        pressure_in_reservoir = atoms_per_molecule*num_molecules_in_reservoir*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_B);
        long wanted_num_molecules = wanted_pressure*volume_in_reservoir/temperature_in_reservoir/atoms_per_molecule;
        // Delta is the actual delta number of particles, but in order to smooth out the particle adding, we only fix 10% each timestep
        long delta = wanted_num_molecules-num_molecules_in_reservoir;

        int delta_num_particles = 0;
        int num_removed_particles = 0;
        if(pressure_in_reservoir<wanted_pressure) {
            int num_add = abs(delta)*0.1;
            // Added particles did not go through boundary
            delta_num_particles = num_add;
            for(int i=0;i<num_add;i++) {
                add_molecule_in_pressure_reservoirs(false);
            }
        } else {
            num_removed_particles = abs(delta)*0.1;
            // Removed particles might have gone through boundary
            delta_num_particles = -num_removed_particles ;
            for(int i=0;i<num_removed_particles;i++) {
                if(remove_molecule_in_pressure_reservoir(false)) { }
                else i--;
            }
        }

        int current_number_of_particles = num_molecules_in_reservoir - delta_num_particles;
        int change_in_particle_count_since_last_measurement = current_number_of_particles - reservoir_b_particle_count;

        // Removed particles are seen as moved to another part of the reservoir, they should still be counted as in the reservoir.
        flux_count += change_in_particle_count_since_last_measurement + num_removed_particles;
        reservoir_b_particle_count = current_number_of_particles;
    }
}
