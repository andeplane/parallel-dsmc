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
    if(myid==0) accelerate();
    move();
    if(myid==0) collide();
    if(myid==0 && settings->maintain_pressure) maintain_pressure();
}

void System::send_molecules_to_slaves() {

    int molecules_per_node = num_molecules/num_nodes;
    for(int node_id=1;node_id<num_nodes;node_id++) {
        int start_index = molecules_per_node*node_id;
        int end_index = start_index + molecules_per_node;
        if(node_id == num_nodes-1) end_index = num_molecules - 1;
        end_index = min(end_index,num_molecules - 1);

        int num_send = end_index - start_index + 1;

        MPI_Send(&num_send,1,MPI_INT,node_id,10,MPI_COMM_WORLD);
        MPI_Send(&r[3*start_index],3*num_send,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD);
        MPI_Send(&v[3*start_index],3*num_send,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD);
        MPI_Send(&r0[3*start_index],3*num_send,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD);
    }

    num_molecules_this_node = molecules_per_node;
}

void System::receive_molecules_from_master() {
    MPI_Status status;

    MPI_Recv(&num_molecules,1,MPI_INT,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(r,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(v,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
    MPI_Recv(r0,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD,&status);
    num_molecules_this_node = num_molecules;
}

void System::send_molecules_to_master() {
    MPI_Send(r,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
    MPI_Send(v,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
    MPI_Send(r0,3*num_molecules,MPI_DOUBLE,0,10,MPI_COMM_WORLD);
}

void System::receive_molecules_from_slaves() {
    MPI_Status status;

    int molecules_per_node = num_molecules/num_nodes;
    for(int node_id=1;node_id<num_nodes;node_id++) {
        int start_index = molecules_per_node*node_id;
        int end_index = start_index + molecules_per_node;
        if(node_id == num_nodes-1) end_index = num_molecules - 1;
        end_index = min(end_index,num_molecules - 1);

        int num_receive = end_index - start_index + 1;

        MPI_Recv(&r[3*start_index],3*num_receive,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD,&status);
        MPI_Recv(&v[3*start_index],3*num_receive,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD,&status);
        MPI_Recv(&r0[3*start_index],3*num_receive,MPI_DOUBLE,node_id,10,MPI_COMM_WORLD,&status);
    }
}

void System::move() {
    timer->start_mpi();
    if(myid==0) send_molecules_to_slaves();
    else receive_molecules_from_master();
    timer->end_mpi();

    timer->start_moving();
    // cout << myid << " will move " << num_molecules_this_node << " molecules." << endl;
    for(int n=0;n<num_molecules_this_node;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }
    timer->end_moving();

    timer->start_mpi();
    if(myid==0) receive_molecules_from_slaves();
    else send_molecules_to_master();
    timer->end_mpi();

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
        r[0] = length[0]*rnd->nextDouble();
        r[1] = length[1]*rnd->nextDouble();
        r[2] = length[2]*rnd->nextDouble();
        if(find_position_in_A)  r[settings->gravity_direction] = reservoir_size*rnd->nextDouble();
        else r[settings->gravity_direction] = (length[settings->gravity_direction] - reservoir_size) + reservoir_size*rnd->nextDouble();

        did_collide = *world_grid->get_voxel(r)>=voxel_type_wall;
    }
}

void System::add_molecule_in_pressure_reservoirs(bool add_in_A) {
    int n = num_molecules;

    v[3*n+0] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    v[3*n+1] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    v[3*n+2] = rnd->nextGauss()*sqrt(temperature/settings->mass);

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

    if(remove_from_A) cell = reservoir_A_cells[ reservoir_A_cells.size()*rnd->nextDouble() ];
    else cell = reservoir_B_cells[ reservoir_B_cells.size()*rnd->nextDouble() ];

    if(cell->num_molecules>0) {
        // Remove this random molecule
        int this_molecule_index_in_cell = cell->num_molecules*rnd->nextDouble();
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
    double max_ke = 0;
    for(int i=0;i<reservoir_A_cells.size();i++) {
        Cell *cell = reservoir_A_cells[i];
        num_molecules_in_reservoir += cell->num_molecules;
        volume_in_reservoir += cell->volume;
        double ke = cell->calculate_kinetic_energy();
        if(ke>max_ke) max_ke=ke;
        kinetic_energy_in_reservoir += ke;
    }

    double temperature_in_reservoir = 2.0/3*kinetic_energy_in_reservoir/(num_molecules_in_reservoir*atoms_per_molecule);

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
    double temperature_in_reservoir = 2.0/3*kinetic_energy_in_reservoir/num_molecules_in_reservoir;

    double temp_over_volume = 0;
    if(volume_in_reservoir>0) {
        temp_over_volume = temperature_in_reservoir/volume_in_reservoir;
        pressure_in_reservoir = atoms_per_molecule*num_molecules_in_reservoir*temp_over_volume;
        double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_B);
        long wanted_num_molecules = wanted_pressure*volume_in_reservoir/temperature_in_reservoir/atoms_per_molecule;
        long delta = wanted_num_molecules-num_molecules_in_reservoir;

        if(pressure_in_reservoir<wanted_pressure) {
            int num_add = abs(delta)*0.1;
            for(int i=0;i<num_add;i++) {
                add_molecule_in_pressure_reservoirs(false);
            }
        } else {
            int num_remove = abs(delta)*0.1;
            for(int i=0;i<num_remove;i++) {
                if(remove_molecule_in_pressure_reservoir(false)) { }
                else i--;
            }
        }
    }
}
