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
#include <topology.h>

void System::step() {
    steps += 1;
    t += dt;
    accelerate();
    move();
    if(topology->num_processors>1) mpi_move();
    timer->start_moving();
    update_molecule_cells();
    timer->end_moving();

    collide();
    if(settings->maintain_pressure) maintain_pressure();
}

void System::mpi_move() {
    timer->start_mpi();

    // MPI_Barrier(MPI_COMM_WORLD);
    for(int i=0; i<topology->num_processors; i++) node_num_new_molecules[i] = 0;

    for(unsigned long n=0; n<num_molecules_local; n++) {
        int node_id = topology->index_from_position(&r[3*n]);
        if(node_id != myid) {
            node_molecule_data[node_id][6*node_num_new_molecules[node_id] + 0] = r[3*n+0];
            node_molecule_data[node_id][6*node_num_new_molecules[node_id] + 1] = r[3*n+1];
            node_molecule_data[node_id][6*node_num_new_molecules[node_id] + 2] = r[3*n+2];
            node_molecule_data[node_id][6*node_num_new_molecules[node_id] + 3] = v[3*n+0];
            node_molecule_data[node_id][6*node_num_new_molecules[node_id] + 4] = v[3*n+1];
            node_molecule_data[node_id][6*node_num_new_molecules[node_id] + 5] = v[3*n+2];
            node_num_new_molecules[node_id]++;
            remove_molecule_from_system(n);
            n--;
        }
    }

    MPI_Status status;

    if(myid>0) {
        // Take all the other nodes first. Node 0 is slower.
        for(int node_id=1; node_id<topology->num_processors; node_id++) {
            if(node_id == myid) continue;
            int num_recieve;
            if(node_id < myid) {
                MPI_Send(&node_num_new_molecules[node_id], 1, MPI_INT, node_id, 100, MPI_COMM_WORLD);
                MPI_Recv(&num_recieve, 1, MPI_INT, node_id, 100, MPI_COMM_WORLD, &status);
                // 6 doubles per molecule
                MPI_Send(&node_molecule_data[node_id][0], 6*node_num_new_molecules[node_id],MPI_DOUBLE,node_id,100, MPI_COMM_WORLD);
                MPI_Recv(&mpi_receive_buffer[0], 6*num_recieve, MPI_DOUBLE, node_id, 100, MPI_COMM_WORLD, &status);
            } else {
                MPI_Recv(&num_recieve, 1, MPI_INT, node_id, 100, MPI_COMM_WORLD, &status);
                MPI_Send(&node_num_new_molecules[node_id], 1, MPI_INT, node_id, 100, MPI_COMM_WORLD);
                // 6 doubles per molecule
                MPI_Recv(&mpi_receive_buffer[0], 6*num_recieve, MPI_DOUBLE, node_id, 100, MPI_COMM_WORLD, &status);
                MPI_Send(&node_molecule_data[node_id][0], 6*node_num_new_molecules[node_id],MPI_DOUBLE,node_id,100, MPI_COMM_WORLD);
            }

            add_molecules_from_mpi(mpi_receive_buffer, num_recieve);
        }

        // Then communicate with node 0
        int num_recieve;
        MPI_Send(&node_num_new_molecules[0], 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
        MPI_Recv(&num_recieve, 1, MPI_INT, 0, 100, MPI_COMM_WORLD, &status);
        // 6 doubles per molecule
        MPI_Send(&node_molecule_data[0][0], 6*node_num_new_molecules[0],MPI_DOUBLE,0,100, MPI_COMM_WORLD);
        MPI_Recv(&mpi_receive_buffer[0], 6*num_recieve, MPI_DOUBLE, 0, 100, MPI_COMM_WORLD, &status);
        add_molecules_from_mpi(mpi_receive_buffer, num_recieve);
    } else {
        for(int node_id=1; node_id<topology->num_processors; node_id++) {
            int num_recieve;
            MPI_Recv(&num_recieve, 1, MPI_INT, node_id, 100, MPI_COMM_WORLD, &status);
            MPI_Send(&node_num_new_molecules[node_id], 1, MPI_INT, node_id, 100, MPI_COMM_WORLD);
            // 6 doubles per molecule
            MPI_Recv(&mpi_receive_buffer[0], 6*num_recieve, MPI_DOUBLE, node_id, 100, MPI_COMM_WORLD, &status);
            MPI_Send(&node_molecule_data[node_id][0], 6*node_num_new_molecules[node_id],MPI_DOUBLE,node_id,100, MPI_COMM_WORLD);
            add_molecules_from_mpi(mpi_receive_buffer, num_recieve);
        }
    }

    num_molecules_global = 0;
    MPI_Allreduce(&num_molecules_local, &num_molecules_global, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD) ;
    timer->end_mpi();
}

void System::move() {
    timer->start_moving();
    for(int n=0;n<num_molecules_local;n++) {
        // mover->move_molecule_cylinder(n,dt,rnd,0);
        mover->move_molecule(n,dt,rnd,0);
        // mover->move_molecule_box(n,dt,rnd,0);
    }
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
    if(settings->flow_direction < 0 || settings->gravity == 0) return;
    timer->start_accelerate();

    int flow_dir = settings->flow_direction;
    double gravity = settings->gravity*dt;

    for(int n=0;n<num_molecules_local;n++) {
        v[3*n+flow_dir] += gravity;
    }
    timer->end_accelerate();
}

void System::update_molecule_cells() {
    for(int n=0;n<num_molecules_local;n++) {
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

void System::add_molecule_to_cell(Cell *cell, const int &molecule_index) {
    cell->add_molecule(molecule_index,molecule_index_in_cell,molecule_cell_index);
    num_molecules_local++;
}

void System::add_molecules_from_mpi(vector<double> &data, const int &num_new_molecules) {
    for(int i=0; i<num_new_molecules; i++) {
        int n = num_molecules_local;
        r[3*n+0] = data[6*i+0];
        r[3*n+1] = data[6*i+1];
        r[3*n+2] = data[6*i+2];
        r0[3*n+0] = data[6*i+0];
        r0[3*n+1] = data[6*i+1];
        r0[3*n+2] = data[6*i+2];
        v[3*n+0] = data[6*i+3];
        v[3*n+1] = data[6*i+4];
        v[3*n+2] = data[6*i+5];

        Cell *cell = all_cells[cell_index_from_position(&r[3*n])];
        cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
        num_molecules_local++;
    }
}

void System::remove_molecule_from_system(const long &molecule_index) {
    long cell_index = molecule_cell_index[molecule_index];
    Cell *cell = all_cells[cell_index];
    cell->remove_molecule(molecule_index,molecule_index_in_cell);

    // Move the last molecule into that memory location
    int last_molecule_index = num_molecules_local-1;

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
    num_molecules_local--;
}

void System::add_molecules_in_inlet_reservoir(Cell *cell, const double &velocity_std_dev, const int &delta_num_molecules) {
    int neighbor_cell_index_vector[3];
    neighbor_cell_index_vector[0] = cell->index_vector[0];
    neighbor_cell_index_vector[1] = cell->index_vector[1];
    neighbor_cell_index_vector[2] = cell->index_vector[2];
    neighbor_cell_index_vector[settings->flow_direction] += 1; // We want the next cell in flow direction
    int neighbor_cell_index = cell_index_from_ijk(neighbor_cell_index_vector[0], neighbor_cell_index_vector[1], neighbor_cell_index_vector[2]);
    vector<double> average_velocity_neighbor_cell = ((Cell*)all_cells[neighbor_cell_index])->update_average_velocity();

    for(int i=0; i<delta_num_molecules; i++) {
        int molecule_index = num_molecules_local;

        v[3*molecule_index+0] = rnd->next_gauss()*velocity_std_dev + average_velocity_neighbor_cell[0];
        v[3*molecule_index+1] = rnd->next_gauss()*velocity_std_dev + average_velocity_neighbor_cell[1];
        v[3*molecule_index+2] = rnd->next_gauss()*velocity_std_dev + average_velocity_neighbor_cell[2];
        find_position_in_cell(cell, &r[3*molecule_index]);
        r0[3*molecule_index+0] = r[3*molecule_index+0];
        r0[3*molecule_index+1] = r[3*molecule_index+1];
        r0[3*molecule_index+2] = r[3*molecule_index+2];

        add_molecule_to_cell(cell, molecule_index);
    }
}

void System::remove_molecules_in_inlet_reservoir(Cell *cell, const int &delta_num_molecules) {
    // delta_num_molecules is a negative number
    for(int i=0; i<abs(delta_num_molecules); i++) {
        int this_molecule_index_in_cell = cell->num_molecules*rnd->next_double();
        long molecule_index = cell->molecules[this_molecule_index_in_cell];
        remove_molecule_from_system(molecule_index);
    }
}

void System::add_molecules_in_outlet_reservoir(Cell *cell, const double &velocity_std_dev, const int &delta_num_molecules) {
    int neighbor_cell_index_vector[3];
    neighbor_cell_index_vector[0] = cell->index_vector[0];
    neighbor_cell_index_vector[1] = cell->index_vector[1];
    neighbor_cell_index_vector[2] = cell->index_vector[2];
    neighbor_cell_index_vector[settings->flow_direction] -= 1; // We want the previous cell in flow direction
    int neighbor_cell_index = cell_index_from_ijk(neighbor_cell_index_vector[0], neighbor_cell_index_vector[1], neighbor_cell_index_vector[2]);
    vector<double> average_velocity_neighbor_cell = ((Cell*)all_cells[neighbor_cell_index])->update_average_velocity();

    for(int i=0; i<delta_num_molecules; i++) {
        int molecule_index = num_molecules_local;

        v[3*molecule_index+0] = rnd->next_gauss()*velocity_std_dev + average_velocity_neighbor_cell[0];
        v[3*molecule_index+1] = rnd->next_gauss()*velocity_std_dev + average_velocity_neighbor_cell[1];
        v[3*molecule_index+2] = rnd->next_gauss()*velocity_std_dev + average_velocity_neighbor_cell[2];
        find_position_in_cell(cell, &r[3*molecule_index]);
        r0[3*molecule_index+0] = r[3*molecule_index+0];
        r0[3*molecule_index+1] = r[3*molecule_index+1];
        r0[3*molecule_index+2] = r[3*molecule_index+2];

        add_molecule_to_cell(cell, molecule_index);
    }
}

void System::remove_molecules_in_outlet_reservoir(Cell *cell, const int &delta_num_molecules) {
    // delta_num_molecules is a negative number
    for(int i=0; i<abs(delta_num_molecules); i++) {
        int this_molecule_index_in_cell = cell->num_molecules*rnd->next_double();
        long molecule_index = cell->molecules[this_molecule_index_in_cell];
        remove_molecule_from_system(molecule_index);
    }
}

void System::maintain_pressure_A() {
    /*
     * Strategy: Keep the number of molecules in a cell equal to that of the neighbor cell towards the center of the system
     *           If adding any new molecules, let the velocity be equal to the average in the neighbor cell plus a MB-distribution.
     */

    double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_A);
    double wanted_density = wanted_pressure / temperature;

    double velocity_std_dev = sqrt(wanted_pressure/density); // Shouldn't this density be the wanted density?
    velocity_std_dev = sqrt(temperature/settings->mass);

    for(int i=0;i<reservoir_A_cells.size();i++) {
        Cell *cell = reservoir_A_cells[i];
        if(cell->volume == 0) continue;
        int num_molecules_in_reservoir = cell->num_molecules;
        int wanted_num_molecules_in_reservoir = wanted_density*cell->volume / atoms_per_molecule;
        int delta_num_molecules = wanted_num_molecules_in_reservoir - num_molecules_in_reservoir;
        // cout << "I want " << wanted_num_molecules_in_reservoir << " molecules, and I have " << num_molecules_in_reservoir << endl;

        if(delta_num_molecules > 0) {
            add_molecules_in_inlet_reservoir(cell, velocity_std_dev, delta_num_molecules);
        } else if(delta_num_molecules < 0) {
            remove_molecules_in_inlet_reservoir(cell, delta_num_molecules);
        }
    }
}

void System::maintain_pressure_B() {
    /*
     * Strategy: Keep the number of molecules in a cell equal to that of the neighbor cell towards the center of the system
     *           If adding any new molecules, let the velocity be equal to the average in the neighbor cell plus a MB-distribution.
     */

    double wanted_pressure = unit_converter->pressure_from_SI(settings->pressure_B);
    double velocity_std_dev = sqrt(wanted_pressure/density); // Shouldn't this density be the wanted density?
    velocity_std_dev = sqrt(temperature/settings->mass);

    for(int i=0;i<reservoir_B_cells.size();i++) {
        Cell *cell = reservoir_B_cells[i];
        if(cell->volume == 0) continue;

        // Find neighbor cell
        int neighbor_cell_index_vector[3];
        neighbor_cell_index_vector[0] = cell->index_vector[0];
        neighbor_cell_index_vector[1] = cell->index_vector[1];
        neighbor_cell_index_vector[2] = cell->index_vector[2];
        neighbor_cell_index_vector[settings->flow_direction] -= 1; // We want the previous cell in flow direction
        int neighbor_cell_index = cell_index_from_ijk(neighbor_cell_index_vector[0], neighbor_cell_index_vector[1], neighbor_cell_index_vector[2]);
        Cell *neighbor_cell = all_cells[neighbor_cell_index];

        int num_molecules_in_reservoir = cell->num_molecules;
        int num_molecules_in_neighbor_cell = neighbor_cell->num_molecules;
        // cout << "Neighbor cell has " << num_molecules_in_neighbor_cell << " molecules, and I have " << num_molecules_in_reservoir << endl;

        int wanted_num_molecules_in_reservoir = num_molecules_in_neighbor_cell;
        int delta_num_molecules = wanted_num_molecules_in_reservoir - num_molecules_in_reservoir;

        if(delta_num_molecules > 0) {
            add_molecules_in_outlet_reservoir(cell, velocity_std_dev, delta_num_molecules);
        } else if(delta_num_molecules < 0) {
            remove_molecules_in_outlet_reservoir(cell, delta_num_molecules);
        }
    }
}

void System::maintain_pressure() {
    timer->start_pressure();
    maintain_pressure_A();
    maintain_pressure_B();
    timer->end_pressure();
}
