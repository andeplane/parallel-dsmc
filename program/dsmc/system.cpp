#include <system.h>

#include <cmath>
#include <fstream>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <unitconverter.h>
#include <dsmc_io.h>
#include <settings.h>
#include <moleculemover.h>
#include <time.h>
#include <dsmctimer.h>
#include <topology.h>
#include <colliderbase.h>
#include <colliderspecular.h>
#include <colliderthermal.h>
#include <collidercercignanilampis.h>
#include <collidermaxwell.h>
#include <cvector.h>

void System::step() {
    for(int n=0; n<num_molecules_local; n++) steps_since_collision[n]++;
    steps += 1;
    t += dt;
    string step_state = " accelerate()";

    try {
        accelerate();
        step_state = " move()";
        move();
        step_state = " mpi_move()";
        if(topology->num_processors>1) mpi_move();
        timer->start_moving();
        step_state = " update_molecule_cells()";
        update_molecule_cells();
        timer->end_moving();

        collide();
    } catch (string ex) {
        cout << "Got exception: " << ex << endl;
        cout << "while doing " << step_state << endl;
    }

    // if(settings->maintain_pressure) maintain_pressure();
}

void System::mpi_move() {
    timer->start_mpi();
    for(int i=0; i<topology->num_processors; i++) node_num_new_molecules.at(i) = 0;

    for(unsigned long n=0; n<num_molecules_local; n++) {
        int node_id = topology->index_from_molecule_index(n);
        if(node_id != myid) {
            node_id = topology->facet_id_to_node_id_list[topology->node_id_to_facet_id_list[node_id]];
            if(node_id < 0 || node_id >= topology->num_processors) cout << "Moved existing molecule to a node that doesn't exists" << endl;
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 0) = r.at(3*n+0);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 1) = r.at(3*n+1);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 2) = r.at(3*n+2);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 3) = v.at(3*n+0);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 4) = v.at(3*n+1);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 5) = v.at(3*n+2);
            node_num_new_molecules.at(node_id)++;
            remove_molecule_from_system(n);
            n--;
        }
    }

    MPI_Status status;
    int num_recieve;
    int num_send;
    for(int dimension = 0; dimension < 3; dimension++) {
        for(int upper = 0; upper<=1; upper++) {
            int facet_index = 2*dimension + upper;
            int node_id = topology->facet_id_to_node_id_list[facet_index];
            if(node_id == myid) continue;
            num_send = node_num_new_molecules.at(node_id);
            if(topology->my_parity[dimension] == 0) {
                MPI_Send(&num_send, 1, MPI_INT, node_id, 10+facet_index, MPI_COMM_WORLD);
                MPI_Recv(&num_recieve, 1, MPI_INT, MPI_ANY_SOURCE, 10+facet_index, MPI_COMM_WORLD, &status);
                // 6 doubles per molecule
                if(num_send) MPI_Send(&node_molecule_data.at(node_id).front(), 6*num_send,mpi_data_type,node_id,100+facet_index, MPI_COMM_WORLD);
                if(num_recieve) MPI_Recv(&mpi_receive_buffer.front(), 6*num_recieve, mpi_data_type, MPI_ANY_SOURCE, 100+facet_index, MPI_COMM_WORLD, &status);
            } else if (topology->my_parity[dimension] == 1){
                MPI_Recv(&num_recieve, 1, MPI_INT, MPI_ANY_SOURCE, 10+facet_index, MPI_COMM_WORLD, &status);
                MPI_Send(&num_send, 1, MPI_INT, node_id, 10+facet_index, MPI_COMM_WORLD);
                // 6 doubles per molecule
                if(num_recieve) MPI_Recv(&mpi_receive_buffer.front(), 6*num_recieve, mpi_data_type, MPI_ANY_SOURCE, 100+facet_index, MPI_COMM_WORLD, &status);
                if(num_send) MPI_Send(&node_molecule_data.at(node_id).front(), 6*num_send,mpi_data_type,node_id,100+facet_index, MPI_COMM_WORLD);
            }

            node_num_new_molecules.at(node_id) = 0; // We have sent everything we wanted to this node now.
            add_molecules_from_mpi(mpi_receive_buffer, num_recieve);
        }

    }

    num_molecules_global = 0;
    timer->end_mpi();

    timer->start_mpi_reduce();
    MPI_Reduce(&num_molecules_local, &num_molecules_global, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD) ;
    timer->end_mpi_reduce();
}

void System::move() {
    timer->start_moving();
    CVector system_length(length[0], length[1], length[2]);
    #pragma simd
    for(int n=0;n<num_molecules_local;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }

    for(int n=0;n<num_molecules_local;n++) {
        mover->apply_periodic_boundary_conditions(n, r, system_length);
    }
    timer->end_moving();
}

void System::collide() {
    timer->start_colliding();
    for(int i=0;i<active_cells.size();i++) {
        Cell *cell = active_cells.at(i);

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
        v.at(3*n+flow_dir) += gravity;
    }
    timer->end_accelerate();
}

void System::update_molecule_cells() {
    for(int n=0;n<num_molecules_local;n++) {
        int cell_index_new = cell_index_map[cell_index_from_position(n)];
        int cell_index_old = cell_index_map[molecule_cell_index.at(n)];

        Cell *new_cell;
        Cell *old_cell;
        try {
            new_cell = active_cells.at(cell_index_new);
            old_cell = active_cells.at(cell_index_old);
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
            cout << "Error in update_molecule_cells(), molecule wants a cell that is outside the range" << endl;
        }

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

void System::add_molecules_from_mpi(vector<data_type> &data, const int &num_new_molecules) {
    for(int i=0; i<num_new_molecules; i++) {
        int node_id = topology->index_from_position(&data[6*i+0]);
        if(node_id != myid) {
            if(node_id < 0 || node_id >= topology->num_processors) throw string("Moved new molecule to a node that doesn't exists");
            node_id = topology->facet_id_to_node_id_list[topology->node_id_to_facet_id_list[node_id]];
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 0) = data.at(6*i+0);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 1) = data.at(6*i+1);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 2) = data.at(6*i+2);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 3) = data.at(6*i+3);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 4) = data.at(6*i+4);
            node_molecule_data.at(node_id).at(6*node_num_new_molecules.at(node_id) + 5) = data.at(6*i+5);
            node_num_new_molecules.at(node_id)++;
            continue;
        }
        int n = num_molecules_local;
        r.at(3*n+0) = data.at(6*i+0);
        r.at(3*n+1) = data.at(6*i+1);
        r.at(3*n+2) = data.at(6*i+2);
        v.at(3*n+0) = data.at(6*i+3);
        v.at(3*n+1) = data.at(6*i+4);
        v.at(3*n+2) = data.at(6*i+5);
        int cell_index = cell_index_map[cell_index_from_position(n)];
        Cell *cell;
        try {
            cell = active_cells.at(cell_index);
        } catch (const std::out_of_range& oor) {
            std::cerr << "Out of Range error: " << oor.what() << '\n';
            cout << "New molecule from another node didn't land inside a cell that exists." << endl;
        }
        cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
        num_molecules_local++;
    }
}

void System::remove_molecule_from_system(const long &molecule_index) {
    Cell *cell = cell_containing_molecule(molecule_index);
    cell->remove_molecule(molecule_index,molecule_index_in_cell);

    // Move the last molecule into that memory location
    int last_molecule_index = num_molecules_local-1;

    if(last_molecule_index==molecule_index) {
        // We are lucky, the molecule to remove is the last one.
        num_molecules_local--;
        return;
    }

    // Remove the last molecule from current cell to re-add it after we moved all data
    cell = cell_containing_molecule(last_molecule_index);
    cell->remove_molecule(last_molecule_index,molecule_index_in_cell);

    // Copy the last molecule into the index of the removed molecule
    memcpy(&r.at(3*molecule_index),&r.at(3*last_molecule_index),3*sizeof(data_type));
    memcpy(&v.at(3*molecule_index),&v.at(3*last_molecule_index),3*sizeof(data_type));

    // Re-add it
    cell->add_molecule(molecule_index,molecule_index_in_cell, molecule_cell_index);

    num_molecules_local--;
}

bool System::validate_number_of_cells() {
    double num_voxels_per_cell_x = double(world_grid->global_nx) / cells_x;
    double num_voxels_per_cell_y = double(world_grid->global_ny) / cells_y;
    double num_voxels_per_cell_z = double(world_grid->global_nz) / cells_z;

    if(fmod(num_voxels_per_cell_x, 1.0) > 0) return false;
    if(fmod(num_voxels_per_cell_y, 1.0) > 0) return false;
    if(fmod(num_voxels_per_cell_z, 1.0) > 0) return false;

    double num_cells_per_cpu_x = double(cells_x)/topology->num_processors_vector[0];
    double num_cells_per_cpu_y = double(cells_y)/topology->num_processors_vector[1];
    double num_cells_per_cpu_z = double(cells_z)/topology->num_processors_vector[2];

    if(fmod(num_cells_per_cpu_x, 1.0) > 0) return false;
    if(fmod(num_cells_per_cpu_y, 1.0) > 0) return false;
    if(fmod(num_cells_per_cpu_z, 1.0) > 0) return false;

    return true;
}

void System::initialize(Settings *settings_, int myid_) {
    myid = myid_;
    settings = settings_;
    MAX_MOLECULE_NUM = settings->max_molecules_per_node;
    timer = new DSMCTimer();
    timer->start_system_initialize();
    steps = 0;
    collisions = 0;
    t = 0;

    init_randoms();
    unit_converter = new UnitConverter();
    length[0] = settings->Lx;
    length[1] = settings->Ly;
    length[2] = settings->Lz;
    one_over_length[0] = 1.0/length[0];
    one_over_length[1] = 1.0/length[1];
    one_over_length[2] = 1.0/length[2];

    temperature      = unit_converter->temperature_from_SI(settings->temperature);;
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);
    most_probable_velocity = sqrt(temperature);  // Most probable initial velocity
    density = unit_converter->number_density_from_SI(settings->density);
    diam = settings->diam;
    dt = settings->dt;
    atoms_per_molecule = settings->atoms_per_molecule;

    if(myid==0) cout << "Initializing system..." << endl;
    cell_length_x = length[0]/(settings->cells_x);
    cell_length_y = length[1]/(settings->cells_y);
    cell_length_z = length[2]/(settings->cells_z);

    cells_x = settings->cells_x;
    cells_y = settings->cells_y;
    cells_z = settings->cells_z;

    num_cells_vector[0] = cells_x;
    num_cells_vector[1] = cells_y;
    num_cells_vector[2] = cells_z;
    int num_cells_total = cells_x*cells_y*cells_z;
    cell_index_map.resize(num_cells_total,-1);

    io = new DSMC_IO(this);
    topology = new Topology(myid, settings->nx, settings->ny, settings->nz, this);
    node_num_new_molecules.resize(topology->num_processors, 0);
    node_molecule_data.resize(topology->num_processors,vector<data_type>(MAX_MPI_DATA,0));

//    node_molecule_data = new data_type*[topology->num_processors];
//    for(int i=0; i<topology->num_processors; i++) { node_molecule_data[i] = new data_type[MAX_MPI_DATA]; }
    MPI_Barrier(MPI_COMM_WORLD);
    if(myid==0) cout << "Loading world..." << endl;
    world_grid = new Grid(settings->ini_file.getstring("world"),this);
    MPI_Barrier(MPI_COMM_WORLD);

    if(myid == 0) {
        cout << "Validating voxel-cell-cpu configuration..." << endl;
        if(validate_number_of_cells()) {
            cout << "Nubmer of cells OK" << endl;
        } else {
            cout << "Error: voxel-cell-cpu configuration not good." << endl;
            cout << "Num voxels x: " << world_grid->global_nx << endl;
            cout << "Num voxels y: " << world_grid->global_ny << endl;
            cout << "Num voxels z: " << world_grid->global_nz << endl;
            cout << "Num cells x: " << cells_x << endl;
            cout << "Num cells y: " << cells_y << endl;
            cout << "Num cells z: " << cells_z << endl;
            cout << "Num processors x: " << topology->num_processors_vector[0] << endl;
            cout << "Num processors y: " << topology->num_processors_vector[1] << endl;
            cout << "Num processors z: " << topology->num_processors_vector[2] << endl;
            exit(1);
        }
    }

    // Calculate cell volume 
    volume = length[0]*length[1]*length[2];
    volume_global = volume*porosity_global;

    // First create all the cells
    if(myid==0) cout << "Creating cells..." << endl;
    setup_cells();

    // Update system volume with the correct porosity
    double volume_per_cpu = volume/topology->num_processors;
    volume = volume_per_cpu*porosity;
    num_molecules_local = density*volume/atoms_per_molecule;

    if(myid==0) {
        int num_molecules_global_calculated = density*volume_global/atoms_per_molecule;
        int num_molecules_per_cpu_calculated = num_molecules_global_calculated/topology->num_processors;

        cout << "Creating/loading " << num_molecules_global_calculated << " molecules (" << num_molecules_per_cpu_calculated << " per cpu)" << endl;

    }
    MPI_Barrier(MPI_COMM_WORLD);
    setup_molecules();
    MPI_Barrier(MPI_COMM_WORLD);
    num_molecules_global = 0;
    timer->start_mpi_reduce();
    MPI_Reduce(&num_molecules_local, &num_molecules_global, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD) ;
    timer->end_mpi_reduce();

    mpi_receive_buffer.resize(MAX_MPI_DATA,0);

    if(myid==0) cout << "Creating surface collider..." << endl;
    double sqrt_wall_temp_over_mass = sqrt(wall_temperature/settings->mass);
    ColliderBase *surface_collider;
    if(settings->surface_interaction_model.compare("thermal") == 0) {
        surface_collider = new ColliderThermal(sqrt_wall_temp_over_mass);
    } else if(settings->surface_interaction_model.compare("specular") == 0) {
        surface_collider = new ColliderSpecular();
    } else if(settings->surface_interaction_model.compare("cercignani_lampis") == 0) {
        surface_collider = new ColliderCercignaniLampis(sqrt_wall_temp_over_mass, this);
    } else if(settings->surface_interaction_model.compare("maxwell") == 0) {
        surface_collider = new ColliderMaxwell(sqrt_wall_temp_over_mass);
    }

    if(myid==0) cout << "Creating molecule mover..." << endl;
    mover = new MoleculeMover();
    mover->initialize(this, surface_collider);

    if(myid==0) {
        double mean_free_path = volume_global/(sqrt(2.0)*M_PI*diam*diam*num_molecules_global*atoms_per_molecule);
        int num_active_cells = num_cells_total* porosity_global;
        printf("done.\n\n");
        printf("%ld molecules (%ld per node)\n",num_molecules_global, num_molecules_global / topology->num_processors);
        printf("%d cells\n",num_cells_total);
        printf("Porosity: %f\n", porosity_global);
        printf("System volume: %f\n",length[0]*length[1]*length[2]);
        printf("Effective system volume: %f\n",volume_global);
        printf("Density: %E\n",unit_converter->number_density_to_SI(density));
        printf("Surface interaction model: %s\n",settings->surface_interaction_model.c_str());
        printf("Mean free path: %.4f \n",mean_free_path);
        printf("Mean free paths per cell: %.2f \n",min( min(length[0]/cells_x/mean_free_path,length[1]/cells_y/mean_free_path), length[2]/cells_z/mean_free_path));
        printf("%ld atoms per molecule\n",(unsigned long)atoms_per_molecule);
        printf("%ld molecules per active cell\n",num_molecules_global/num_active_cells);

        printf("dt = %f\n\n",dt);
        cout << endl;
    }

    timer->end_system_initialize();
}

inline void System::find_position(const int &index) {
    bool did_collide = true;
    bool is_inside = false;

    while(did_collide || !is_inside) {
        r.at(3*index + 0) = topology->origin[0] + topology->length[0]*rnd->next_double();
        r.at(3*index + 1) = topology->origin[1] + topology->length[1]*rnd->next_double();
        r.at(3*index + 2) = topology->origin[2] + topology->length[2]*rnd->next_double();

        did_collide = *world_grid->get_voxel(r,index)>=voxel_type_wall;
        is_inside = topology->is_position_inside(r, index);
    }
}

inline int System::cell_index_from_ijk(const int &i, const int &j, const int &k) {
    return i*cells_y*cells_z + j*cells_z + k;
}

int System::cell_index_from_position(const int &index) {
    int i = r.at(3*index + 0)/length[0]*cells_x;
    int j = r.at(3*index + 1)/length[1]*cells_y;
    int k = r.at(3*index + 2)/length[2]*cells_z;

    return cell_index_from_ijk(i,j,k);
}

void System::setup_molecules() {
    if(num_molecules_local > MAX_MOLECULE_NUM) {
        char error[1000];
        sprintf(error, "Node %d wants to create %ld molecules (MAX_MOLECULE_NUM=%d)",myid, num_molecules_local, MAX_MOLECULE_NUM);
        cout << error << endl;
        exit(1);
    }
    r.resize(3*MAX_MOLECULE_NUM);
    v.resize(3*MAX_MOLECULE_NUM);

    molecule_index_in_cell.resize(MAX_MOLECULE_NUM,-1);
    molecule_cell_index.resize(MAX_MOLECULE_NUM,-1);
    steps_since_collision.resize(MAX_MOLECULE_NUM,0);

    if(settings->load_previous_state) {
        io->load_state_from_file_binary();
        return;
    }

    double sqrt_temp_over_mass = sqrt(temperature/settings->mass);
    for(int n=0; n<num_molecules_local; n++) {
        v.at(3*n + 0) = rnd->next_gauss()*sqrt_temp_over_mass;
        v.at(3*n + 1) = rnd->next_gauss()*sqrt_temp_over_mass;
        v.at(3*n + 2) = rnd->next_gauss()*sqrt_temp_over_mass;
        find_position(n);
        int cell_index = cell_index_map[cell_index_from_position(n)];
        Cell *cell = active_cells.at(cell_index);
        cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
    }
}

void System::setup_cells() {
    int global_voxel_origin_x = topology->index_vector[0]*world_grid->nx_per_cpu;
    int global_voxel_origin_y = topology->index_vector[1]*world_grid->ny_per_cpu;
    int global_voxel_origin_z = topology->index_vector[2]*world_grid->nz_per_cpu;
    vector<Cell *> temp_cell_vector;
    for(int k=0;k<world_grid->nz_per_cpu;k++) {
        int c_z = (float)(k + global_voxel_origin_z)/world_grid->global_nz*cells_z;
        for(int j=0;j<world_grid->ny_per_cpu;j++) {
            int c_y = (float)(j + global_voxel_origin_y)/world_grid->global_ny*cells_y;
            for(int i=0;i<world_grid->nx_per_cpu;i++) {
                int c_x = (float)(i + global_voxel_origin_x)/world_grid->global_nx*cells_x;
                int cell_index = cell_index_from_ijk(c_x,c_y,c_z);
                Cell *cell;// = all_cells[cell_index];
                if(cell_index_map[cell_index] == -1) {
                    cell_index_map[cell_index] = temp_cell_vector.size(); // The current size is the index of this cell

                    cell = new Cell(this);
                    cell->index = cell_index;

                    cell->vr_max = 3*most_probable_velocity;
                    cell->Lx = cell_length_x; cell->Ly = cell_length_y; cell->Lz = cell_length_z;
                    cell->origin[0] = c_x*cell_length_x; cell->origin[1] = c_y*cell_length_y; cell->origin[2] = c_z*cell_length_z;
                    cell->index_vector[0] = c_x; cell->index_vector[1] = c_y; cell->index_vector[2] = c_z;
                    temp_cell_vector.push_back(cell);
                    cell->is_reservoir = false;
                } else {
                    int mapped_cell_index = cell_index_map[cell_index];
                    cell = temp_cell_vector.at(mapped_cell_index);
                }

                cell->total_pixels++;
                CVector voxel_index_adjusted_for_origin = CVector(i,j,k) + world_grid->voxel_origin;
                cell->pixels += *world_grid->get_voxel(voxel_index_adjusted_for_origin)<voxel_type_wall;
            }
        }
    }

    // Reset cell index map
    int num_cells_total = cells_x*cells_y*cells_z;
    for(int i=0; i<num_cells_total; i++) cell_index_map[i] = -1;

    double temp_volume = 0;
    // Update cell index map and put cells into active_cells instead. We only keep cells with positive volume.
    for(int i=0;i<temp_cell_vector.size();i++) {
        Cell *cell = temp_cell_vector.at(i);
        cell->vr_max = 3*most_probable_velocity;
        cell->update_volume();
        if(cell->volume > 0) {
            temp_volume += cell->volume;
            int cell_index = cell->index;
            cell_index_map[cell_index] = active_cells.size();
            active_cells.push_back(cell);
        } else {
            delete cell;
        }
    }

    temp_cell_vector.clear();
}

void System::init_randoms() {
    long seed = settings->seed;
    if(seed == 0)  seed = time(NULL);
    seed = -abs(seed + 1242461*(myid+1));
    rnd = new Random(seed, settings->alpha_n, settings->alpha_t);
}

void System::count_reservoir_particles() {
    return;
//    reservoir_b_particle_count = 0;

//    for(int i=0;i<reservoir_B_cells.size();i++) {
//        Cell *cell = reservoir_B_cells[i];
//        reservoir_b_particle_count += cell->num_molecules;
//    }
}

