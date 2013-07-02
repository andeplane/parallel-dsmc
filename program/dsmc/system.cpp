#include <system.h>

#include <math.h>
#include <fstream>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <unitconverter.h>
#include <dsmc_io.h>
#include <time.h>
#include <dsmctimer.h>
#include <mpi.h>
#include <dsmctimer.h>
#include <moleculemover.h>
#include <settings.h>

void System::set_topology() {
    /*----------------------------------------------------------------------
    Defines a logical network topology.  Prepares a neighbor-node ID table,
    nn, & a shift-vector table, sv, for internode message passing.  Also
    prepares the node parity table, myparity.
    ----------------------------------------------------------------------*/
    length[0] = settings->Lx;
    length[1] = settings->Ly;
    length[2] = settings->Lz;
    num_processors[0] = settings->nodes_x;
    num_processors[1] = settings->nodes_y;
    num_processors[2] = settings->nodes_z;
    num_nodes = num_processors[0]*num_processors[1]*num_processors[2];
    node_index[0] = myid/(num_processors[1]*num_processors[2]);
    node_index[1] = (myid/num_processors[2]) % num_processors[1];
    node_index[2] = myid%num_processors[2];
    num_cells_per_node[0] = settings->cells_per_node_x;
    num_cells_per_node[1] = settings->cells_per_node_y;
    num_cells_per_node[2] = settings->cells_per_node_z;

    for(int a=0; a<3; a++) node_length[a] = length[a] / num_processors[0];
    for(int a=0; a<3; a++) num_cells[a] = num_cells_per_node[a]*num_processors[a];
    for(int a=0; a<3; a++) cell_length[a] = length[a]/(num_cells[a]);
    for(int a=0; a<3; a++) origo[a] = (double)node_index[a] * node_length[a];

    num_cells_total = num_cells[0]*num_cells[1]*num_cells[2];

    /* Integer vectors to specify the six neighbor nodes */
    int integer_vector[6][3] = {
        {-1,0,0}, {1,0,0}, {0,-1,0}, {0,1,0}, {0,0,-1}, {0,0,1}
    };

    int k1[3];

    for (int n=0; n<6; n++) {
        for (int a=0; a<3; a++) {
            k1[a] = (node_index[a]+integer_vector[n][a]+num_processors[a])%num_processors[a];
        }

        /* Scalar neighbor ID, nn */
        neighbor_nodes[n] = k1[0]*num_processors[1]*num_processors[2]+k1[1]*num_processors[2]+k1[2];
        /* Shift vector, sv */
        for (int a=0; a<3; a++) shift_vector[n][a] = node_length[a]*integer_vector[n][a];
    }
}

void System::sync_mpi_initialize() {
    unsigned long num_active_cells_local = active_cells.size();

    MPI_Allreduce(&num_molecules_local,&num_molecules_global,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&porosity,&porosity_global,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&num_active_cells_local,&num_active_cells,1,MPI_UNSIGNED_LONG,MPI_SUM,MPI_COMM_WORLD);
    porosity_global /= num_nodes;
}

void System::initialize(Settings *settings_, int myid_) {
    myid = myid_;
    settings = settings_;
    timer = new DSMCTimer();
    timer->start_system_initialize();
    steps = 0;
    collisions = 0;
    t = 0;

    init_randoms();
    unit_converter = new UnitConverter();
    set_topology();

    reservoir_size = length[settings->gravity_direction]*settings->reservoir_fraction/2;
    temperature      = unit_converter->temperature_from_SI(settings->temperature);
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);

    density = unit_converter->number_density_from_SI(settings->density);
    diam = settings->diam;
    atoms_per_molecule = settings->atoms_per_molecule;
    if(myid==0) cout << "Initializing system..." << endl;
    io = new DSMC_IO(this);
    if(myid==0) cout << "Loading world..." << endl;

    world_grid = new Grid(settings->ini_file.getstring("world"),this);

    if(myid==0) cout << "Creating cells..." << endl;
    setup_cells();
    if(myid==0) cout << "Calculating porosity..." << endl;
    calculate_porosity();
    volume = length[0]*length[1]*length[2]*porosity;
    num_molecules_local = density*volume/atoms_per_molecule;

    if(myid==0) cout << "Creating/loading molecules..." << endl;
    setup_molecules();

    mpi_receive_buffer = new double[9*MAX_MOLECULE_NUM];
    mpi_send_buffer    = new double[9*MAX_MOLECULE_NUM];

    mpv = sqrt(temperature);  // Most probable initial velocity
    dt = settings->dt;

    if(myid==0) cout << "Updating cell volume..." << endl;
    update_cell_volume();

    if(myid==0) cout << "Creating molecule mover..." << endl;
    mover = new MoleculeMover();
    mover->initialize(this);

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*num_molecules_global*atoms_per_molecule);

    sync_mpi_initialize();
    if(myid==0) {
        printf("done.\n\n");
        printf("%ld molecules\n",num_molecules_global);
        printf("%ld (%ld inactive) cells\n",num_active_cells,num_cells_total - num_active_cells);
        printf("Porosity: %f\n",porosity);
        printf("System volume: %f\n",length[0]*length[1]*length[2]);
        printf("Effective system volume: %f\n",volume);
        printf("Mean free path: %.4f \n",mfp);
        printf("Mean free paths per cell: %.2f \n",min( min(length[0]/num_cells[0]/mfp,length[1]/num_cells[1]/mfp), length[2]/num_cells[2]/mfp));
        printf("%ld atoms per molecule\n",(unsigned long)atoms_per_molecule);
        printf("%d molecules per active cell\n",int(num_molecules_global/num_active_cells));

        printf("dt = %f\n\n",dt);
        cout << endl;
    }

    timer->end_system_initialize();
}
int collis_lols = 0;

inline void System::find_position(double *r) {
    bool did_collide = true;
    while(did_collide) {
        collis_lols++;
        r[0] = origo[0] + node_length[0]*rnd->nextDouble();
        r[1] = origo[1] + node_length[1]*rnd->nextDouble();
        r[2] = origo[2] + node_length[2]*rnd->nextDouble();

        did_collide = *world_grid->get_voxel(r)>=voxel_type_wall;
    }
}

inline int System::cell_index_from_ijk(const int &i, const int &j, const int &k) {
    return i*num_cells_per_node[1]*num_cells_per_node[2] + j*num_cells_per_node[2] + k;
}

int System::cell_index_from_position(double *r) {
    int i = (r[0] - origo[0])/node_length[0]*num_cells_per_node[0];
    int j = (r[1] - origo[1])/node_length[1]*num_cells_per_node[1];
    int k = (r[2] - origo[2])/node_length[2]*num_cells_per_node[2];

    return cell_index_from_ijk(i,j,k);
}

void System::update_cell_volume() {
    active_cells.reserve(all_cells.size());

    for(unsigned long i=0;i<all_cells.size();i++) {
        Cell *cell = all_cells[i];
        cell->vr_max = 3*mpv;
        cell->update_volume();
        if(cell->volume>0) active_cells.push_back(cell);
    }
}

void System::set_initial_positions() {
    for(unsigned long n=0; n<num_molecules_local; n++) {
        for(int a=0; a<3; a++) {
            r0[3*n+a] = r[3*n+a] + origo[a]; // Initial position
        }
    }
}

void System::setup_molecules() {
    r  = new double[3*MAX_MOLECULE_NUM];
    v  = new double[3*MAX_MOLECULE_NUM];
    r0 = new double[3*MAX_MOLECULE_NUM];

    molecule_index_in_cell = new unsigned long[MAX_MOLECULE_NUM];
    molecule_cell_index    = new unsigned long[MAX_MOLECULE_NUM];

    if(settings->load_previous_state) {
        io->load_state_from_file_binary();
        set_initial_positions();
        return;
    }

    double sqrt_temp_over_mass = sqrt(temperature/settings->mass);

    for(unsigned long n=0;n<num_molecules_local;n++) {
        v[3*n+0] = rnd->nextGauss()*sqrt_temp_over_mass;
        v[3*n+1] = rnd->nextGauss()*sqrt_temp_over_mass;
        v[3*n+2] = rnd->nextGauss()*sqrt_temp_over_mass;
        find_position(&r[3*n]);

        set_initial_positions();
        Cell *cell = all_cells[cell_index_from_position(&r[3*n])];
        cell->add_molecule(n,this->molecule_index_in_cell,this->molecule_cell_index);
    }
    io->save_state_to_file_binary();
}

void System::setup_cells() {
    all_cells.reserve(num_cells_total);
    int idx[3];

    for(idx[0]=0;idx[0]<num_cells_per_node[0];idx[0]++) {
        for(idx[1]=0;idx[1]<num_cells_per_node[1];idx[1]++) {
            for(idx[2]=0;idx[2]<num_cells_per_node[2];idx[2]++) {
                Cell *cell = new Cell(this);

                cell->index = cell_index_from_ijk(idx[0],idx[1],idx[2]);
                cell->vr_max = 3*mpv;
                all_cells.push_back(cell);
                cell->is_reservoir = false;

                if( (double)idx[settings->gravity_direction]/num_cells_per_node[settings->gravity_direction] < settings->reservoir_fraction/2) {
                    reservoir_A_cells.push_back(cell);
                    cell->is_reservoir = true;
                } else if( (double)idx[settings->gravity_direction]/num_cells_per_node[settings->gravity_direction] > (1-settings->reservoir_fraction/2)) {
                    reservoir_B_cells.push_back(cell);
                    cell->is_reservoir = true;
                }
            }
        }
    }
}

void System::calculate_porosity() {
    int filled_pixels = 0;
    int all_pixels = 0;

    int i_start = float(node_index[0])*num_cells_per_node[0]/num_cells[0]*world_grid->Nx;
    int i_end   = float(node_index[0]+1)*num_cells_per_node[0]/num_cells[0]*world_grid->Nx;

    int j_start = float(node_index[1])*num_cells_per_node[1]/num_cells[1]*world_grid->Ny;
    int j_end   = float(node_index[1]+1)*num_cells_per_node[1]/num_cells[1]*world_grid->Ny;

    int k_start = float(node_index[2])*num_cells_per_node[2]/num_cells[2]*world_grid->Nz;
    int k_end   = float(node_index[2]+1)*num_cells_per_node[2]/num_cells[2]*world_grid->Nz;

    int cell_index, c_x, c_y, c_z;

    for(int k=k_start;k<k_end;k++) {
        for(int j=j_start;j<j_end;j++) {
            for(int i=i_start;i<i_end;i++) {
                c_x = i*num_cells_per_node[0]/(float)world_grid->Nx;
                c_y = j*num_cells_per_node[1]/(float)world_grid->Ny;
                c_z = k*num_cells_per_node[2]/(float)world_grid->Nz;
                cell_index = cell_index_from_ijk(c_x,c_y,c_z);
                Cell *c = all_cells[cell_index];

                c->total_pixels++;
                all_pixels++;

                c->pixels += *world_grid->get_voxel(i,j,k)<voxel_type_wall;
                filled_pixels += *world_grid->get_voxel(i,j,k)<voxel_type_wall;
            }
        }
    }

    porosity = (float)filled_pixels / all_pixels;
}

void System::init_randoms() {
    long seed = time(NULL);
    seed = 1;
    rnd = new Random(-seed);
}

void System::step() {
    // steps += 1;
    // t += dt;
    // accelerate();
    // move();
    // collide();
    // if(settings->maintain_pressure) maintain_pressure();
}

void System::move() {
    timer->start_moving();
    for(unsigned long n=0;n<num_molecules_local;n++) {
        mover->move_molecule(n,dt,rnd,0);
    }

    timer->end_moving();

    timer->start_moving();
    if(myid==0) update_molecule_cells();
    timer->end_moving();
}

void System::collide() {
    timer->start_colliding();

    for(unsigned long i=0;i<active_cells.size();i++) {
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

    for(unsigned long n=0;n<num_molecules_local;n++) {
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
    int n = num_molecules_local;

    v[3*n+0] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    v[3*n+1] = rnd->nextGauss()*sqrt(temperature/settings->mass);
    v[3*n+2] = rnd->nextGauss()*sqrt(temperature/settings->mass);

    find_position_in_reservoirs(&r[3*n],add_in_A);
    r0[3*n+0] = r[3*n+0];
    r0[3*n+1] = r[3*n+1];
    r0[3*n+2] = r[3*n+2];

    Cell *cell = all_cells[cell_index_from_position(&r[3*n])];
    cell->add_molecule(n,molecule_index_in_cell,molecule_cell_index);
    num_molecules_local++;
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
        return true;
    } else return false;
}

void System::update_molecule_cells() {
    for(unsigned long n=0;n<num_molecules_local;n++) {
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

    for(unsigned long i=0;i<reservoir_A_cells.size();i++) {
        Cell *cell = reservoir_A_cells[i];
        num_molecules_in_reservoir += cell->num_molecules;
        volume_in_reservoir += cell->volume;
        kinetic_energy_in_reservoir += cell->calculate_kinetic_energy();
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

    for(unsigned long i=0;i<reservoir_B_cells.size();i++) {
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
