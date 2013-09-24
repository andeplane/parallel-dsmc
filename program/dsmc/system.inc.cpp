#include <system.h>
#include <mpi.h>
#include <dsmctimer.h>
#include <moleculemover.h>
#include <settings.h>
#include <colliderbase.h>
#include <colliderspecular.h>
#include <colliderthermal.h>
#include <collidercercignanilampis.h>
#include <collidermaxwell.h>

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
    length[0] = settings->Lx;
    length[1] = settings->Ly;
    length[2] = settings->Lz;

    reservoir_size = length[settings->gravity_direction]*settings->reservoir_fraction/2;

    temperature      = unit_converter->temperature_from_SI(settings->temperature);;
    wall_temperature = unit_converter->temperature_from_SI(settings->wall_temperature);
    mpv = sqrt(temperature);  // Most probable initial velocity
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
    io = new DSMC_IO(this);

    cout << "Loading world..." << endl;

    world_grid = new Grid(settings->ini_file.getstring("world"),this);

    // First create all the cells
    cout << "Creating cells..." << endl;
    setup_cells();
    // Calculate porosity based on the world grid
    cout << "Calculating porosity..." << endl;
    calculate_porosity();
    // Calculate cell volume 
    volume = length[0]*length[1]*length[2];
    cout << "Updating cell volume..." << endl;
    update_cell_volume();
    // Update system volume with the correct porosity
    volume = length[0]*length[1]*length[2]*porosity;
    num_molecules = density*volume/atoms_per_molecule;

    cout << "Creating/loading molecules..." << endl;
    setup_molecules();

    mpi_receive_buffer = new double[9*MAX_MOLECULE_NUM];

    mfp = volume/(sqrt(2.0)*M_PI*diam*diam*num_molecules*atoms_per_molecule);
    
    cout << "Creating surface collider..." << endl;
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

    cout << "Creating molecule mover..." << endl;
    mover = new MoleculeMover();
    mover->initialize(this, surface_collider);

    int number_of_cells = all_cells.size();
    int number_of_cells_all = cells_x*cells_y*cells_z;

    printf("done.\n\n");
    printf("%d molecules\n",num_molecules);
    printf("%d (%d inactive) cells\n",number_of_cells,number_of_cells_all - number_of_cells);
    printf("Porosity: %f\n",porosity);
    printf("System volume: %f\n",length[0]*length[1]*length[2]);
    printf("Effective system volume: %f\n",volume);
    printf("Surface interaction model: %s\n",settings->surface_interaction_model.c_str());
    printf("Mean free path: %.4f \n",mfp);
    printf("Mean free paths per cell: %.2f \n",min( min(length[0]/cells_x/mfp,length[1]/cells_y/mfp), length[2]/cells_z/mfp));
    printf("%ld atoms per molecule\n",(unsigned long)atoms_per_molecule);
    printf("%d molecules per active cell\n",num_molecules/number_of_cells);

    printf("dt = %f\n\n",dt);
    cout << endl;
    timer->end_system_initialize();
}

inline void System::find_position(double *r) {
    bool did_collide = true;
    while(did_collide) {
        r[0] = length[0]*rnd->next_double();
        r[1] = length[1]*rnd->next_double();
        r[2] = length[2]*rnd->next_double();

        // did_collide = *world_grid->get_voxel(r)>=voxel_type_wall;
        double cylinder_center_x = length[0]*0.5;
        double cylinder_center_y = length[1]*0.5;
        double dx = r[0] - cylinder_center_x;
        double dy = r[1] - cylinder_center_y;
        double dr2 = dx*dx + dy*dy;
        did_collide = dr2 >= CYLINDER_RADIUS_SQUARED*0.9;
    }
}

inline int System::cell_index_from_ijk(const int &i, const int &j, const int &k) {
    return i*cells_y*cells_z + j*cells_z + k;
}

int System::cell_index_from_position(double *r) {
    int i = r[0]/length[0]*cells_x;
    int j = r[1]/length[1]*cells_y;
    int k = r[2]/length[2]*cells_z;

    return cell_index_from_ijk(i,j,k);
}

void System::update_cell_volume() {
    active_cells.reserve(all_cells.size());

    for(int i=0;i<all_cells.size();i++) {
        Cell *cell = all_cells[i];
        cell->vr_max = 3*mpv;
        cell->update_volume();
        if(cell->volume>0) active_cells.push_back(cell);
    }
}

void System::setup_molecules() {
    r = new double[3*MAX_MOLECULE_NUM];

    molecule_index_in_cell = new unsigned long[MAX_MOLECULE_NUM];
    molecule_cell_index    = new unsigned long[MAX_MOLECULE_NUM];

    v = new double[3*MAX_MOLECULE_NUM];
    r0 = new double[3*MAX_MOLECULE_NUM];

    if(settings->load_previous_state) {
        io->load_state_from_file_binary();
        return;
    }

    double sqrt_temp_over_mass = sqrt(temperature/settings->mass);

    for(int n=0;n<num_molecules;n++) {
        v[3*n+0] = rnd->next_gauss()*sqrt_temp_over_mass;
        v[3*n+1] = rnd->next_gauss()*sqrt_temp_over_mass;
        v[3*n+2] = rnd->next_gauss()*sqrt_temp_over_mass;
        find_position(&r[3*n]);
        r0[3*n+0] = r[3*n+0];
        r0[3*n+1] = r[3*n+1];
        r0[3*n+2] = r[3*n+2];
        Cell *cell = all_cells[cell_index_from_position(&r[3*n])];
        cell->add_molecule(n,this->molecule_index_in_cell,this->molecule_cell_index);
    }
}

void System::setup_cells() {
    int num_cells = cells_x*cells_y*cells_z;
    all_cells.reserve(num_cells);
    int idx[3];
    int count = 0;

    for(idx[0]=0;idx[0]<cells_x;idx[0]++) {
        for(idx[1]=0;idx[1]<cells_y;idx[1]++) {
            for(idx[2]=0;idx[2]<cells_z;idx[2]++) {
                Cell *cell = new Cell(this);

                cell->index = cell_index_from_ijk(idx[0],idx[1],idx[2]);
                cell->vr_max = 3*mpv;
                all_cells.push_back(cell);
                cell->is_reservoir = false;

                if( (double)idx[settings->gravity_direction]/num_cells_vector[settings->gravity_direction] < settings->reservoir_fraction/2) {
                    reservoir_A_cells.push_back(cell);
                    cell->is_reservoir = true;
                } else if( (double)idx[settings->gravity_direction]/num_cells_vector[settings->gravity_direction] > (1-settings->reservoir_fraction/2)) {
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

    int i_start = 0; //float(my_vector_index[0])*settings->cells_per_node_x/cells_x*world_grid->Nx;
    int i_end   = world_grid->Nx; //float(my_vector_index[0]+1)*settings->cells_per_node_x/cells_x*world_grid->Nx;

    int j_start = 0; //float(my_vector_index[1])*settings->cells_per_node_y/cells_y*world_grid->Ny;
    int j_end   = world_grid->Ny; //float(my_vector_index[1]+1)*settings->cells_per_node_y/cells_y*world_grid->Ny;

    int k_start = 0; //float(my_vector_index[2])*settings->cells_per_node_z/cells_z*world_grid->Nz;
    int k_end   = world_grid->Nz; //float(my_vector_index[2]+1)*settings->cells_per_node_z/cells_z*world_grid->Nz;
    int cell_index, c_x, c_y, c_z;

    for(int k=k_start;k<k_end;k++) {
        for(int j=j_start;j<j_end;j++) {
            for(int i=i_start;i<i_end;i++) {
                c_x = i*cells_x/(float)world_grid->Nx;
                c_y = j*cells_y/(float)world_grid->Ny;
                c_z = k*cells_z/(float)world_grid->Nz;
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
    seed = 2;
    rnd = new Random(-seed, settings->alpha_n, settings->alpha_t);
}

void System::count_reservoir_particles() {
    reservoir_b_particle_count = 0;

    for(int i=0;i<reservoir_B_cells.size();i++) {
        Cell *cell = reservoir_B_cells[i];
        reservoir_b_particle_count += cell->num_molecules;
    }
}
