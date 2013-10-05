#include <statisticssampler.h>

#include <cell.h>
#include <math.h>
#include <unitconverter.h>
#include <dsmc_io.h>
#include <system.h>
#include <mpi.h>
#include <settings.h>
#include <moleculemover.h>
#include <colliderbase.h>

int num_bins_per_dimension;
int num_bins;
double *vel;
int *count_;

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;
    settings = system->settings;
    temperature_sampled_at = -1;
    kinetic_energy_sampled_at = -1;
    velocity_distribution_sampled_at = -1;
    flux_sampled_at = -1;
    permeability_sampled_at = -1;
    density_sampled_at = -1;
    linear_density_sampled_at = -1;
    permeability = 0;
    flux = 0;
    num_samples = 0;

    num_bins_per_dimension = settings->velocity_bins;
    num_bins = num_bins_per_dimension*num_bins_per_dimension;
    vel = new double[num_bins];
    count_ = new int[num_bins];

    memset((void*)count_,0,num_bins*sizeof(int));
    memset((void*)vel,0,num_bins*sizeof(double));
}

void StatisticsSampler::sample() {
    if(settings->statistics_interval && system->steps % settings->statistics_interval != 0) return;
    double t_in_nano_seconds = system->unit_converter->time_to_SI(system->t)*1e9;

//    if(settings->velocity_profile_type.compare("other") == 0) {
//        sample_velocity_distribution();
//    } else if(settings->velocity_profile_type.compare("cylinder") == 0) {
//        sample_velocity_distribution_cylinder();
//    } else if(settings->velocity_profile_type.compare("box") == 0) {
//        sample_velocity_distribution_box();
//    }
    sample_temperature();
//    sample_density();
//    sample_linear_density();
//    sample_permeability();

    double kinetic_energy_per_molecule = kinetic_energy / (system->num_molecules_global*system->atoms_per_molecule);

    if(system->myid==0) {

        fprintf(system->io->energy_file, "%f %f %f\n",t_in_nano_seconds, system->unit_converter->energy_to_eV(kinetic_energy_per_molecule), system->unit_converter->temperature_to_SI(temperature));

//        if(settings->maintain_pressure) {
//            fprintf(system->io->flux_file, "%f %ld\n",t_in_nano_seconds, system->flux_count);
//        } else {
//            fprintf(system->io->flux_file, "%f %ld\n",t_in_nano_seconds, system->mover->count_periodic[settings->flow_direction]);
//        }

//        fprintf(system->io->permeability_file, "%f %E\n",t_in_nano_seconds, system->unit_converter->permeability_to_SI(permeability));

        double pressure = system->num_molecules_global*system->atoms_per_molecule / system->volume * temperature;
        cout << system->steps << "   t=" << t_in_nano_seconds << "   T=" << system->unit_converter->temperature_to_SI(temperature) << "   Collisions: " <<  system->collisions <<   "   Wall collisions: " << system->mover->surface_collider->num_collisions << "   Pressure: " << system->unit_converter->pressure_to_SI(pressure) <<  "   Molecules: " << system->num_molecules_global << endl ;
        fprintf(system->io->pressure_file, "%f %E\n",t_in_nano_seconds, pressure);
        fprintf(system->io->num_molecules_file, "%f %ld\n",t_in_nano_seconds, system->num_molecules_global);
    }

    num_samples++;
}

void StatisticsSampler::sample_kinetic_energy() {
    if(system->steps == kinetic_energy_sampled_at) return;
    kinetic_energy_sampled_at = system->steps;

    double kinetic_energy_local = 0;

    for(unsigned int i=0;i<system->num_molecules_local;i++) {
        kinetic_energy_local += (system->v[3*i+0]*system->v[3*i+0] + system->v[3*i+1]*system->v[3*i+1] + system->v[3*i+2]*system->v[3*i+2]);
    }
    kinetic_energy_local *= 0.5*settings->mass*system->atoms_per_molecule;
    kinetic_energy = 0;
    MPI_Reduce(&kinetic_energy_local, &kinetic_energy,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

}

void StatisticsSampler::sample_temperature() {
    if(system->steps == temperature_sampled_at) return;
    temperature_sampled_at = system->steps;

    sample_kinetic_energy();
    double kinetic_energy_per_molecule = kinetic_energy / (system->num_molecules_global*system->atoms_per_molecule);
    temperature = 2.0/3*kinetic_energy_per_molecule;
}

void StatisticsSampler::sample_flux() {
    if(system->steps == flux_sampled_at) return;
    flux_sampled_at = system->steps;
    // Either, the gas is pressure driven or gravity driven. We measure flux different in these two cases.
    if(settings->maintain_pressure) {
        flux = system->flux_count / system->t;
    } else {
        flux = system->mover->count_periodic[settings->flow_direction] / system->t;
    }

}

void StatisticsSampler::sample_permeability() {
    if(system->steps == permeability_sampled_at) return;
    permeability_sampled_at = system->steps;
    if(settings->flow_direction<0) return;

    sample_flux();
    double volume_per_molecule = system->volume / system->num_molecules_local;
    double viscosity_dsmc_units = system->unit_converter->viscosity_from_SI(settings->viscosity);
    double volume_flux = flux*volume_per_molecule;
    double L = system->length[settings->flow_direction];
    double mass_density = system->density*settings->mass;

    double area = 1;
    for(int a=0;a<3;a++) {
        if(a != settings->flow_direction) area *= system->length[a];
    }
    double pressure_in_reservoir_a = system->unit_converter->pressure_from_SI(settings->pressure_A);
    double pressure_in_reservoir_b = system->unit_converter->pressure_from_SI(settings->pressure_B);

    if(settings->maintain_pressure) {
        // Expression from Darcy's law of gases
        permeability = 2*pressure_in_reservoir_b*volume_flux*L*viscosity_dsmc_units / (area * (pressure_in_reservoir_a*pressure_in_reservoir_a - pressure_in_reservoir_b*pressure_in_reservoir_b));
    } else {
        permeability = volume_flux*L*viscosity_dsmc_units / (mass_density*system->length[settings->flow_direction]*settings->gravity*area);
    }
}

void StatisticsSampler::sample_velocity_distribution_cylinder() {
    int N = this->system->settings->velocity_bins;
    double center_x = system->length[0]/2;
    double center_y = system->length[1]/2;

    double *v_of_r = new double[N];
    int *v_of_r_count = new int[N];
    memset(v_of_r,0,N*sizeof(double));
    memset(v_of_r_count,0,N*sizeof(int));

    double dr_max = sqrt(system->length[0]*system->length[0]+system->length[1]*system->length[1]);

    for(int i=0;i<system->num_molecules_local;i++) {
        double dx = system->r[3*i+0] - center_x;
        double dy = system->r[3*i+1] - center_y;
        double dr = sqrt(dx*dx + dy*dy);
        int v_of_r_index = N*dr/dr_max;
        if(v_of_r_index>=N) continue;
        
        // double v_norm = sqrt(system->v[3*i+2]*system->v[3*i+2] + system->v[3*i+1]*system->v[3*i+1] + system->v[3*i+0]*system->v[3*i+0]);
        double vz = system->v[3*i+2];

        v_of_r[v_of_r_index] += vz;
        v_of_r_count[v_of_r_index]++;
    }

    for(int i=0;i<N;i++) {
        if(v_of_r_count[i]>0) v_of_r[i] /= v_of_r_count[i];
        fprintf(system->io->velocity_file,"%f ",system->unit_converter->velocity_to_SI(v_of_r[i]));
        // fprintf(system->io->velocity_file,"%d ", v_of_r_count[i]);
    }
    fprintf(system->io->velocity_file,"\n");
    delete v_of_r;
    delete v_of_r_count;

}

void StatisticsSampler::sample_velocity_distribution_box() {
    int N = this->system->settings->velocity_bins;

    double *v_of_y = new double[N];
    int *v_of_y_count = new int[N];
    memset(v_of_y,0,N*sizeof(double));
    memset(v_of_y_count,0,N*sizeof(int));

    for(int i=0;i<system->num_molecules_local;i++) {
        double y = system->r[3*i+1];
        int v_of_y_index = N*(y/system->length[1]);

        // double v_norm = sqrt(system->v[3*i+2]*system->v[3*i+2] + system->v[3*i+1]*system->v[3*i+1] + system->v[3*i+0]*system->v[3*i+0]);
        double vz = system->v[3*i+2];

        if(v_of_y_index >= N) continue;
        v_of_y[v_of_y_index] += vz;
        v_of_y_count[v_of_y_index]++;
    }

    for(int i=0;i<N;i++) {
        if(v_of_y_count[i]>0) v_of_y[i] /= v_of_y_count[i];
        fprintf(system->io->velocity_file,"%f ",system->unit_converter->velocity_to_SI(v_of_y[i]));
    }
    fprintf(system->io->velocity_file,"\n");
    delete v_of_y;
    delete v_of_y_count;
}

void StatisticsSampler::sample_density() {
    if(system->steps == density_sampled_at) return;
    density_sampled_at = system->steps;

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_x = system->r[3*i+0] / system->length[0]*num_bins_per_dimension;
        int bin_y = system->r[3*i+1] / system->length[1]*num_bins_per_dimension;
        int bin_z = system->r[3*i+2] / system->length[2]*num_bins_per_dimension;

        int index = bin_x*num_bins_per_dimension + bin_y;

        count_[index]++;
    }
}

void StatisticsSampler::sample_linear_density() {
    this->sample_temperature();
    int num_bins = this->system->settings->velocity_bins;

    int *linear_density_count = new int[num_bins];
    memset(linear_density_count,0,num_bins*sizeof(int));

    for(int i=0;i<system->num_molecules_local;i++) {
        double z = system->r[3*i+2];
        int bin_index = num_bins*(z/system->length[2]);
        linear_density_count[bin_index]++;
    }

    double volume_per_bin = system->volume / num_bins;

    for(int i=0;i<num_bins;i++) {
        double density = system->atoms_per_molecule*linear_density_count[i]/volume_per_bin;
        double pressure = density*temperature; // Ideal gas law

        fprintf(system->io->linear_density_file,"%E ", system->unit_converter->number_density_to_SI(density));
        fprintf(system->io->linear_pressure_file,"%E ",system->unit_converter->pressure_to_SI(pressure));
    }

    fprintf(system->io->linear_density_file,"\n");
    fprintf(system->io->linear_pressure_file,"\n");
    delete linear_density_count;
}

void StatisticsSampler::sample_velocity_distribution() {
    if(system->steps == velocity_distribution_sampled_at) return;
    velocity_distribution_sampled_at = system->steps;
    density_sampled_at = system->steps;

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_x = system->r[3*i+0]*system->one_over_length[0]*num_bins_per_dimension;
        int bin_y = system->r[3*i+1]*system->one_over_length[1]*num_bins_per_dimension;
        int bin_z = system->r[3*i+2]*system->one_over_length[2]*num_bins_per_dimension;

        // int index = bin_x*num_bins_per_dimension + bin_y;
        int index = bin_x*num_bins_per_dimension + bin_y;

        vel[index] += system->v[3*i+2];
        count_[index]++;
    }
}

void StatisticsSampler::finalize() {
    for(int i=0;i<num_bins;i++) {
        if(settings->velocity_profile_type.compare("other") == 0) {
            // If we sample the 2d field, the averages are done in cpp code
            if(count_[i]>0) vel[i] /= count_[i];
            fprintf(system->io->velocity_file,"%f ",system->unit_converter->velocity_to_SI(vel[i]));
        }

        fprintf(system->io->density_file,"%f ",(double)count_[i] / this->num_samples);
    }

    fprintf(system->io->velocity_file,"\n");
    fclose(system->io->energy_file);
    fclose(system->io->velocity_file);
    fclose(system->io->flux_file);
    fclose(system->io->permeability_file);
}

/*
void StatisticsSampler::calculate_velocity_profile() {
    if(!print_velocity_profile) return;
    if(++velocity_profile_samples % sample_every_n) return;

    int N = 100;
	Molecule *molecule;
    vec velocities = zeros<vec>(N,1);
    vec velocity_count= zeros<vec>(N,1);

	int n;


    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        n = N*molecule->r(1)/system->height;
        velocities(n) += molecule->v(0);
        velocity_count(n)++;
	}

    for(n=0;n<N;n++) {
        if(velocity_count(n)>0)
            velocities(n) /= velocity_count(n);
    }
	
	for(int n=0;n<N;n++) 
        fprintf(velocity_file,"%f ",velocities(n));

    fprintf(velocity_file,"\n");
}

double StatisticsSampler::get_temperature() {
    if(!print_temperature) return 0;

    return temperature_sum/(temperature_samples/ini->getint("sample_every_n"));
}

void StatisticsSampler::calculate_velocity_field() {
    if(!print_velocity_field) return;

    Cell *c;
    for(int j=0;j<system->settings->cells_y;j++) {
        for(int i=0;i<system->settings->cells_x;i++) {
            c = system->cells[i][j];

            fprintf(velocity_field_file_x,"%.5f ",c->momentum(0));
            fprintf(velocity_field_file_y,"%.5f ",c->momentum(1));
        }

        fprintf(velocity_field_file_x,"\n");
        fprintf(velocity_field_file_y,"\n");
    }
}

double StatisticsSampler::calculate_diffusion_constant() {
    Molecule *molecule;
    double r_squared = 0;
    for(int i=0;i<system->N;i++) {
        molecule = system->molecules[i];
        r_squared += molecule->squared_distance_from_initial_position();
    }
    r_squared /= system->N;

    double diffusion_constant = r_squared/6/system->t;
    return diffusion_constant;
}


mat StatisticsSampler::calculate_global_pressure_tensor() {
    vector<Molecule*> molecules;
    for(int n=0;n<system->N;n++)
        molecules.push_back(system->molecules[n]);

    return calculate_pressure_tensor(molecules);
}

vector<mat> StatisticsSampler::calculate_local_pressure_tensor() {
    vector<mat> pressure_tensors;

    Cell *c;
    for(int i=0;i<system->settings->cells_x;i++) {
        for(int j=0;j<system->settings->cells_y;j++) {
            c = system->cells[i][j];
            vector<Molecule*> molecules;
            for(int n=0;n<c->particles;n++)
                molecules.push_back(system->molecules[c->particle_indices[n]]);
            pressure_tensors.push_back(calculate_pressure_tensor(molecules));
        }
    }

    return pressure_tensors;
}

mat StatisticsSampler::calculate_pressure_tensor(vector<Molecule*> molecules) {
    mat p = zeros<mat>(3,3);
    vec v = zeros<vec>(3,1);
    double rho = 0;
    int N = molecules.size();
    if(!N) return p;

    Molecule *m;

    for(int n=0;n<molecules.size();n++) {
        m = molecules[n];
        v += m->mass*m->v/N;
        rho += m->mass/N;
        for(int i=0;i<3;i++) {
            for(int j=0;j<3;j++) {
                p(i,j) += m->v(i)*m->v(j)*m->mass/N;
            }
        }
    }
    v /= rho;

    for(int i=0;i<3;i++) {
        for(int j=0;j<3;j++) {
            p(i,j) -= v(i)*v(j)*rho;
        }
    }

    return p;
}

void StatisticsSampler::update_cell_statistics() {
    Cell *c; Molecule *m;
    for(int i=0;i<system->settings->cells_x;i++) {
        for(int j=0;j<system->settings->cells_y;j++) {
            c = system->cells[i][j];
            vector<Molecule*> molecules;
            vec momentum = zeros<vec>(3,1);
            double energy = 0;

            for(int n=0;n<c->particles;n++) {
                m = system->molecules[c->particle_indices[n]];

                energy += 0.5*dot(m->v,m->v)*m->mass*m->atoms;
                momentum += m->v*m->mass*m->atoms;

                molecules.push_back(m);
            }
            c->update_pressure_tensor(calculate_pressure_tensor(molecules));
            c->update_momentum(momentum);
            c->update_energy(energy);
            c->update_temperature(2*energy/(3*c->particles));
        }
    }
}
*/
