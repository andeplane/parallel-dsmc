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
#include <dsmctimer.h>
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

    num_bins = settings->sampling_bins;

    if(settings->velocity_profile_type.compare("area") == 0) {
        velocity_distribution.resize(num_bins*num_bins,0);
        velocity_distribution_count.resize(num_bins*num_bins,0);
    } else {
        velocity_distribution.resize(num_bins,0);
        velocity_distribution_count.resize(num_bins,0);
    }

    count_across_channel.resize(num_bins,0);
    kinetic_energy_across_channel.resize(num_bins,0);
}

void StatisticsSampler::sample() {
    if(settings->statistics_interval && system->steps % settings->statistics_interval != 0) return;
    double t_in_nano_seconds = system->unit_converter->time_to_SI(system->t)*1e9;
    sample_velocity_distribution();

    sample_temperature();
    sample_stats_across_channel();
    sample_permeability();

    double kinetic_energy_per_molecule = kinetic_energy / (system->num_molecules_global*system->atoms_per_molecule);
    collisions = 0;
    wall_collisions = 0;
    system->timer->start_mpi_reduce();
    MPI_Reduce(&system->collisions,&collisions,1,MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&system->mover->surface_collider->num_collisions,&wall_collisions,1,MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);
    system->timer->end_mpi_reduce();


    if(system->myid==0) {
        double pressure = system->num_molecules_global*system->atoms_per_molecule / system->volume_global * temperature;
        cout << system->steps << "   t=" << t_in_nano_seconds << "   T=" << system->unit_converter->temperature_to_SI(temperature) << "   Collisions: " <<  collisions <<   "   Wall collisions: " << wall_collisions << "   Pressure: " << system->unit_converter->pressure_to_SI(pressure) <<  "   Molecules: " << system->num_molecules_global << endl ;
        fprintf(system->io->energy_file, "%f %f %f\n",t_in_nano_seconds, system->unit_converter->energy_to_eV(kinetic_energy_per_molecule), system->unit_converter->temperature_to_SI(temperature));
        fprintf(system->io->pressure_file, "%f %E\n",t_in_nano_seconds, pressure);
        fprintf(system->io->flux_file, "%f %E\n",t_in_nano_seconds, flux);
        fprintf(system->io->num_molecules_file, "%f %ld\n",t_in_nano_seconds, system->num_molecules_global);
        fprintf(system->io->temperature_file, "%f %f\n",t_in_nano_seconds, system->unit_converter->temperature_to_SI(temperature));
        fprintf(system->io->permeability_file, "%f %E\n",t_in_nano_seconds, system->unit_converter->permeability_to_SI(permeability));
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
    system->timer->start_mpi_reduce();
    MPI_Reduce(&kinetic_energy_local, &kinetic_energy,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    system->timer->end_mpi_reduce();

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
    double flux_local = 0;
    double elapsed_time_this_run = system->t - system->t0;
    // Either, the gas is pressure driven or gravity driven. We measure flux different in these two cases.
    if(settings->maintain_pressure) {
        flux_local = system->flux_count / elapsed_time_this_run;
    } else {
        flux_local = system->mover->count_periodic[settings->flow_direction] / elapsed_time_this_run;
    }

    MPI_Reduce(&flux_local, &flux, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

void StatisticsSampler::sample_permeability() {
    if(system->steps == permeability_sampled_at) return;
    permeability_sampled_at = system->steps;
    if(settings->flow_direction<0) return;

    sample_flux();
    if(system->myid == 0) {
        double volume_per_molecule = system->volume / system->num_molecules_global;
        double viscosity_dsmc_units = system->unit_converter->viscosity_from_SI(settings->viscosity);
        double volume_flux = flux*volume_per_molecule;
        double L = system->length[settings->flow_direction];
        double mass_density = system->density*settings->mass;

        double area = 1;
        for(int a=0;a<3;a++) {
            if(a != settings->flow_direction) area *= system->length[a];
        }

        if(settings->maintain_pressure) {
            // Expression from Darcy's law of gases
            double pressure_in_reservoir_a = system->unit_converter->pressure_from_SI(settings->pressure_A);
            double pressure_in_reservoir_b = system->unit_converter->pressure_from_SI(settings->pressure_B);
            permeability = 2*pressure_in_reservoir_b*volume_flux*L*viscosity_dsmc_units / (area * (pressure_in_reservoir_a*pressure_in_reservoir_a - pressure_in_reservoir_b*pressure_in_reservoir_b));
        } else {
            double pressure_in_reservoir_b = system->density*system->temperature;
            double pressure_in_reservoir_a = pressure_in_reservoir_b + mass_density*settings->gravity*system->length[settings->flow_direction];
            permeability = 2*pressure_in_reservoir_b*volume_flux*L*viscosity_dsmc_units / (area * (pressure_in_reservoir_a*pressure_in_reservoir_a - pressure_in_reservoir_b*pressure_in_reservoir_b));
            // permeability = volume_flux*L*viscosity_dsmc_units / (mass_density*system->length[settings->flow_direction]*settings->gravity*area);
        }
    }
}

void StatisticsSampler::sample_stats_across_channel() {
    for(int i=0;i<system->num_molecules_local;i++) {
        double pos_across_the_channel = system->r[3*i + 1];
        int bin_index = num_bins*(pos_across_the_channel/system->length[1]);
        count_across_channel[bin_index]++;
        kinetic_energy_across_channel[bin_index] += 0.5*settings->mass*(system->v[3*i+0]*system->v[3*i+0] + system->v[3*i+1]*system->v[3*i+1] + system->v[3*i+2]*system->v[3*i+2]);
    }
}

void StatisticsSampler::sample_velocity_distribution() {
    if(system->steps == velocity_distribution_sampled_at) return;
    velocity_distribution_sampled_at = system->steps;
    density_sampled_at = system->steps;

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_x = system->r[3*i+0]*system->one_over_length[0]*num_bins;
        int bin_y = system->r[3*i+1]*system->one_over_length[1]*num_bins;
        int bin_z = system->r[3*i+2]*system->one_over_length[2]*num_bins;

        int index = bin_y;

        if(settings->velocity_profile_type.compare("area") == 0) {
            index = bin_x*num_bins + bin_y;
        }

        velocity_distribution[index] += system->v[3*i+2];
        velocity_distribution_count[index]++;
    }
}

void StatisticsSampler::gather_velocity_distribution() {
    double *velocity_distribution_global = new double[velocity_distribution.size()];
    unsigned long *velocity_distribution_count_global = new unsigned long[velocity_distribution.size()];
    memset(velocity_distribution_global,0,velocity_distribution.size());
    memset(velocity_distribution_count_global,0,velocity_distribution.size());

    system->timer->start_mpi_reduce();
    MPI_Reduce(&velocity_distribution[0],velocity_distribution_global,velocity_distribution.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&velocity_distribution_count[0],velocity_distribution_count_global,velocity_distribution.size(),MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    system->timer->end_mpi_reduce();

    if(system->myid==0) {
        for(int i=0; i<velocity_distribution.size(); i++) {
            // If we sample the 2d field, the averages are done in cpp code
            if(velocity_distribution_count_global[i]>0) velocity_distribution_global[i] /= velocity_distribution_count_global[i];
            fprintf(system->io->velocity_file,"%f ",system->unit_converter->velocity_to_SI(velocity_distribution_global[i]));
            fprintf(system->io->density_file,"%f ",(double)velocity_distribution_count_global[i] / num_samples);
        }

        fprintf(system->io->velocity_file,"\n");
    }

    delete velocity_distribution_global;
    delete velocity_distribution_count_global;
}

void StatisticsSampler::gather_stats_across_channel() {
    sample_temperature();
    vector<double> kinetic_energy_across_channel_global;
    vector<unsigned long> count_across_channel_global;
    kinetic_energy_across_channel_global.resize(num_bins, 0);
    count_across_channel_global.resize(num_bins, 0);

    system->timer->start_mpi_reduce();
    MPI_Reduce(&kinetic_energy_across_channel[0],&kinetic_energy_across_channel_global[0],kinetic_energy_across_channel_global.size(),MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&count_across_channel[0],&count_across_channel_global[0],count_across_channel_global.size(),MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    system->timer->end_mpi_reduce();

    if(system->myid == 0) {
        double volume_per_bin = system->volume_global / count_across_channel_global.size();

        for(int i=0; i<count_across_channel_global.size(); i++) {
            double temperature_in_bin = 0;
            if(count_across_channel_global[i] > 0) temperature_in_bin = kinetic_energy_across_channel_global[i]*2.0/(3*count_across_channel_global[i]);
            double density = system->atoms_per_molecule*count_across_channel_global[i]/volume_per_bin;
            double pressure = density*temperature_in_bin; // Ideal gas law

            fprintf(system->io->linear_density_file,"%E ", system->unit_converter->number_density_to_SI(density));
            fprintf(system->io->linear_pressure_file,"%E ",system->unit_converter->pressure_to_SI(pressure));
            fprintf(system->io->linear_temperature_file,"%E ",system->unit_converter->temperature_to_SI(temperature_in_bin));
        }

        fprintf(system->io->linear_density_file,"\n");
        fprintf(system->io->linear_pressure_file,"\n");
        fprintf(system->io->linear_temperature_file,"\n");
    }
}

void StatisticsSampler::finalize() {
    gather_velocity_distribution();
    gather_stats_across_channel();
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
