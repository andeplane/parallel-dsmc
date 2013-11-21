#include <system.h>
#include <statisticssampler.h>
#include <statisticalproperty.h>
#include <dsmc_io.h>
#include <settings.h>
#include <dsmctimer.h>
#include <colliderbase.h>
#include <moleculemover.h>

using std::iterator;
using std::vector;

StatisticsSampler::StatisticsSampler(System *system_) {
    system = system_;

    energy = new MeasureEnergy(system->io->energy_file,system->myid,system->settings->statistics_interval);
    temperature = new MeasureTemperature(system->io->temperature_file,system->myid, system->settings->statistics_interval, energy);
    flux = new MeasureFlux(system->io->flux_file,system->myid,system->settings->statistics_interval);
    pressure = new MeasurePressure(system->io->pressure_file,system->myid,system->settings->statistics_interval, temperature);
    permeability = new MeasurePermeability(system->io->permeability_file, system->myid, system->settings->statistics_interval, flux);
    count = new MeasureCount(system->io->num_molecules_file, system->myid, system->settings->statistics_interval, system->settings->sampling_bins);
    velocity = new MeasureVelocityDistributionPoiseuille(system->io->velocity_file, system->myid, system->settings->statistics_interval,system->settings->sampling_bins, count);

    statistical_properties.push_back(energy);
    statistical_properties.push_back(temperature);
    statistical_properties.push_back(flux);
    statistical_properties.push_back(pressure);
    statistical_properties.push_back(permeability);
    statistical_properties.push_back(count);
    statistical_properties.push_back(velocity);
}


void StatisticsSampler::sample() {
    system->timer->start_sample();

    for(int i=0; i<statistical_properties.size(); i++) {
        StatisticalProperty *value = statistical_properties[i];
        value->update(system);
    }

    if(system->steps % system->settings->statistics_interval > 0) return;

    unsigned long collisions = 0;
    unsigned long wall_collisions = 0;

    system->timer->start_mpi_reduce();
    MPI_Reduce(&system->collisions,&collisions,1,MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);
    MPI_Reduce(&system->mover->surface_collider->num_collisions,&wall_collisions,1,MPI_UNSIGNED_LONG, MPI_SUM,0, MPI_COMM_WORLD);
    system->timer->end_mpi_reduce();

    if(system->myid==0) {
        cout << system->steps << "   t=" << system->t_in_nano_seconds() << "ns   T=" << system->unit_converter->temperature_to_SI(temperature->get_current_value()) << "K   Collisions: " <<  collisions <<   "   Wall collisions: " << wall_collisions << "   Pressure: " << system->unit_converter->pressure_to_SI(pressure->get_current_value()) <<  "Pa   Molecules: " << system->num_molecules_global << endl ;
    }

    system->timer->end_sample();
}

void StatisticsSampler::finalize() {
    for(int i=0; i<statistical_properties.size(); i++) {
        StatisticalProperty *value = statistical_properties[i];
        value->finalize(system->unit_converter);
    }
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
