#include <statisticalproperty.h>
#include <system.h>
#include <dsmctimer.h>
#include <unitconverter.h>
#include <moleculemover.h>
#include <settings.h>

StatisticalProperty::StatisticalProperty(int myid_, int interval_, FILE *file_) :
    myid(myid_),
    interval(interval_),
    last_sample(0),
    file(file_)
{ }

MeasureEnergy::MeasureEnergy(FILE *file_, int myid_, int interval_) :
    StatisticalProperty(myid_,interval_, file_)
{

}

void MeasureEnergy::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    double kinetic_energy_local = 0;
    double kinetic_energy_global = 0;

    for(unsigned int i=0;i<system->num_molecules_local;i++) {
        kinetic_energy_local += (system->v[3*i+0]*system->v[3*i+0] + system->v[3*i+1]*system->v[3*i+1] + system->v[3*i+2]*system->v[3*i+2]);
    }
    kinetic_energy_local *= 0.5*system->settings->mass*system->atoms_per_molecule;

    system->timer->start_mpi_reduce();
    MPI_Reduce(&kinetic_energy_local, &kinetic_energy_global,1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    system->timer->end_mpi_reduce();

    if(myid==0) {
        fprintf(file, "%f %f\n",system->t_in_nano_seconds(), system->unit_converter->energy_to_eV(kinetic_energy_global));
        value.add_value(kinetic_energy_global);
    }
}

double MeasureEnergy::get_current_value() {
    return value.get_current_value()[0];

}

void MeasureEnergy::finalize(UnitConverter *unit_converter) {

}

MeasureTemperature::MeasureTemperature(FILE *file_, int myid_, int interval_, MeasureEnergy *energy_) :
    StatisticalProperty(myid_,interval_, file_),
    energy(energy_)
{

}

void MeasureTemperature::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    energy->update(system);

    if(myid==0) {
        double temperature = 2.0/3.0*energy->get_current_value() / (system->num_molecules_global*system->atoms_per_molecule);

        fprintf(file, "%f %f\n",system->t_in_nano_seconds(), system->unit_converter->temperature_to_SI(temperature));
        value.add_value(temperature);
    }
}

double MeasureTemperature::get_current_value() {
    return value.get_current_value()[0];
}

void MeasureTemperature::finalize(UnitConverter *unit_converter) {

}

MeasureFlux::MeasureFlux(FILE *file_, int myid_, int interval_) :
    StatisticalProperty(myid_,interval_, file_)
{

}

void MeasureFlux::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    double flux_local = 0;
    double flux_global = 0;
    double elapsed_time_this_run = system->t - system->t0;

    flux_local = system->mover->count_periodic[system->settings->flow_direction] / elapsed_time_this_run;
    MPI_Reduce(&flux_local, &flux_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myid==0) {
        fprintf(file, "%f %f\n",system->t_in_nano_seconds(), flux_global);
        value.add_value(flux_global);
    }
}

double MeasureFlux::get_current_value() {
    return value.get_current_value()[0];
}

void MeasureFlux::finalize(UnitConverter *unit_converter) {

}

MeasurePressure::MeasurePressure(FILE *file_, int myid_, int interval_, MeasureTemperature *temperature_) :
    StatisticalProperty(myid_, interval_, file_),
    temperature(temperature_)
{

}

void MeasurePressure::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    temperature->update(system);

    if(myid==0) {
        double pressure = system->num_molecules_global*system->atoms_per_molecule / system->volume_global * temperature->get_current_value();
        fprintf(file, "%f %f\n",system->t_in_nano_seconds(), system->unit_converter->pressure_to_SI(pressure));
        value.add_value(pressure);
    }
}

double MeasurePressure::get_current_value() {
    return value.get_current_value()[0];
}

void MeasurePressure::finalize(UnitConverter *unit_converter) {

}

MeasurePermeability::MeasurePermeability(FILE *file_, int myid_, int interval_, MeasureFlux *flux_) :
    StatisticalProperty(myid_, interval_, file_),
    flux(flux_)
{

}

void MeasurePermeability::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    flux->update(system);

    if(myid==0) {
        double volume_per_molecule = system->volume_global / system->num_molecules_global;
        double viscosity_dsmc_units = system->unit_converter->viscosity_from_SI(system->settings->viscosity);
        double volume_flux = flux->get_current_value()*volume_per_molecule;
        double L = system->length[system->settings->flow_direction];
        double mass_density = system->density*system->settings->mass;

        double area = 1;
        for(int a=0;a<3;a++) {
            if(a != system->settings->flow_direction) area *= system->length[a];
        }

        double pressure_in_reservoir_b = system->density*system->temperature;
        double pressure_in_reservoir_a = pressure_in_reservoir_b + mass_density*system->settings->gravity*system->length[system->settings->flow_direction];
        double permeability = 2*pressure_in_reservoir_b*volume_flux*L*viscosity_dsmc_units / (area * (pressure_in_reservoir_a*pressure_in_reservoir_a - pressure_in_reservoir_b*pressure_in_reservoir_b));
        value.add_value(permeability);
    }
}

double MeasurePermeability::get_current_value() {
    return value.get_current_value()[0];
}

void MeasurePermeability::finalize(UnitConverter *unit_converter) {
    if(myid!=0) return;

    fprintf(file, "%E\n",unit_converter->permeability_to_darcy(value.get_current_value()[0]));
}

MeasureCount::MeasureCount(FILE *file_, int myid_, int interval_, int number_of_bins) :
    StatisticalProperty(myid_, interval_, file_)
{
    this->resize(number_of_bins);
}

void MeasureCount::resize(int number_of_bins) {
    value.resize(number_of_bins);
}

void MeasureCount::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    vector<unsigned long> count(value.number_of_bins, 0);

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_y = system->r[3*i+1]*system->one_over_length[1]*value.number_of_bins;

        int index = bin_y;
        count[index]++;
    }

    value.add_value(count);
    count.clear();
}

vector<unsigned long> MeasureCount::get_current_value() {
    return value.get_current_value();
}

MeasureVelocityDistributionPoiseuille::MeasureVelocityDistributionPoiseuille(FILE *file_, int myid_, int interval_, int number_of_bins, MeasureCount *count_) :
    StatisticalProperty(myid_, interval_, file_)
{
    count = count_;
    resize(number_of_bins);
}

void MeasureVelocityDistributionPoiseuille::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    vector<double> velocity_distribution(value.number_of_bins,0);

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_y = system->r[3*i+1]*system->one_over_length[1]*value.number_of_bins;

        int index = bin_y;

        velocity_distribution[index] += system->v[3*i+2];
    }

    value.add_value(velocity_distribution);
    velocity_distribution.clear();
}

vector<double> MeasureVelocityDistributionPoiseuille::get_current_value() {
    vector<double> values_global(value.number_of_bins,0);
    vector<unsigned long> count_values_global(value.number_of_bins,0);

    vector<double> values(value.current_value);
    vector<unsigned long> count_values = count->get_current_value();


    MPI_Reduce(&values[0],&values_global[0],value.number_of_bins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&count_values[0],&count_values_global[0],value.number_of_bins,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);

    if(myid==0) {
        for(int i=0; i<values_global.size(); i++) {
            values_global[i] /= max(count_values_global[i],(unsigned long)1); // Normalize
        }
    }

    count_values_global.clear();
    values.clear();
    count_values_global.clear();

    return values_global;
}

void MeasureVelocityDistributionPoiseuille::finalize(UnitConverter *unit_converter) {
    vector<double> velocity_distribution = get_current_value();
    if(myid!=0) return;
    for(int i=0; i<velocity_distribution.size(); i++) {
        fprintf(file,"%E ", unit_converter->velocity_to_SI(velocity_distribution[i]));
    }
}

void MeasureVelocityDistributionPoiseuille::resize(int number_of_bins) {
    value.resize(number_of_bins);
}
