#include <system.h>
#include <statisticalproperty.h>
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
        fprintf(file, "%f %E\n",system->t_in_nano_seconds(), system->unit_converter->energy_to_eV(kinetic_energy_global));
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
        double temperature = 2.0/3.0*energy->get_current_value() / system->get_number_of_atoms_global();

        fprintf(file, "%f %f\n",system->t_in_nano_seconds(), system->unit_converter->temperature_to_SI(temperature));
        value.add_value(temperature);
    }
}

double MeasureTemperature::get_current_value() {
    return value.get_current_value()[0];
}

void MeasureTemperature::finalize(UnitConverter *unit_converter) {

}

MeasureNumberFlowRate::MeasureNumberFlowRate(FILE *file_, int myid_, int interval_) :
    StatisticalProperty(myid_,interval_, file_)
{

}

void MeasureNumberFlowRate::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    double flux_local = 0;
    double flux_global = 0;

    flux_local = system->mover->count_periodic[system->settings->flow_direction] / system->get_elapsed_time_this_run();
    MPI_Reduce(&flux_local, &flux_global, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if(myid==0) {
        fprintf(file, "%f %E\n",system->t_in_nano_seconds(), flux_global);
        value.add_value(flux_global);
    }
}

double MeasureNumberFlowRate::get_current_value() {
    return value.get_current_value()[0];
}

void MeasureNumberFlowRate::finalize(UnitConverter *unit_converter) {

}



MeasureVolumetricFlowRate::MeasureVolumetricFlowRate(FILE *file_, int myid_, int interval_, MeasureNumberFlowRate *number_flow_rate_) :
    StatisticalProperty(myid_,interval_, file_),
    number_flow_rate(number_flow_rate_)
{

}

void MeasureVolumetricFlowRate::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;
    number_flow_rate->update(system);

    if(myid==0) {
        // Global values are already calculated in flux
        double volume_per_particle = system->volume_global / system->num_molecules_global;
        double volumetric_flow_rate = number_flow_rate->get_current_value() * volume_per_particle;

        fprintf(file, "%f %E\n",system->t_in_nano_seconds(), volumetric_flow_rate);
        value.add_value(volumetric_flow_rate);
    }
}

double MeasureVolumetricFlowRate::get_current_value() {
    return value.get_current_value()[0];
}

void MeasureVolumetricFlowRate::finalize(UnitConverter *unit_converter) {

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
        double pressure = system->get_number_of_atoms_global() / system->volume_global * temperature->get_current_value();
        fprintf(file, "%f %E\n",system->t_in_nano_seconds(), system->unit_converter->pressure_to_SI(pressure));
        value.add_value(pressure);
    }
}

double MeasurePressure::get_current_value() {
    return value.get_current_value()[0];
}

void MeasurePressure::finalize(UnitConverter *unit_converter) {

}

MeasurePermeability::MeasurePermeability(FILE *file_, int myid_, int interval_, MeasureVolumetricFlowRate *volumetric_flow_rate_, MeasurePressure *pressure_) :
    StatisticalProperty(myid_, interval_, file_),
    volumetric_flow_rate(volumetric_flow_rate_),
    pressure(pressure_)
{

}

void MeasurePermeability::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    volumetric_flow_rate->update(system);
    pressure->update(system);

    if(myid==0) {
        double viscosity_dsmc_units = system->unit_converter->viscosity_from_SI(system->settings->viscosity);
        double volumetric_flow_rate_value = volumetric_flow_rate->get_current_value();
        double L = system->length[system->settings->flow_direction];
        double mass_density = system->density*system->settings->mass;

        double area = 1;
        for(int a=0;a<3;a++) {
            if(a != system->settings->flow_direction) area *= system->length[a];
        }

        // double average_pressure = pressure->get_current_value();
        // double pressure_in_reservoir_a = average_pressure;
        // double pressure_in_reservoir_b = pressure_in_reservoir_a - mass_density*system->settings->gravity*system->length[system->settings->flow_direction];
        // double permeability = 2*average_pressure*volumetric_flow_rate_value*L*viscosity_dsmc_units / (area * (pressure_in_reservoir_a*pressure_in_reservoir_a - pressure_in_reservoir_b*pressure_in_reservoir_b));
        double permeability = volumetric_flow_rate_value*viscosity_dsmc_units / (area * system->density * system->settings->gravity);
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

vector<unsigned long> MeasureCount::get_sum() {
    return value.get_sum();
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
    count->update(system);

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_y = system->r[3*i+1]*system->one_over_length[1]*value.number_of_bins;

        int index = bin_y;

        velocity_distribution[index] += system->v[3*i+2];
    }

    value.add_value(velocity_distribution);
    velocity_distribution.clear();
}

vector<double> MeasureVelocityDistributionPoiseuille::get_average() {
    vector<double> values_global(value.number_of_bins,0);
    vector<unsigned long> count_values_global(value.number_of_bins,0);

    vector<double> values(value.get_sum());
    vector<unsigned long> count_values = count->get_sum();

    MPI_Reduce(&values[0],&values_global[0],value.number_of_bins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&count_values[0],&count_values_global[0],value.number_of_bins,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);

    if(myid==0) {
        for(int i=0; i<values_global.size(); i++) {
            values_global[i] /= max(count_values_global[i],(unsigned long)1); // Normalize
        }
    }

    count_values_global.clear();
    values.clear();
    count_values.clear();

    return values_global;
}

void MeasureVelocityDistributionPoiseuille::finalize(UnitConverter *unit_converter) {
    vector<double> velocity_distribution = get_average();
    if(myid!=0) return;
    for(int i=0; i<velocity_distribution.size(); i++) {
        fprintf(file,"%E ", unit_converter->velocity_to_SI(velocity_distribution[i]));
    }
}

void MeasureVelocityDistributionPoiseuille::resize(int number_of_bins) {
    value.resize(number_of_bins);
}

MeasureTemperatureDistribution::MeasureTemperatureDistribution(FILE *file_, int myid_, int interval_, int number_of_bins, MeasureCount *count_) :
    StatisticalProperty(myid_, interval_, file_)
{
    count = count_;
    resize(number_of_bins);
}

void MeasureTemperatureDistribution::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    vector<double> temperature_distribution(value.number_of_bins,0);
    count->update(system);

    for(int i=0; i<system->num_molecules_local; i++) {
        int bin_y = system->r[3*i+1]*system->one_over_length[1]*value.number_of_bins;
        int index = bin_y;

        temperature_distribution[index] += 0.5*system->settings->mass*system->v[3*i+2]*system->v[3*i+2];
    }

    value.add_value(temperature_distribution);
    temperature_distribution.clear();
}

vector<double> MeasureTemperatureDistribution::get_average() {
    vector<double> values_global(value.number_of_bins,0);
    vector<unsigned long> count_values_global(value.number_of_bins,0);

    vector<double> values(value.get_sum());
    vector<unsigned long> count_values = count->get_sum();

    MPI_Reduce(&values[0],&values_global[0],value.number_of_bins,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&count_values[0],&count_values_global[0],value.number_of_bins,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);

    if(myid==0) {
        for(int i=0; i<values_global.size(); i++) {
            values_global[i] /= max(count_values_global[i],(unsigned long)1); // Normalize
            values_global[i] *= 2.0/3.0; // Energies are already measured per particle
        }
    }

    count_values_global.clear();
    values.clear();
    count_values.clear();

    return values_global;
}

void MeasureTemperatureDistribution::finalize(UnitConverter *unit_converter) {
    vector<double> temperature_distribution = get_average();

    if(myid!=0) return;
    for(int i=0; i<temperature_distribution.size(); i++) {
        fprintf(file,"%E ", unit_converter->temperature_to_SI(temperature_distribution.at(i)));
    }
}

void MeasureTemperatureDistribution::resize(int number_of_bins) {
    value.resize(number_of_bins);
}

MeasurePressureDistribution::MeasurePressureDistribution(FILE *file_, int myid_, int interval_, int number_of_bins, MeasureCount *count_, MeasureTemperatureDistribution *temperature_distribution_) :
    StatisticalProperty(myid_, interval_, file_)
{
    count = count_;
    temperature_distribution = temperature_distribution_;
    resize(number_of_bins);
}

void MeasurePressureDistribution::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    count->update(system);
    temperature_distribution->update(system);
}

vector<double> MeasurePressureDistribution::get_average(System *system) {
    vector<double> values_global(value.number_of_bins,0);
    vector<unsigned long> count_values_global(value.number_of_bins,0);

    vector<unsigned long> count_values = count->get_sum();
    vector<double> temperature_values = temperature_distribution->get_average();

    MPI_Reduce(&count_values[0],&count_values_global[0],value.number_of_bins,MPI_UNSIGNED_LONG,MPI_SUM,0,MPI_COMM_WORLD);

    if(myid==0) {
        double volume_per_bin = system->volume_global / value.number_of_bins; //
        for(int i=0; i<values_global.size(); i++) {
            double density = count_values_global[i] / volume_per_bin;
            double temperature = temperature_values[i];

            values_global[i] = density*temperature;
        }
    }

    count_values_global.clear();
    count_values.clear();
    temperature_values.clear();

    return values_global;
}

void MeasurePressureDistribution::finalize(System *system) {
    vector<double> pressure_distribution = get_average(system);

    if(myid!=0) return;
    for(int i=0; i<pressure_distribution.size(); i++) {
        fprintf(file,"%E ", system->unit_converter->pressure_to_SI(pressure_distribution.at(i)));
    }
}

void MeasurePressureDistribution::resize(int number_of_bins) {
    value.resize(number_of_bins);
}

MeasureVelocityDistribution::MeasureVelocityDistribution(FILE *file_, int myid_, int interval_, int number_of_bins, System *system) :
    StatisticalProperty(myid_, interval_, file_)
{
    resize(number_of_bins);
    std_dev = sqrt(system->settings->mass/system->temperature);
    v_min = -2*std_dev;
    v_max = 2*std_dev;
    bin_size = (v_max - v_min) / value.number_of_bins;
}

void MeasureVelocityDistribution::update(System *system) {
    if((system->steps % interval) || system->steps == last_sample) return;
    last_sample = system->steps;

    vector<long> velocity_distribution(value.number_of_bins,0);

    for(unsigned int i=0; i<system->num_molecules_local; i++) {
        int bin_x = (system->v[3*i + 0] + v_min) / bin_size;
        int bin_y = (system->v[3*i + 1] + v_min) / bin_size;
        int bin_z = (system->v[3*i + 2] + v_min) / bin_size;

        if(bin_x < 0 || bin_x > value.number_of_bins-1) continue;

        velocity_distribution[bin_x]++;
    }

    value.add_value(velocity_distribution);
    velocity_distribution.clear();
}

void MeasureVelocityDistribution::finalize(UnitConverter *unit_converter) {
    vector<long> values_global(value.number_of_bins,0);
    vector<long> values(value.get_sum());
    MPI_Reduce(&values[0],&values_global[0],value.number_of_bins,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
    values.clear();

    if(myid!=0) return;

    for(int i=0; i<value.number_of_bins; i++) {
        double v = v_min + i*bin_size;

        fprintf(file,"%f %ld\n", unit_converter->velocity_to_SI(v),values_global[i]);
    }


}

void MeasureVelocityDistribution::resize(int number_of_bins) {
    value.resize(number_of_bins);
}
