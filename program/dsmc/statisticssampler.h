#pragma once

class System;
class Settings;
#include <stdio.h>
#include <fstream>
#include <cinifile.h>

class StatisticsSampler {

private:
	System *system;
    Settings *settings;
    long temperature_sampled_at;
    long kinetic_energy_sampled_at;
    long velocity_distribution_sampled_at;
    long density_sampled_at;
    long linear_density_sampled_at;
    long flux_sampled_at;
    long permeability_sampled_at;
    long num_samples;

    int num_bins;
    vector<double> velocity_distribution;
    vector<unsigned long> velocity_distribution_count;
    vector<double> kinetic_energy_across_channel;
    vector<unsigned long> count_across_channel;

public:
    StatisticsSampler(System *system);
    void sample();
    void sample_kinetic_energy();
    void sample_temperature();
    void sample_velocity_distribution();
    void sample_flux();
    void sample_permeability();
    void sample_stats_across_channel();
    void gather_velocity_distribution();
    void gather_stats_across_channel();
    void finalize();

    unsigned long collisions;
    unsigned long wall_collisions;
    double kinetic_energy;
    double temperature;
    double flux;
    double permeability;
};
