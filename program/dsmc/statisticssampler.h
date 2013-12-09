#pragma once

#include <stdio.h>
#include <fstream>
#include <cinifile.h>
#include <map>
#include <string>
#include <statisticalproperty.h>
#include <vector>

class System;
class Settings;
using std::vector;

class StatisticsSampler {

private:
	System *system;
    vector<StatisticalProperty*> statistical_properties;
public:
    StatisticsSampler(System *system);
    void sample();
    void finalize();
    MeasureEnergy *energy;
    MeasureTemperature *temperature;
    MeasurePressure *pressure;
    MeasureFlux *flux;
    MeasureVolumetricFlowRate *volumetric_flow_rate;
    MeasurePermeability *permeability;
    MeasureCount *count;
    MeasureVelocityDistributionPoiseuille *velocity;
    MeasureTemperatureDistribution *temperature_distribution;
    MeasurePressureDistribution *pressure_distribution;
};
