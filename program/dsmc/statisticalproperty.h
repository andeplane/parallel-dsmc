#pragma once
#include <algorithm>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <statisticalvalue.h>

class UnitConverter;
class System;
using std::max;

class StatisticalProperty
{
protected:
    FILE *file;
    int myid;
    int interval;
    int last_sample;
public:

    StatisticalProperty(int myid_, int interval_, FILE *file_);
    virtual void update(System *system) = 0;
    virtual void finalize(UnitConverter *unit_converter) {}
    virtual void finalize(System *system) {}
    virtual void resize(int number_of_bins) {}
};

class MeasureEnergy : public StatisticalProperty {
public:
    StatisticalValue<double> value;
    MeasureEnergy(FILE *file_, int myid_, int interval_);

    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    double get_current_value();
};

class MeasureTemperature : public StatisticalProperty {
protected:
    MeasureEnergy *energy;
public:
    StatisticalValue<double> value;
    MeasureTemperature(FILE *file_, int myid_, int interval_, MeasureEnergy *energy_);
    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    double get_current_value();
};

class MeasureFlux : public StatisticalProperty {
public:
    StatisticalValue<double> value;
    MeasureFlux(FILE *file_, int myid_, int interval_);
    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    double get_current_value();
};

class MeasurePressure : public StatisticalProperty {
protected:
    MeasureTemperature *temperature;
public:
    StatisticalValue<double> value;
    MeasurePressure(FILE *file_, int myid_, int interval_, MeasureTemperature *temperature_);
    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    double get_current_value();
};

class MeasurePermeability : public StatisticalProperty {
protected:
    MeasureFlux *flux;
public:
    StatisticalValue<double> value;
    MeasurePermeability(FILE *file_, int myid_, int interval_, MeasureFlux *flux_);
    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    double get_current_value();
};

class MeasureCount : public StatisticalProperty {
public:
    StatisticalValue<unsigned long> value;
    MeasureCount(FILE *file_, int myid_, int interval_, int bins);
    virtual void update(System *system);
    vector<unsigned long> get_current_value();
    vector<unsigned long> get_sum();
    virtual void resize(int number_of_bins);
};

class MeasureVelocityDistributionPoiseuille : public StatisticalProperty {
protected:
    MeasureCount *count;
public:
    StatisticalValue<double> value;
    MeasureVelocityDistributionPoiseuille(FILE *file_, int myid_, int interval_, int bins, MeasureCount *count);
    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    vector<double> get_average();
    virtual void resize(int number_of_bins);
};

class MeasureTemperatureDistribution : public StatisticalProperty {
protected:
    MeasureCount *count;
public:
    StatisticalValue<double> value;
    MeasureTemperatureDistribution(FILE *file_, int myid_, int interval_, int bins, MeasureCount *count);
    virtual void update(System *system);
    virtual void finalize(UnitConverter *unit_converter);
    vector<double> get_average();
    virtual void resize(int number_of_bins);
};

class MeasurePressureDistribution : public StatisticalProperty {
protected:
    MeasureCount *count;
    MeasureTemperatureDistribution *temperature_distribution;
public:
    StatisticalValue<double> value;
    MeasurePressureDistribution(FILE *file_, int myid_, int interval_, int bins, MeasureCount *count_, MeasureTemperatureDistribution *temperature_distribution_);
    virtual void update(System *system);
    virtual void finalize(System *system);
    vector<double> get_average(System *system);
    virtual void resize(int number_of_bins);
};
