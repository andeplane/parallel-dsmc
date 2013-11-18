#pragma once
#include <algorithm>
#include <math.h>
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
    virtual void finalize(UnitConverter *unit_converter) = 0;
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
