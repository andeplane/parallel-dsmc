#pragma once

class System;

class DSMCTimer
{
public:
    DSMCTimer();
    double t0;

    double moving_t0;
    double moving;
    double moving_global;
    void start_moving();
    void end_moving();
    double fraction_moving();

    double colliding_t0;
    double colliding;
    double colliding_global;
    void start_colliding();
    void end_colliding();
    double fraction_colliding();

    double mpi_t0;
    double mpi;
    double mpi_global;
    void start_mpi();
    void end_mpi();
    double fraction_mpi();

    double system_initialize_t0;
    double system_initialize;
    double system_initialize_global;
    void start_system_initialize();
    void end_system_initialize();
    double fraction_system_initialize();

    double io_t0;
    double io;
    double io_global;
    void start_io();
    void end_io();
    double fraction_io();

    double sample_t0;
    double sample;
    double sample_global;
    void start_sample();
    void end_sample();
    double fraction_sample();

    double accelerate_t0;
    double accelerate;
    double accelerate_global;
    void start_accelerate();
    void end_accelerate();
    double fraction_accelerate();

    double pressure_t0;
    double pressure;
    double pressure_global;
    void start_pressure();
    void end_pressure();
    double fraction_pressure();

    void gather_all_nodes(System *system);
};
