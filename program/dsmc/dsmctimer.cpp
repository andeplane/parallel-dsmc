#include <dsmctimer.h>
#include <mpi.h>
#include <system.h>
#include <topology.h>
#include <dsmc_io.h>
#include <settings.h>
#include <iostream>
#include <fstream>
using std::endl;

DSMCTimer::DSMCTimer() {
    t0 = MPI_Wtime();
    io = 0;
    moving = 0;
    colliding = 0;
    mpi = 0;
    mpi_reduce = 0;
    sample = 0;
    accelerate = 0;
    pressure = 0;
    system_initialize = 0;
}

void DSMCTimer::start_moving() {
    moving_t0 = MPI_Wtime();
}

void DSMCTimer::end_moving() {
    moving += MPI_Wtime() - moving_t0;
}

double DSMCTimer::fraction_moving() {
    double t1 = MPI_Wtime();
    return moving_global/(t1-t0);
}

void DSMCTimer::start_colliding() {
    colliding_t0 = MPI_Wtime();
}

void DSMCTimer::end_colliding() {
    colliding += MPI_Wtime() - colliding_t0;
}

double DSMCTimer::fraction_colliding() {
    double t1 = MPI_Wtime();
    return colliding_global/(t1-t0);
}

void DSMCTimer::start_mpi_reduce() {
    mpi_reduce_t0 = MPI_Wtime();
}

void DSMCTimer::end_mpi_reduce() {
    mpi_reduce += MPI_Wtime() - mpi_reduce_t0;
}

double DSMCTimer::fraction_mpi_reduce() {
    double t1 = MPI_Wtime();
    return mpi_reduce_global/(t1-t0);
}

void DSMCTimer::start_mpi() {
    mpi_t0 = MPI_Wtime();
}

void DSMCTimer::end_mpi() {
    mpi += MPI_Wtime() - mpi_t0;
}

double DSMCTimer::fraction_mpi() {
    double t1 = MPI_Wtime();
    return mpi_global/(t1-t0);
}

void DSMCTimer::start_system_initialize() {
    system_initialize_t0 = MPI_Wtime();
}

void DSMCTimer::end_system_initialize() {
    system_initialize += MPI_Wtime() - system_initialize_t0;
}

double DSMCTimer::fraction_system_initialize() {
    double t1 = MPI_Wtime();
    return system_initialize_global/(t1-t0);
}

void DSMCTimer::start_io() {
    io_t0 = MPI_Wtime();
}

void DSMCTimer::end_io() {
    io += MPI_Wtime() - io_t0;
}

double DSMCTimer::fraction_io() {
    double t1 = MPI_Wtime();
    return io_global/(t1-t0);
}

void DSMCTimer::start_sample() {
    sample_t0 = MPI_Wtime();
}

void DSMCTimer::end_sample() {
    sample += MPI_Wtime() - sample_t0;
}

double DSMCTimer::fraction_sample() {
    double t1 = MPI_Wtime();
    return sample_global/(t1-t0);
}

void DSMCTimer::start_accelerate() {
    accelerate_t0 = MPI_Wtime();
}

void DSMCTimer::end_accelerate() {
    accelerate += MPI_Wtime() - accelerate_t0;
}

double DSMCTimer::fraction_accelerate() {
    double t1 = MPI_Wtime();
    return accelerate_global/(t1-t0);
}

void DSMCTimer::start_pressure() {
    pressure_t0 = MPI_Wtime();
}

void DSMCTimer::end_pressure() {
    pressure += MPI_Wtime() - pressure_t0;
}

double DSMCTimer::fraction_pressure() {
    double t1 = MPI_Wtime();
    return pressure_global/(t1-t0);
}

void DSMCTimer::gather_all_nodes(System &system) {
    colliding_global = 0;
    moving_global = 0;
    mpi_global = 0;
    mpi_reduce_global = 0;
    system_initialize_global = 0;
    accelerate_global = 0;
    pressure_global = 0;
    sample_global = 0;
    io_global = 0;

    MPI_Reduce(&colliding,&colliding_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&moving,&moving_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&mpi,&mpi_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&system_initialize,&system_initialize_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&accelerate,&accelerate_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&pressure,&pressure_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&sample,&sample_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&io,&io_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&mpi_reduce,&mpi_reduce_global,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    colliding_global /= system.topology->num_processors;
    mpi_reduce_global /= system.topology->num_processors;
    moving_global /= system.topology->num_processors;
    mpi_global /= system.topology->num_processors;
    system_initialize_global /= system.topology->num_processors;
    accelerate_global /= system.topology->num_processors;
    pressure_global /= system.topology->num_processors;
    sample_global /= system.topology->num_processors;
    io_global /= system.topology->num_processors;
}

void DSMCTimer::save_to_file(System &system) {
    colliding_global = colliding;
    moving_global = moving;
    mpi_global = mpi;
    system_initialize_global = system_initialize;
    accelerate_global = accelerate;
    pressure_global = pressure;
    sample_global = sample;
    io_global = io;

    double fraction_total = (fraction_moving() + fraction_colliding() + fraction_io() + fraction_mpi() + fraction_sample() + fraction_accelerate() + fraction_pressure() + fraction_system_initialize() + fraction_mpi_reduce());
    double time_total = system_initialize + moving + colliding + io + mpi + sample + accelerate + pressure + mpi_reduce;

    double total_time = MPI_Wtime() - t0;
    system.io->time_statistics_file->precision(2);
    *system.io->time_statistics_file << "Program finished after " << (long)total_time << " seconds. Time analysis:" << endl;
    *system.io->time_statistics_file << fixed
         << "      Moving            : " << system.timer->moving << " s ( " << 100*fraction_moving() << "% )" <<  endl
         << "      Colliding         : " << system.timer->colliding << " s ( " << 100*fraction_colliding() << "% )" <<  endl
         << "      Pressure          : " << system.timer->pressure << " s ( " << 100*fraction_pressure() << "% )" <<  endl
         << "      Accelerating      : " << system.timer->accelerate << " s ( " << 100*fraction_accelerate() << "% )" <<  endl
         << "      Sampling          : " << system.timer->sample << " s ( " << 100*fraction_sample() << "% )" <<  endl
         << "      Disk IO           : " << system.timer->io << " s ( " << 100*fraction_io() << "% )" <<  endl
         << "      System initialize : " << system.timer->system_initialize << " s ( " << 100*fraction_system_initialize() << "% )" <<  endl
         << "      MPI communication : " << system.timer->mpi << " s ( " << 100*fraction_mpi() << "% )" <<  endl
         << "      MPI reduce        : " << system.timer->mpi_reduce << " s ( " << 100*fraction_mpi_reduce() << "% )" <<  endl << endl
         << "      TOTAL             : " << time_total << " s ( " << 100*fraction_total << "% )" <<  endl;
    *system.io->time_statistics_file << endl << system.settings->timesteps / total_time << " timesteps / second. " << endl;
    *system.io->time_statistics_file << (double)system.num_molecules_local*system.settings->timesteps / (1000*total_time) << "k atom-timesteps / second. (this processor)" << endl;
}
