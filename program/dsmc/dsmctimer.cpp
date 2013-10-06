#include <dsmctimer.h>
#include <mpi.h>
#include <system.h>
#include <topology.h>

DSMCTimer::DSMCTimer() {
    t0 = MPI_Wtime();
    io = 0;
    moving = 0;
    colliding = 0;
    mpi = 0;
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
    colliding_global /= system.topology->num_processors;
    moving_global /= system.topology->num_processors;
    mpi_global /= system.topology->num_processors;
    system_initialize_global /= system.topology->num_processors;
    accelerate_global /= system.topology->num_processors;
    pressure_global /= system.topology->num_processors;
    sample_global /= system.topology->num_processors;
    io_global /= system.topology->num_processors;
}
