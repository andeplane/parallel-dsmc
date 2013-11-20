#include <iostream>
#include <cmath>
#include "time.h"
#include <fstream>
#include <system.h>
#include <statisticssampler.h>
#include <unitconverter.h>
#include <settings.h>
#include <dsmc_io.h>
#include <mpi.h>
#include <dsmctimer.h>
#include <colliderbase.h>
#include <moleculemover.h>
#define VERSION "1.0.08"
using namespace std;

int main(int args, char* argv[]) {
    int num_nodes, myid;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double t_start = MPI_Wtime();
    if(myid == 0) cout << endl << "DSMC BINARY SOLVER VERSION " << VERSION << endl << endl;

    Settings *settings = new Settings("dsmc.ini");
    if(num_nodes != settings->nx*settings->ny*settings->nz) {
        if(myid==0) cout << "Error, wrong number of CPUs, aborting!" << endl;
        exit(1);
    }

    System system;
    try {
        system.initialize(settings, myid);
    } catch (string exception) {
        cout << "Exception thrown: " << exception << endl;
    }

    ifstream to_continue("Tocontinue");
    double t = 0;
    unsigned long steps = 0;
    unsigned long collisions = 0;
    unsigned long wall_collisions = 0;

    if(settings->load_previous_state) {
        to_continue >> t;
        to_continue >> steps;
        to_continue >> collisions;
        to_continue >> wall_collisions;
    }

    to_continue.close();
    system.t = t;
    system.t0 = t;
    system.steps = steps;
    system.collisions = collisions;
    system.mover->surface_collider->num_collisions = wall_collisions;

    StatisticsSampler *sampler = new StatisticsSampler(&system);

    if(myid==0) {
        for(int i=0;i<settings->timesteps;i++) {
            system.io->save_state_to_movie_file();
            system.step();

            system.timer->start_sample();
            sampler->sample();
            system.timer->end_sample();
        }        
    } else {
        for(int i=0;i<settings->timesteps;i++) {
            system.io->save_state_to_movie_file();
            system.step();

            system.timer->start_sample();
            sampler->sample();
            system.timer->end_sample();
        }
    }
    system.io->save_state_to_file_binary();

    system.timer->start_sample();
    sampler->finalize();
    system.timer->end_sample();

    system.timer->save_to_file(system);
    system.timer->gather_all_nodes(system);

    system.io->finalize();

    if(myid==0) {
            double system_initialize_percentage = system.timer->fraction_system_initialize();
            double fraction_moving = system.timer->fraction_moving();
            double fraction_colliding = system.timer->fraction_colliding();
            double fraction_io = system.timer->fraction_io();
            double fraction_mpi = system.timer->fraction_mpi();
            double fraction_sample = system.timer->fraction_sample();
            double fraction_accelerate = system.timer->fraction_accelerate();
            double fraction_pressure = system.timer->fraction_pressure();
            double fraction_system_initialize = system.timer->fraction_system_initialize();
            double fraction_mpi_reduce = system.timer->fraction_mpi_reduce();

            double fraction_total = (fraction_moving + fraction_colliding + fraction_io + fraction_mpi + fraction_sample + fraction_accelerate + fraction_pressure + fraction_system_initialize + fraction_mpi_reduce);
            double time_total = system.timer->system_initialize + system.timer->moving + system.timer->colliding + system.timer->io + system.timer->mpi + system.timer->sample + system.timer->accelerate + system.timer->pressure + system.timer->mpi_reduce;

            double total_time = MPI_Wtime() - t_start;
            cout.precision(2);
            cout << endl << "Program finished after " << (long)total_time << " seconds. Time analysis:" << endl;
            cout << fixed
                 << "      Moving            : " << system.timer->moving_global << " s ( " << 100*fraction_moving << "% )" <<  endl
                 << "      Colliding         : " << system.timer->colliding_global << " s ( " << 100*fraction_colliding << "% )" <<  endl
                 << "      Pressure          : " << system.timer->pressure_global << " s ( " << 100*fraction_pressure << "% )" <<  endl
                 << "      Accelerating      : " << system.timer->accelerate_global << " s ( " << 100*fraction_accelerate << "% )" <<  endl
                 << "      Sampling          : " << system.timer->sample_global << " s ( " << 100*fraction_sample << "% )" <<  endl
                 << "      Disk IO           : " << system.timer->io_global << " s ( " << 100*fraction_io << "% )" <<  endl
                 << "      System initialize : " << system.timer->system_initialize_global << " s ( " << 100*system_initialize_percentage << "% )" <<  endl
                 << "      MPI communication : " << system.timer->mpi_global << " s ( " << 100*fraction_mpi << "% )" <<  endl
                 << "      MPI reduce        : " << system.timer->mpi_reduce_global << " s ( " << 100*fraction_mpi_reduce << "% )" <<  endl << endl
                 << "      TOTAL             : " << time_total << " s ( " << 100*fraction_total << "% )" <<  endl;
            cout << endl << settings->timesteps / total_time << " timesteps / second. " << endl;
            cout << (double)system.num_molecules_global*settings->timesteps / (1000*total_time) << "k atom-timesteps / second. " << endl;
            cout << (double)system.num_molecules_global*settings->timesteps / (1000*total_time*num_nodes) << "k atom-timesteps / second (per node). " << endl;
            cout << total_time/3600.0*num_nodes << " cpu hours. " << endl;
            ofstream to_continue_write("Tocontinue");
            to_continue_write << system.t << " " << system.steps << " " << system.collisions << " " << system.mover->surface_collider->num_collisions;
            to_continue_write.close();
        }


    MPI_Finalize();

    return 0;
}
