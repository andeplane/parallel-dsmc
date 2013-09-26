#include <iostream>
#include "math.h"
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

using namespace std;

int main(int args, char* argv[]) {
    int num_nodes, myid;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &num_nodes);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    double t_start = MPI_Wtime();

    Settings *settings = new Settings("dsmc.ini");
    System system;

    system.initialize(settings, myid);
    
    ifstream to_continue("Tocontinue");
    double t = 0;
    unsigned long steps = 0;
    unsigned long collisions = 0;
    to_continue >> t;
    to_continue >> steps;
    to_continue >> collisions;
    to_continue.close();
    system.t = t;
    system.steps = steps;
    system.collisions = collisions;

    StatisticsSampler *sampler;

    if(myid==0) {
        sampler = new StatisticsSampler(&system);
        for(int i=0;i<settings->timesteps;i++) {
            system.io->save_state_to_movie_file();
            system.step();

            system.timer->start_sample();
            sampler->sample();
            system.timer->end_sample();
        }
        system.io->save_state_to_file_binary();
        sampler->finalize();
        system.io->finalize();
        cout << "Num wall collisions = " << system.mover->surface_collider->num_collisions << endl;

    } else {
        for(int i=0;i<settings->timesteps;i++) {
            system.step();
        }
    }

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

            double fraction_total = fraction_moving + fraction_colliding + fraction_io + fraction_mpi + fraction_sample + fraction_accelerate + fraction_pressure + fraction_system_initialize;
            double time_total = system.timer->system_initialize + system.timer->moving + system.timer->colliding + system.timer->io + system.timer->mpi + system.timer->sample + system.timer->accelerate + system.timer->pressure;

            double total_time = MPI_Wtime() - t_start;
            cout.precision(2);
            cout << endl << "Program finished after " << (long)total_time << " seconds. Time analysis:" << endl;
            cout << fixed
                 << "      Moving            : " << system.timer->moving << " s ( " << 100*fraction_moving << "% )" <<  endl
                 << "      Colliding         : " << system.timer->colliding << " s ( " << 100*fraction_colliding << "% )" <<  endl
                 << "      Pressure          : " << system.timer->pressure << " s ( " << 100*fraction_pressure << "% )" <<  endl
                 << "      Accelerating      : " << system.timer->accelerate << " s ( " << 100*fraction_accelerate << "% )" <<  endl
                 << "      Sampling          : " << system.timer->sample << " s ( " << 100*fraction_sample << "% )" <<  endl
                 << "      Disk IO           : " << system.timer->io << " s ( " << 100*fraction_io << "% )" <<  endl
                 << "      System initialize : " << system.timer->system_initialize << " s ( " << 100*system_initialize_percentage << "% )" <<  endl
                 << "      MPI communication : " << system.timer->mpi << " s ( " << 100*fraction_mpi << "% )" <<  endl << endl
                 << "      TOTAL             : " << time_total << " s ( " << 100*fraction_total << "% )" <<  endl;
            cout << endl << settings->timesteps / total_time << " timesteps / second. " << endl;
            cout << (double)system.num_molecules*settings->timesteps / (1000*total_time) << "k atom-timesteps / second. " << endl;
            cout << (double)system.num_molecules*settings->timesteps / (1000*total_time*num_nodes) << "k atom-timesteps / second (per node). " << endl;
            ofstream to_continue_write("Tocontinue");
            to_continue_write << system.t << " " << system.steps << " " << system.collisions;
            to_continue_write.close();
        }


    MPI_Finalize();

    return 0;
}
