#include <moleculemover.h>
#include <system.h>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <grid.h>
#include <settings.h>
#include <colliderbase.h>
#include <cvector.h>

MoleculeMover::MoleculeMover()
{

}

double sqrt_wall_temp_over_mass = 0;

void MoleculeMover::initialize(System *system_, ColliderBase *surface_collider_) {
    surface_collider = surface_collider_;
    system = system_;
    voxels = system->world_grid->voxels;
    grid = system->world_grid;
    sqrt_wall_temp_over_mass = sqrt(system->wall_temperature/system->settings->mass);
    count_periodic[0] = 0;
    count_periodic[1] = 0;
    count_periodic[2] = 0;
}

void MoleculeMover::move_molecules(double dt, Random *rnd) {
    for(int n=0;n<system->num_molecules_local;n++) {
        move_molecule(n,dt,rnd,0);
    }
}

int idx = 0;

void MoleculeMover::do_move(data_type *r, data_type *v, double dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;
}

int sign(double a) {
    return (a<0) ? -1 : 1;
}

void MoleculeMover::apply_periodic_boundary_conditions(int &molecule_index, vector<data_type> &r, const CVector &system_length) {
        if(r[3*molecule_index + 0] >= system_length.x)  { r[3*molecule_index + 0] -= system_length.x; count_periodic[0]++; }
        else if(r[3*molecule_index + 0] < 0)         { r[3*molecule_index + 0] += system_length.x; count_periodic[0]--; }

        if(r[3*molecule_index + 1] >= system_length.y) { r[3*molecule_index + 1] -= system_length.y; count_periodic[1]++;}
        else if(r[3*molecule_index + 1] < 0)         { r[3*molecule_index + 1] += system_length.y; count_periodic[1]--; }

        if(r[3*molecule_index + 2] >= system_length.z) { r[3*molecule_index + 2] -= system_length.z; count_periodic[2]++; }
        else if(r[3*molecule_index + 2] < 0)         { r[3*molecule_index + 2] += system_length.z; count_periodic[2]--; }
}

void MoleculeMover::move_molecule(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;
    data_type *r = &system->r.at(3*molecule_index);
    data_type *v = &system->v.at(3*molecule_index);

    int idx = grid->get_index_of_voxel(r);

    if(voxels[idx] != voxel_type_empty && depth==0) {
        cout << "We have fucked up SUPER BIG TIME with molecule " << molecule_index << " at timestep " << system->steps << endl;
        exit(1);
    }

    do_move(r, v, tau);
    idx = grid->get_index_of_voxel(r);

    // We now have three possible outcomes
    if(voxels[idx] >= voxel_type_wall) {
        // We hit a wall. First, move back to find
        while(voxels[idx] != voxel_type_boundary) {
            if(voxels[idx] == voxel_type_wall) {
                tau /= 2;
                do_move(r, v, -tau); // Move back
                idx = grid->get_index_of_voxel(r);
            } else {
                dt -= tau;
                if(dt > 1e-5 && depth < 10) {
                    move_molecule(molecule_index,dt,rnd,depth+1);
                }
                return;
            }
        }

        if(voxels[idx] == voxel_type_empty) {
            cout << "Error 1: We hit a wall, but managed to skip the surface voxel. Decrease your timestep." << endl;
            exit(1);
        }

        int collision_voxel_index = idx;

        while(voxels[idx] == voxel_type_boundary) {
            collision_voxel_index = idx;
            do_move(r, v, -tau); // Move back

            tau = grid->get_time_until_collision(r, v, tau, collision_voxel_index); // Time until collision with voxel boundary
            do_move(r, v, tau); // Move over there
            idx = grid->get_index_of_voxel(r);
        }

        // We're not at the boundary anymore, so we can move over here and do happy colliding
        dt -= tau;

        if(voxels[collision_voxel_index] != voxel_type_boundary) {
            cout << "Error 2: We hit a wall, but managed to skip the surface voxel. Decrease your timestep." << endl;
            exit(1);
        }

        system->steps_since_collision[molecule_index] = 0;
        surface_collider->collide(rnd, v, &grid->normals[3*collision_voxel_index], &grid->tangents1[3*collision_voxel_index], &grid->tangents2[3*collision_voxel_index]);
    } else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule(molecule_index,dt,rnd,depth+1);
    }
}

