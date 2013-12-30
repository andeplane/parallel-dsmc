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

inline void MoleculeMover::do_move(data_type *r, data_type *v, double dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;
}

int sign(double a) {
    return (a<0) ? -1 : 1;
}

int i_t = 2102;
int i_m_index = 9159;

void MoleculeMover::move_molecule(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;
    data_type *r = &system->r[3*molecule_index];
    data_type *v = &system->v[3*molecule_index];

    int idx = grid->get_index_of_voxel(r);

    if(voxels[idx] != voxel_type_empty) {
        cout << "Warning, an invalid situation occured, maybe too large timestep?" << endl;
        system->find_position(molecule_index);
        return;
    }

    do_move(r, v, tau);
    idx = grid->get_index_of_voxel(r);
    int count = 0;
    // We now have three possible outcomes
    if(voxels[idx] >= voxel_type_wall) {
        // We hit a wall. First, move back to find boundry
        while(voxels[idx] != voxel_type_boundary) {
            if(++count > 100) {
                system->find_position(molecule_index);
                cout << "Warning, we found an infinite loop with molecule " << molecule_index << " at timestep " << system->steps << "on a wall, index " << idx << endl;
                return;
            }

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
        count = 0;
        while(voxels[idx] == voxel_type_boundary) {
            if(++count > 100) {
                system->find_position(molecule_index);
                cout << "Warning, we found an infinite loop with molecule " << molecule_index << " at timestep " << system->steps << "on a surface, index " << idx << endl;
                return;
            }
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
        bool print_details = i_m_index == molecule_index && i_t == system->steps && false;
        surface_collider->collide(rnd, v, &grid->normals[3*collision_voxel_index], &grid->tangents1[3*collision_voxel_index], &grid->tangents2[3*collision_voxel_index], print_details);
    } else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule(molecule_index,dt,rnd,depth+1);
    }
}

void MoleculeMover::move_molecule_cylinder(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;
    data_type *r = &system->r[3*molecule_index];
    data_type *v = &system->v[3*molecule_index];
    double x0 = r[0] - system->length[0]/2.0;
    double y0 = r[1] - system->length[1]/2.0;

    do_move(r,v,tau);

    double dx = r[0] - system->length[0]/2.0;
    double dy = r[1] - system->length[1]/2.0;
    double dr2 = dx*dx + dy*dy;
    if(dr2 > CYLINDER_RADIUS_SQUARED) {
        double a = v[0]*v[0] + v[1]*v[1];
        double b = 2*x0*v[0] + 2*y0*v[1];
        double c = x0*x0 + y0*y0 - CYLINDER_RADIUS_SQUARED;
        double t0 = (-b + sqrt(b*b - 4*a*c)) / (2*a);
        double t1 = (-b - sqrt(b*b - 4*a*c)) / (2*a);
        double correct_time = max(t0,t1);
        double time_to_move_back = tau - correct_time;
        do_move(r,v,-time_to_move_back);
        float *normals = new float[3];
        float *tangents1 = new float[3];
        float *tangents2 = new float[3];

        normals[0] = -(r[0] - system->length[0]/2.0);
        normals[1] = -(r[1] - system->length[1]/2.0);
        normals[2] = 0;
        double norm = sqrt(normals[0]*normals[0] + normals[1]*normals[1]);
        normals[0] /= norm;
        normals[1] /= norm;

        tangents1[0] = 0; tangents1[1] = 0; tangents1[2] = 1;

        tangents2[0] = tangents1[1]*normals[2] - tangents1[2]*normals[1];
        tangents2[1] = tangents1[2]*normals[0] - tangents1[0]*normals[2];
        tangents2[2] = tangents1[0]*normals[1] - tangents1[1]*normals[0];

        // Normalize
        norm = sqrt(tangents2[0]*tangents2[0] + tangents2[1]*tangents2[1] + tangents2[2]*tangents2[2]);

        tangents2[0] /= norm; tangents2[1] /= norm; tangents2[2] /= norm;

        surface_collider->collide(rnd, v, normals, tangents1, tangents2, false);
        move_molecule_cylinder(molecule_index,time_to_move_back,rnd,depth+1);
    }
}
