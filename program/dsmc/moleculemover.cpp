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

void MoleculeMover::do_move(double *r, double *v, double dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;
}

inline int get_index_of_voxel(double *r, const double &nx_div_lx,const double &ny_div_ly,const double &nz_div_lz, const int &Nx, const int &NyNx) {
    int i =  r[0]*nx_div_lx;
    int j =  r[1]*ny_div_ly;
    int k =  r[2]*nz_div_lz;

    return i + j*Nx + k*NyNx;
}

int sign(double a) {
    return (a<0) ? -1 : 1;
}

void MoleculeMover::move_molecule_box(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;
    
    double *r = &system->r[3*molecule_index];
    double *v = &system->v[3*molecule_index];
    
    do_move(r,v,tau);
    double y = r[1];
    
    double system_center_x = system->length[0]*0.5;
    double system_center_y = system->length[1]*0.5;

    double upper_wall = system_center_y*(1+BOX_FRACTION);
    double lower_wall = system_center_y*(1-BOX_FRACTION);

    int collided = 0;
    if(r[1] >= upper_wall) collided = 1;
    if(r[1] <= lower_wall) collided = 2;
    double dy = 0;

    if(collided > 0) {
        // We collided, move back
        do_move(r,v,-tau);
        // Calculate distance to the wall
        if(collided == 1) {
            // We will collide in the positive direction
            dy = upper_wall - r[1] - 1e-5;
        } else {
            dy = lower_wall - r[1] + 1e-5;
        }

        double time_until_collision = dy / v[1];
        do_move(r,v,time_until_collision);
        dt -= time_until_collision;

        float normal[3];
        normal[0] = 0;
        normal[1] = (collided == 1) ? -1 : 1;
        normal[2] = 0;

        float tangent1[3];
        tangent1[0] = 0;
        tangent1[1] = 0;
        tangent1[2] = 1;

        float tangent2[3];
        
        tangent2[0] = 1;
        tangent2[1] = 0;
        tangent2[2] = 0;

        surface_collider->collide(rnd, v, &normal[0], &tangent1[0], &tangent2[0]);
    } else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule_box(molecule_index,dt,rnd,depth+1);
    }

}

void MoleculeMover::move_molecule_cylinder(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;

    double *r = &system->r[3*molecule_index];
    double *v = &system->v[3*molecule_index];

    do_move(r,v,tau);

    double cylinder_center_x = system->length[0]*0.5;
    double cylinder_center_y = system->length[1]*0.5;

    double dx = r[0] - cylinder_center_x; // Moved origin to center of circle
    double dy = r[1] - cylinder_center_y;
    double dr2 = dx*dx + dy*dy;

    if(dr2>=CYLINDER_RADIUS_SQUARED) {
        // We collided, move back
        do_move(r,v,-tau);

        dx = r[0] - cylinder_center_x;
        dy = r[1] - cylinder_center_y;
        // Time of collision
        double a = v[0]*v[0] + v[1]*v[1];
        double b = 2*(dx*v[0] + dy*v[1]);
        double c = -CYLINDER_RADIUS_SQUARED + dx*dx + dy*dy;
        tau = (-b + sqrt(b*b - 4*a*c))/(2*a);
        if(tau < 0) tau = (-b - sqrt(b*b - 4*a*c))/(2*a);
        tau -= 1e-5;

        do_move(r,v,tau);
        dt -= tau;

        float normal[3];
        normal[0] = cylinder_center_x - r[0];
        normal[1] = cylinder_center_y - r[1];
        normal[2] = 0;
        double norm = sqrt(normal[0]*normal[0] + normal[1]*normal[1]);
        normal[0] /= norm;
        normal[1] /= norm;

        float tangent1[3];
        tangent1[0] = 0;
        tangent1[1] = 0;
        tangent1[2] = 1;

        float tangent2[3];
        // t2 = n x t1
        tangent2[0] = tangent1[1]*normal[2] - tangent1[2]*normal[1];
        tangent2[1] = tangent1[2]*normal[0] - tangent1[0]*normal[2];
        tangent2[2] = tangent1[0]*normal[1] - tangent1[1]*normal[0];

        // Normalize
        norm = sqrt(tangent2[0]*tangent2[0] + tangent2[1]*tangent2[1] + tangent2[2]*tangent2[2]);

        tangent2[0] /= norm;
        tangent2[1] /= norm;
        tangent2[2] /= norm;

        surface_collider->collide(rnd, v, &normal[0], &tangent1[0], &tangent2[0]);
    } else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule_cylinder(molecule_index,dt,rnd,depth+1);
    }
}

void MoleculeMover::apply_periodic_boundary_conditions(int &molecule_index, double *r, const CVector &system_length) {
        if(r[3*molecule_index + 0] > system_length.x)  { r[3*molecule_index + 0] -= system_length.x; count_periodic[0]++; }
        else if(r[3*molecule_index + 0] < 0)         { r[3*molecule_index + 0] += system_length.x; count_periodic[0]--; }

        if(r[3*molecule_index + 1] > system_length.y) { r[3*molecule_index + 1] -= system_length.y; count_periodic[1]++;}
        else if(r[3*molecule_index + 1] < 0)         { r[3*molecule_index + 1] += system_length.y; count_periodic[1]--; }

        if(r[3*molecule_index + 2] > system_length.z) { r[3*molecule_index + 2] -= system_length.z; count_periodic[2]++; }
        else if(r[3*molecule_index + 2] < 0)         { r[3*molecule_index + 2] += system_length.z; count_periodic[2]--; }
}

void MoleculeMover::move_molecule(int &molecule_index, double dt, Random *rnd, int depth) {
//    double nx_div_lx = grid->global_nx*system->one_over_length[0];
//    double ny_div_ly = grid->global_ny*system->one_over_length[1];
//    double nz_div_lz = grid->global_nz*system->one_over_length[2];
//    int nx = grid->Nx;
//    int nynx = grid->Nx*grid->Ny;

    double tau = dt;
    double r_prime[3];
    double *r = &system->r[3*molecule_index];
    double *v = &system->v[3*molecule_index];

    int idx = grid->get_index_of_voxel(r);
    if(voxels[idx] != voxel_type_empty) {
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
            // r_prime[0] = r[0] - v[0]*tau; r_prime[1] = r[1] - v[1]*tau; r_prime[2] = r[2] - v[2]*tau; // Move back, but don't care about periodic boundary conditions
            do_move(r, v, -tau); // Move back
            tau = grid->get_time_until_collision(r, v, collision_voxel_index); // Time until collision with voxel boundary
            do_move(r, v, tau); // Move over there
            idx = grid->get_index_of_voxel(r);
        }

        // We're not at the boundary anymore, so we can move over here and do happy colliding
        dt -= tau;

        if(voxels[collision_voxel_index] != voxel_type_boundary) {
            cout << "Error 2: We hit a wall, but managed to skip the surface voxel. Decrease your timestep." << endl;
            exit(1);
        }

        surface_collider->collide(rnd, v, &grid->normals[3*collision_voxel_index], &grid->tangents1[3*collision_voxel_index], &grid->tangents2[3*collision_voxel_index]);
    } else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule(molecule_index,dt,rnd,depth+1);
    }
}
