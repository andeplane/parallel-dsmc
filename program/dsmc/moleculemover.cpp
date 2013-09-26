#include <moleculemover.h>
#include <system.h>
#include <cell.h>
#include <grid.h>
#include <random.h>
#include <grid.h>
#include <settings.h>
#include <colliderbase.h>

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
    for(int n=0;n<system->num_molecules;n++) {
        move_molecule(n,dt,rnd,0);
    }
}

void MoleculeMover::do_move(double *r, double *v, double *r0, const double &dt) {
    r[0] += v[0]*dt;
    r[1] += v[1]*dt;
    r[2] += v[2]*dt;

    if(r[0] > system->length[0])  { r[0] -= system->length[0]; r0[0] -= system->length[0]; count_periodic[0]++; }
    else if(r[0] < 0)         { r[0] += system->length[0]; r0[0] += system->length[0]; count_periodic[0]--; }

    if(r[1] > system->length[1]) { r[1] -= system->length[1]; r0[1] -= system->length[1]; count_periodic[1]++;}
    else if(r[1] < 0)         { r[1] += system->length[1]; r0[1] += system->length[1]; count_periodic[1]--; }

    if(r[2] > system->length[2]) { r[2] -= system->length[2]; r0[2] -= system->length[2]; count_periodic[2]++; }
    else if(r[2] < 0)         { r[2] += system->length[2]; r0[2] += system->length[2]; count_periodic[2]--; }
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

void MoleculeMover::move_molecule_cylinder(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;

    double *r = &system->r[3*molecule_index];
    double *r0 = &system->r0[3*molecule_index];
    double *v = &system->v[3*molecule_index];

    do_move(r,v,r0,tau);

    double cylinder_center_x = system->length[0]*0.5;
    double cylinder_center_y = system->length[1]*0.5;

    double dx = r[0] - cylinder_center_x; // Moved origin to center of circle
    double dy = r[1] - cylinder_center_y;
    double dr2 = dx*dx + dy*dy;

    if(dr2>=CYLINDER_RADIUS_SQUARED) {
        // We collided, move back
        do_move(r,v,r0,-tau);

        dx = r[0] - cylinder_center_x;
        dy = r[1] - cylinder_center_y;
        // Time of collision
        double a = v[0]*v[0] + v[1]*v[1];
        double b = 2*(dx*v[0] + dy*v[1]);
        double c = -CYLINDER_RADIUS_SQUARED + dx*dx + dy*dy;
        tau = (-b + sqrt(b*b - 4*a*c))/(2*a);
        if(tau < 0) tau = (-b - sqrt(b*b - 4*a*c))/(2*a);
        tau *= 0.99;

        do_move(r,v,r0,tau);
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

void MoleculeMover::move_molecule(int &molecule_index, double dt, Random *rnd, int depth) {
    double tau = dt;

    double *r = &system->r[3*molecule_index];
    double *r0 = &system->r0[3*molecule_index];
    double *v = &system->v[3*molecule_index];

    do_move(r,v,r0,dt);

    double nx_div_lx = grid->Nx/system->length[0];
    double ny_div_ly = grid->Ny/system->length[1];
    double nz_div_lz = grid->Nz/system->length[2];
    int nx = grid->Nx;
    int nynx = grid->Nx*grid->Ny;

    int idx = get_index_of_voxel(r,nx_div_lx,ny_div_ly,nz_div_lz,nx,nynx);

    // We have to calculate time until collision
    if(voxels[idx]>=voxel_type_wall) { // Is wall
        if(voxels[idx]!=voxel_type_boundary) { // Not boundary
            int count = 0;
            while(true) {
                idx = get_index_of_voxel(r,nx_div_lx,ny_div_ly,nz_div_lz,nx,nynx);
                if(voxels[idx]==voxel_type_boundary) {
                    dt -= tau;
                    while(*grid->get_voxel(r)>=voxel_type_wall) {
                        dt += 0.1*tau;
                        do_move(r,v,r0,-0.1*tau);
                    }
                    break;
                }

                if(voxels[idx]>=voxel_type_wall) {
                    do_move(r,v,r0,-tau);
                    tau /= 2;
                }
                else {
                    dt -= tau;
                }

                if(++count > 100) {
                    cout << "TROUBLE 1" << endl;
                    exit(0);
                }

                do_move(r,v,r0,tau);
            }
        }
        else {
            int count = 0;
            while(*grid->get_voxel(r)>=voxel_type_wall) {
                dt += 0.1*tau;
                do_move(r,v,r0,-0.1*tau);
                if(++count > 100) {
                    cout << "TROUBLE 2" << endl;
                    exit(0);
                }
            }
        }

        surface_collider->collide(rnd, v, &grid->normals[3*idx], &grid->tangents1[3*idx], &grid->tangents2[3*idx]);
    }
    else dt = 0;

    if(dt > 1e-5 && depth < 10) {
        move_molecule(molecule_index,dt,rnd,depth+1);
    }
}
