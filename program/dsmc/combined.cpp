#include <mpi.h>
#include <solver.h>

#include <visualizer.h>
#include <copengl.h>
#include <ctexture.h>
#include <cinifile.h>
#include <marchingcubes.h>
#include <cvector.h>
#include <complexgeometry.h>
#include <moviedata.h>
#include <testshader.h>
#include <copengl.h>
#include <camera.h>
#include <system.h>

using std::cout;
using std::endl;

using namespace std;

int main(int args, char* argv[]) {
    int num_processors, myid;
    MPI_Init(&args,&argv) ;
    MPI_Comm_size(MPI_COMM_WORLD, &num_processors);
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);
    
    CIniFile ini;
    ini.load("geometry_settings.ini");
    int screen_width = ini.getint("screen_width");
    int screen_height = ini.getint("screen_height");
    double Lx = ini.getdouble("Lx");
    double Ly = ini.getdouble("Ly");
    double Lz = ini.getdouble("Lz");
    string type = ini.getstring("type");
    ComplexGeometry cg;

    if(type.compare("box") == 0) {
        cg.create_box(ini);
    } else if(type.compare("poiseuille") == 0) {
        cg.create_poiseuille(ini);
    } else if(type.compare("sphere") == 0) {
        cg.create_sphere(ini);
    } else if(type.compare("perlin") == 0) {
        cg.create_perlin_geometry(ini);
    } else if(type.compare("empty") == 0) {
        cg.create_empty_space(ini);
    } else if(type.compare("diamond_square") == 0) {
        cg.create_diamond_square(ini);
    } else if(type.compare("cylinders") == 0) {
        cg.create_cylinders(ini);
    } else if(type.compare("sinus") == 0) {
        cg.create_sinus(ini);
    } else if(type.compare("random_walk") == 0) {
        cg.create_random_walk(ini);
    } else if(type.compare("packed_spheres") == 0) {
        cg.create_packed_spheres(ini);
    }

    if(ini.getbool("save_file")) cg.save_to_file(ini);

    cout << "Controls: Use WSAD and the mouse to move around!" << endl;
    CVector system_length = CVector(10*Lx, 10*Ly, 10*Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, ini.getdouble("marching_cubes_threshold"), false);

    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
    c.build_vbo();

    CTexture sphere(v.opengl);
    // sphere.create_sphere1("balle",512);
    sphere.load_png("sphere2.png","balle");

    Solver *solver = new Solver(num_processors, myid);
    TestShader shader;
    shader.COpenGLPointer = v.opengl;
    try {
        shader.Initialize("Balle");
    } catch (string ex) {
        cout << ex << endl;
    }

    float t = 0;
    float dt = 0.1;
    float omega = 0.1;

    while(true) {
        v.render_begin();
        CVector balle = CVector(sin(omega*t),cos(omega*t),1.0);
        CVector lightpos = v.opengl->camera->position;
        shader.lightpos = lightpos;
        // shader.targetdir = (v.opengl->camera->target - v.opengl->camera->position).normalize();
        shader.targetdir = v.opengl->camera->target.normalize();

        glBegin(GL_POINTS);
        glVertex3f(lightpos.x, lightpos.y, lightpos.z);
        glEnd();
        
        shader.Start();
        if(v.opengl->bool1) c.render_vbo();
        shader.End();
        
        if(v.opengl->bool2) sphere.render_billboards(solver->system.r, solver->system.v, solver->system.steps_since_collision, solver->system.num_molecules_local, 10.0);
        v.render_end();
        solver->step();
        if(!v.is_running) break;
        v.update_window_title();
        t += dt;
    }

    return 0;
}
