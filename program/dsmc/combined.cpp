#include <defines.h>
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
    
    cout << "Controls: Use WSAD and the mouse to move around!" << endl;
    CIniFile ini;
    ini.load("settings.ini");
    int screen_width = ini.getint("screen_width");
    int screen_height = ini.getint("screen_height");
    double Lx = ini.getdouble("Lx");
    double Ly = ini.getdouble("Ly");
    double Lz = ini.getdouble("Lz");
    double scale = ini.getdouble("scale");
    double threshold = ini.getdouble("threshold");

    ComplexGeometry cg;
    cg.create_box(128,128,128,0.9,true,1);
    // cg.create_sphere(300, 300, 300, 0.9, true, true, 1);
    // cg.save_to_file("perlin.bin");

    scale *= 10; // To map system size to particles
    CVector system_length = CVector(scale*Lx, scale*Ly,scale*Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, threshold*1.1, false);

    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
    c.build_vbo();

    string state_folder = "/projects/master/code/base_code";
    CTexture sphere(v.opengl);
    sphere.create_sphere1("balle",512);

    Solver *solver = new Solver(num_processors, myid);
    TestShader shader;
    shader.COpenGLPointer = v.opengl;
    try {
        shader.Initialize("Balle");
    } catch (string ex) {
        cout << ex << endl;
    }

    CVector lightpos = CVector(1,1,1);//v.opengl->camera->position;
    vector<int> balle;
    balle.resize(100000,0);
    while(true) {
        v.render_begin();
        shader.lightpos = lightpos;
        shader.targetdir = v.opengl->camera->target;
        
        shader.Start();
        if(v.opengl->bool1) c.render_vbo();
        shader.End();
        
        sphere.render_billboards(solver->system.r, solver->system.v, solver->system.steps_since_collision, solver->system.num_molecules_local, scale);
        // sphere.render_billboards(solver->system.r, solver->system.v, balle, solver->system.num_molecules_local, 30.0);
        v.render_end();
        solver->step();
        if(!v.is_running) break;
    }

    return 0;
}
