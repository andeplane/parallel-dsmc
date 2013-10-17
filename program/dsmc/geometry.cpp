#include <defines.h>

#ifdef OPENGL
#include <visualizer.h>
#include <copengl.h>
#include <ctexture.h>
#include <complexgeometry.h>
#include <testshader.h>
#include <camera.h>
#endif
#include <cinifile.h>
#include <marchingcubes.h>
#include <cvector.h>

using std::cout;
using std::endl;

int main(int argc, char **argv)
{
    cout << "Controls: Use WSAD and the mouse to move around!" << endl;
    CIniFile ini;
    ini.load("settings.ini");
    int screen_width = ini.getint("screen_width");
    int screen_height = ini.getint("screen_height");
    double Lx = ini.getdouble("Lx");
    double Ly = ini.getdouble("Ly");
    double Lz = ini.getdouble("Lz");
    double threshold = ini.getdouble("threshold");

    ComplexGeometry cg;
    // cg.create_perlin_geometry(100, 100, 100, 1,1,1,3, threshold, false, 1);
    cg.create_sphere(100, 100, 100, 0.8, false, true, 1);
    cg.save_to_file("./sphere_world",CVector(2,2,2));

    #ifdef OPENGL
    CVector system_length = CVector(10*Lx, 10*Ly, 10*Lz);
    MarchingCubes c;
    c.create_marching_cubes_from_complex_geometry(cg, system_length, threshold, false);
    char *window_title = new char[1000];
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]", 60.0);
    Visualizer v(screen_width, screen_height, string(window_title), false, 0.1);
    c.build_vbo();

    TestShader shader;
    shader.COpenGLPointer = v.opengl;
    try {
        shader.Initialize("Balle");
    } catch (string ex) {
        cout << ex << endl;
    }

    CVector lightpos = CVector(1,1,1);//v.opengl->camera->position;

    while(true) {
         v.render_begin();
         shader.lightpos = lightpos;
         shader.targetdir = v.opengl->camera->target;

         shader.Start();
         if(v.opengl->bool1) c.render_vbo();
         shader.End();
         v.render_end();
        if(!v.is_running) break;
    }
    #endif

    return 0;
}
