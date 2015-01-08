#include <complexgeometry.h>
#include <cinifile.h>
#include <cvector.h>

using std::cout;
using std::endl;

int main(int argc, char **argv)
{
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
    } else if(type.compare("distancetoatom") == 0) {
        cg.create_from_binary_distancetoatom(ini);
    }

    if(ini.getbool("create_border")) cg.create_border();

    if(ini.getbool("save_file")) cg.save_to_file(ini);

    return 0;
}
