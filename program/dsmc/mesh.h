#pragma once
#include <vector>
#include <map>
#include <string>

#ifdef OPENGL
#include <cshaders.h>
#include <testshader.h>
#include <GL/glfw.h>      // Include OpenGL Framework library
#endif
using std::vector;
using std::map;
using std::string;

class CVector;
using std::vector;

class Mesh {
public:
    #ifdef OPENGL
    GLuint vbo_buffers[3];
    TestShader testshader;
    #endif
	vector<float> vertices;
	vector<int>   indices;
	vector<float> normals;
	vector<float> colors;
    bool is_initialized;
	int num_vertices;
    bool is_vbo_built;

    Mesh();
    void verify_initialized();
	void set_vertex(const int &index, float value);
	void set_color(const int &index, float value);
	void set_normal(const int &index, float value);
	void add_vertex(CVector &v);
	void add_color(CVector &v, float alpha);
	void add_normal(CVector &v);
	void generate_smooth_normals(map<string, vector<int> > &vertex_map );
    void initialize(unsigned int num_reserved_vertices);
    #ifdef OPENGL
    void render_triangles();
    void build_vbo();
    void render_vbo();
    void enable_blend(bool inverse);
    void disable_blend();
    #endif
};
