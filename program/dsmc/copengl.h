#pragma once
#include <GL/glew.h>
#include <GL/glfw.h>      // Include OpenGL Framework library
#include <oglshader.h>
#include <string>
#include <fpsmanager.hpp> // Include our FpsManager class

using std::string;

class CShaderContainer;
class FpsManager;
class Camera;
class CVector;

class COpenGL {
public:   
    static const int MIPMAP = 0;
    static const int NOMIPMAP = 1;
    static const int MIPMAPALPHA = 2;

    CShaderContainer shadercontainer;

    GLint window_width, window_height; 
    GLint mid_window_x, mid_window_y;
    double aspect_ratio;
    bool full_screen;
    
    Camera *camera;
    string window_title;

    // Create a FPS manager that locks to 60fps and updates the window title with stats every 3 seconds
    FpsManager fps_manager;
    GLfloat perspective;
    GLfloat field_of_view;
    GLfloat near;
    GLfloat far;
    bool bool1;
    bool bool2;
    bool bool3;
    bool bool4;
    bool bool5;

    void initialize(int w, int h, string window_title_, GLFWkeyfun cbfun, GLFWmouseposfun, bool full_screen, double camera_speed);
    void pop();
    void push();
    void init_GL();
    void set_window_title(string title);
    void set_standard_light();
    // Fra nicolaas
    void buffer2texture(GLuint texture, int w, int h, int mipmap);
    void SetOrthographicProjection();
    void ResetPerspectiveProjection();
    CVector coord_to_ray(double px, double py);
    // End fra nicolaas

    COpenGL() { }
}; 
