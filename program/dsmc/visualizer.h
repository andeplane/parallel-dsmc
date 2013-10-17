#pragma once
#include <string>
#include <copengl.h>

using std::string;

class Visualizer
{
public:
    bool is_running;
    COpenGL *opengl;
    char window_title[1000];

    Visualizer(int width, int height, string window_title, bool full_screen, double camera_speed);
    void render_begin();
    void render_end();
    void update_window_title();
};
