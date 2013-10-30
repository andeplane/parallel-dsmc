#include <visualizer.h>
#include <copengl.h>
#include <camera.h>

COpenGL *ogl;

// Function to draw our scene
void Visualizer::render_begin()
{
    is_running = !glfwGetKey(GLFW_KEY_ESC) && glfwGetWindowParam(GLFW_OPENED);
    opengl->fps_manager.enforce_fps();
    opengl->camera->move();

    // Clear the screen and depth buffer
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    // float c = 179/255.0;
    // glClearColor(c,c,c, 0.0f);

    // Reset the matrix
    glMatrixMode(GL_MODELVIEW);
    glDrawBuffer(GL_BACK_RIGHT);
    glLoadIdentity();

    // Move the camera to our location in space
    glRotatef(opengl->camera->get_rot_x(), 1.0f, 0.0f, 0.0f); // Rotate our camera on the x-axis (looking up and down)
    glRotatef(opengl->camera->get_rot_y(), 0.0f, 1.0f, 0.0f); // Rotate our camera on the  y-axis (looking left and right)

    // Translate the ModelView matrix to the position of our camera - everything should now be drawn relative
    // to this position!
    glTranslatef( -opengl->camera->position.x, -opengl->camera->position.y, -opengl->camera->position.z );
    opengl->set_standard_light();
}

void Visualizer::render_end() {
    // ----- Stop Drawing Stuff! ------
    glfwSwapBuffers(); // Swap the buffers to display the scene (so we don't have to watch it being drawn!)
}


// Callback function to handle mouse movements
void handle_mouse_move(int mouse_x, int mouse_y) {
    ogl->camera->handle_mouse_move(mouse_x, mouse_y);
}

// Callback function to handle keypresses
void handle_keypress(int theKey, int theAction) {
    // If a key is pressed, toggle the relevant key-press flag
    if (theAction == GLFW_PRESS)
    {
        switch (theKey)
        {
        case 287:
            ogl->camera->holding_shift = true;
            break;
        case 'W':
            ogl->camera->holding_forward = true;
            break;
        case 'S':
            ogl->camera->holding_backward = true;
            break;
        case 'A':
            ogl->camera->holding_left_strafe = true;
            break;
        case 'D':
            ogl->camera->holding_right_strafe = true;
            break;
        case '1':
            ogl->bool1 = !ogl->bool1;
            break;
        case '2':
            ogl->bool2 = !ogl->bool2;
            break;
        case '3':
            ogl->bool3 = !ogl->bool3;
            break;
        case '4':
            ogl->bool4 = !ogl->bool4;
            break;
        case '5':
            ogl->bool5 = !ogl->bool5;
            break;
        }
    }
    else // If a key is released, toggle the relevant key-release flag
    {
        switch (theKey)
        {
        case 287:
            ogl->camera->holding_shift = false;
            break;
        case 'W':
            ogl->camera->holding_forward = false;
            break;
        case 'S':
            ogl->camera->holding_backward = false;
            break;
        case 'A':
            ogl->camera->holding_left_strafe = false;
            break;
        case 'D':
            ogl->camera->holding_right_strafe = false;
            break;
        }
    }
}

Visualizer::Visualizer(int width, int height, string window_title, bool full_screen, double camera_speed)
{
    is_running = true;
    ogl = new COpenGL();
    ogl->initialize(width, height, window_title, handle_keypress, handle_mouse_move, full_screen, camera_speed);
    opengl = ogl;
    GLenum error = glewInit();
}

void Visualizer::update_window_title() {
    double fps = opengl->fps_manager.average_fps;
    // Calculate the current time in pico seconds to show in the title bar
    sprintf(window_title, "DSMC Geometry Visualizer (DSMCGV) - [%.2f fps]",fps);
    opengl->set_window_title(string(window_title));
}
