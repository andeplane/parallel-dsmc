#include <copengl.h>
#include <cmath.h>
#include <cvector.h>
#include <camera.h>       // Include our Camera header so we can work with Camera objects
#include <cmath>

// Frame counter and window settings variables
int redBits    = 8, greenBits = 8,    blueBits    = 8;
int alphaBits  = 8, depthBits = 24,   stencilBits = 0;
int stereoBits = 2;

void COpenGL::initialize(int w, int h, string window_title_, GLFWkeyfun handle_keypress, GLFWmouseposfun handle_mouse_move, bool full_screen, double camera_speed) {
    // glfwWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
    // glfwWindowHint(GLFW_OPENGL_VERSION_MINOR, 2);
    // glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    // glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    window_title = window_title_;

    window_width = w;
    window_height = h;
    mid_window_x = window_width/2;
    mid_window_y = window_height/2;
    aspect_ratio = (GLfloat)window_width/GLfloat(window_height);

    field_of_view = 60.0f;            // Define our field of view (i.e. how quickly foreshortening occurs)
    near        = 0.01f;             // The near (Z Axis) point of our viewing frustum (default 2.0f)
    far         = 1500.0f;          // The far  (Z Axis) point of our viewing frustum (default 1500.0f)
    bool1 = true;
    bool2 = true;
    bool3 = true;
    bool4 = true;
    bool5 = true;
    // Initialise GLFW
    if (!glfwInit() )
    {
        std::cout << "Failed to initialise GLFW!" << endl;
        glfwTerminate();
        exit(1);
    }

    // Create a window
    // GLFW_FULLSCREEN / GLFW_WINDOW
    if(full_screen) {
        GLFWvidmode desktop;
        glfwGetDesktopMode( &desktop );
        window_width = desktop.Width;
        window_height = desktop.Height;
        redBits = desktop.RedBits;
        greenBits = desktop.GreenBits;
        blueBits = desktop.BlueBits;

        mid_window_x = window_width/2;
        mid_window_y = window_height/2;
        aspect_ratio = (GLfloat)window_width/GLfloat(window_height);

        if( !glfwOpenWindow(window_width, window_height, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_FULLSCREEN))
        {
            std::cout << "Failed to open window!" << std::endl;
            glfwTerminate();
            exit(1);
        }

        glfwSetWindowPos(0,0);
    } else {
        if( !glfwOpenWindow(window_width, window_height, redBits, greenBits, blueBits, alphaBits, depthBits, stencilBits, GLFW_WINDOW))
        {
            std::cout << "Failed to open window!" << std::endl;
            glfwTerminate();
            exit(1);
        }

        glfwSetWindowPos(0,0);
    }

    // Call our initGL function to set up our OpenGL options
    init_GL();

    cout << "Supported shading language: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << endl;
    cout << "OpenGL version: " << glGetString(GL_VERSION) << endl;

    // Instantiate our pointer to a Camera object providing it the size of the window
    camera = new Camera(window_width, window_height, camera_speed);

    // Set the mouse cursor to the centre of our window
    glfwSetMousePos(mid_window_x, mid_window_y);

    // // Specify the function which should execute when a key is pressed or released
    glfwSetKeyCallback(handle_keypress);

    // // Specify the function which should execute when the mouse is moved
    glfwSetMousePosCallback(handle_mouse_move);
}

void COpenGL::set_standard_light() {
    return;
    float a = 100;
    GLfloat ambient[4] = {0.05f, 0.05f, 0.05f, 1.0f};
    GLfloat diffuse[4] = {0.5f, 0.5f, 0.5f, 1.0f};
    // GLfloat position[4] = {a*cos(t), a*sin(t*0.342), a*sin(t*0.442), 1.0f};
    GLfloat position[4] = {camera->position.x, camera->position.y, camera->position.z, 1.0f};
    GLfloat shininess[1] = {10.0f};
    GLfloat specular[4] = {0.2f, 0.2f, 0.2f, 1.0f};

    glDisable(GL_LIGHTING);

    glEnable(GL_LIGHTING);
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position);
    glEnable(GL_LIGHT0);
    glEnable(GL_COLOR_MATERIAL);
    glColor4f(1.0, 1.0, 0.0, 1.0);
/*
    glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, front_mat_shininess);
    glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, front_mat_specular);
    glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, front_mat_diffuse);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_COLOR_MATERIAL);*/
}

void COpenGL::set_window_title(string window_title_) {
  window_title = window_title_;
  glfwSetWindowTitle(window_title.c_str());
  glfwPollEvents();
}

void COpenGL::init_GL() {
    // ----- GLFW Settings -----
    
    glfwDisable(GLFW_MOUSE_CURSOR); // Hide the mouse cursor
 
    glfwSwapInterval(0);            // Disable vsync
 
    // ----- Window and Projection Settings -----
 
    // Set the window title
    glfwSetWindowTitle(window_title.c_str());
 
    // Setup our viewport to be the entire size of the window
    glViewport(0, 0, (GLsizei)window_width, (GLsizei)window_height);
 
    // Change to the projection matrix, reset the matrix and set up our projection
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
 
    // The following code is a fancy bit of math that is eqivilant to calling:
    // gluPerspective(fieldOfView / 2.0f, width / height, near, far);
    // We do it this way simply to avoid requiring glu.h
    // GLfloat aspectRatio = (window_width > window_height)? float(window_width)/float(window_height) : float(window_height)/float(window_width);
    GLfloat fH = tan( float(field_of_view / 360.0f * 3.14159f) ) * near;
    GLfloat fW = fH * aspect_ratio;
    glFrustum(-fW, fW, -fH, fH, near, far);
 
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
 
    // ----- OpenGL settings -----
 
    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);              // Set out clear colour to black, full alpha
    glEnable(GL_DEPTH_TEST);                           // Enable the depth buffer
    glClearDepth(1.0f);                                // Clear the entire depth of the depth buffer
    glDepthFunc(GL_LEQUAL);                            // Set our depth function to overwrite if new value less than or equal to current value
    glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST); // Ask for nicest perspective correction
    glLineWidth(2.0f);                                 // Set a 'chunky' line width
    glShadeModel(GL_SMOOTH);
}

void COpenGL::push() {
  glPushMatrix();
}

void COpenGL::pop() {
  glPopMatrix();
}

CVector COpenGL::coord_to_ray(double px, double py) {
   double P = field_of_view / 360.0 * 2*3.14159; // convert to radians
   float pm[16];     // to get viewport matrix
   CVector res, dir;
   
   double x = px + window_width/2;
   double y = py + window_height/2;
   // modifiers 
   double mmx = tan(P*0.5)*(1.0-(x*2)/(double)window_width)* (window_width/(double)window_height);
   double mmy = tan(P*0.5)*-(1.0-(y*2)/(double)window_height);
   
   glGetFloatv(GL_MODELVIEW_MATRIX, pm);
   // find position in viewspace
   dir = CVector(mmx,mmy,1);
   res = (dir.gl_mat_mul(pm)).normalize();
   
   return res;       
}

// switch to a texture class.. much better
void COpenGL::buffer2texture(GLuint texture, int w, int h, int type) {

  glEnable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0);
  glBindTexture(GL_TEXTURE_2D, texture);
  if (type==MIPMAP)
    type = MIPMAPALPHA;

  if (type==NOMIPMAP) {
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,0,0,w,h,0);
  }
  if (type==MIPMAPALPHA) {
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);
    glCopyTexImage2D(GL_TEXTURE_2D,0,GL_RGBA,0,0,w,h,0);
    unsigned char *t = new unsigned char[w*h*4];

    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );
    // glTexSubImage2D(GL_TEXTURE_2D,0,0,0,w,h,GL_RGBA, GL_UNSIGNED_BYTE,t);
    glGetTexImage(GL_TEXTURE_2D,0,GL_RGBA, GL_UNSIGNED_BYTE,t);
    //  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
    gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA ,w, h, GL_RGBA, GL_UNSIGNED_BYTE, t  );
    delete[] t;
  }

   glDisable(GL_TEXTURE_2D);

}

void COpenGL::SetOrthographicProjection() {
  glMatrixMode(GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();
  gluOrtho2D(0, window_width, 0, window_height);
  glScalef(1, -1, 1);
  glTranslatef(0, -window_height, 0);
  glMatrixMode(GL_MODELVIEW);
}

void COpenGL::ResetPerspectiveProjection() {
  glMatrixMode(GL_PROJECTION);
  glPopMatrix();
  glMatrixMode(GL_MODELVIEW);
}
