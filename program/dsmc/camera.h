#pragma once
 
#include <iostream>
#include <cmath>         // Used only for sin() and cos() functions
 
#include <GL/glfw.h>      // Include OpenGL Framework library for the GLFW_PRESS constant only!
 
#include <cvector.h>


class Camera
{
    protected:
        // Camera rotation
        CVector rotation;
        // Camera movement speed. When we call the move() function on a camera, it moves using these speeds
        CVector speed;
 
        double movement_speed_factor; // Controls how fast the camera moves
        double pitch_sensitivity;    // Controls how sensitive mouse movements affect looking up and down
        double yaw_sensitivity;      // Controls how sensitive mouse movements affect looking left and right
 
        // Window size in pixels and where the midpoint of it falls
        int window_width;
        int window_height;
        int window_mid_x;
        int window_mid_y;
 
        // Method to set some reasonable default values. For internal use by the class only.
        void init_camera(double speed);
 
    public:
        static const double TO_RADS; // The value of 1 degree in radians
 
        // Holding any keys down?
        bool holding_forward;
        bool holding_shift;
        bool holding_backward;
        bool holding_left_strafe;
        bool holding_right_strafe;
        
        CVector target;
        CVector position;
 
        // Constructors
        Camera(float window_width_, float window_height_, double camera_speed);
 
        // Destructor
        ~Camera();
 
        // Mouse movement handler to look around
        void handle_mouse_move(int mouse_x, int mouse_y);
 
        // Method to convert an angle in degress to radians
        double to_rads(const double &angle_in_degrees) const;
 
        // Method to move the camera based on the current direction
        void move();
 
        // --------------------------------- Inline methods ----------------------------------------------
 
        // Setters to allow for change of vertical (pitch) and horizontal (yaw) mouse movement sensitivity
        float get_pitch_sensitivity()            { return pitch_sensitivity;  }
        void  set_pitch_sensitivity(float value) { pitch_sensitivity = value; }
        float get_yaw_sensitivity()              { return yaw_sensitivity;    }
        void  set_yaw_sensitivity(float value)   { yaw_sensitivity   = value; }
 
        // Rotation getters
        CVector get_rotation() const { return rotation;        }
        double get_rot_x()           const { return rotation.x; }
        double get_rot_y()           const { return rotation.y; }
        double get_rot_z()           const { return rotation.z; }
};
