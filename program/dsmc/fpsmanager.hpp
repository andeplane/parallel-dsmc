#pragma once

#include <iostream>
#include <string>
#include <sstream>
 
#ifndef __glfw_h_
    #include <GL/glfw.h> // We need GLFW for this, so let's check for it- although it'd be a doddle to convert to non-GLFW using code.
#endif

using std::string;
using std::cout;
using std::endl;

class FpsManager
{
 
    public:
        double frame_start_time;        // Frame start time
        double frame_end_time;          // Frame end time
        double frame_duration;          // How many milliseconds between the last frame and this frame
 
        double target_fps;              // The desired FPS to run at (i.e. maxFPS)
        double current_fps;             // The current FPS value
        int    frame_count;             // How many frames have been drawn s
 
        double target_frame_duration;   // How many milliseconds each frame should take to hit a target FPS value (i.e. 60fps = 1.0 / 60 = 0.016ms)
        double sleep_duration;          // How long to sleep if we're exceeding the target frame rate duration
 
        double last_report_time;        // The timestamp of when we last reported
        double report_interval;         // How often to update the FPS value
 
        std::string windowTitle;        // Window title to update view GLFW
 
        bool verbose;                   // Whether or not to output FPS details to the console or update the window
        double average_fps;
 
        // Limit the minimum and maximum target FPS value to relatively sane values
        const double MIN_TARGET_FPS = 20.0;
        const double MAX_TARGET_FPS = 60.0; // If you set this above the refresh of your monitor and enable VSync it'll break! Be aware!
 
        // Private method to set relatively sane defaults. Called by constructors before overwriting with more specific values as required.
        void init(double target_fps_)
        {
            set_target_fps(target_fps_);
 
            frame_count     = 0;
 
            current_fps     = 0.0;
            average_fps    = 60.0;
            sleep_duration  = 0.0;
 
            frame_start_time = glfwGetTime();
            frame_end_time   = frame_start_time + 1;
            frame_duration  = 1;
 
            last_report_time = frame_start_time;
            report_interval = 0.2f;
        }
 
    public:
 
        FpsManager()
        {
            init(60.0);
        }

        // Single parameter constructor - just set a desired framerate and let it go.
        // Note: No FPS reporting by default, although you can turn it on or off later with the setVerbose(true/false) method
        FpsManager(int target_fps_)
        {
            init(target_fps_);
        }
 
        // Two parameter constructor which sets a desired framerate and a reporting interval in seconds
        FpsManager(int target_fps_, double report_interval_)
        {
            init(target_fps_);

            set_report_interval(report_interval_);
        }
 
        // Three parameter constructor which sets a desired framerate, how often to report, and the window title to append the FPS to
        FpsManager(int target_fps_, float report_interval_)
        {
            init(target_fps_); // If you specify a window title it's safe to say you want the FPS to update there ;)
 
            set_report_interval(report_interval_);
        }
 
        // Getter and setter for the targetFps property
        int get_target_fps()
        {
            return target_fps;
        }
        void set_target_fps(int target_fps_)
        {
            // Make at least some attempt to sanitise the target FPS...
            if (target_fps_ < MIN_TARGET_FPS)
            {
                target_fps_ = MIN_TARGET_FPS;
                cout << "Limiting FPS rate to legal minimum of " << MIN_TARGET_FPS << " frames per second." << std::endl;
            }
            if (target_fps_ > MAX_TARGET_FPS)
            {
                target_fps_ = MAX_TARGET_FPS;
                cout << "Limiting FPS rate to legal maximum of " << MAX_TARGET_FPS << " frames per second." << std::endl;
            }
 
            // ...then set it and calculate the target duration of each frame at this framerate
            target_fps = target_fps_;
            target_frame_duration = 1.0 / target_fps;
        }
 
        double get_frame_duration() { return frame_duration; } // Returns the time it took to complete the last frame in milliseconds
 
        // Setter for the report interval (how often the FPS is reported) - santises input.
        void set_report_interval(float report_interval_)
        {
            // Ensure the time interval between FPS checks is sane (low cap = 0.1s, high-cap = 10.0s)
            // Negative numbers are invalid, 10 fps checks per second at most, 1 every 10 secs at least.
            if (report_interval_ < 0.1)
            {
                report_interval_ = 0.1;
            }
            if (report_interval_ > 10.0)
            {
                report_interval_ = 10.0;
            }
            report_interval = report_interval_;
        }

        // Method to force our application to stick to a given frame rate and return how long it took to process a frame
        double enforce_fps()
        {
            // Get the current time
            frame_end_time = glfwGetTime();

            // Calculate how long it's been since the frame_start_time was set (at the end of this method)
            frame_duration = frame_end_time - frame_start_time;

            if (report_interval != 0.0f)
            {
                // Calculate and display the FPS every specified time interval
                if ((frame_end_time - last_report_time) > report_interval)
                {
                    // Update the last report time to be now
                    last_report_time = frame_end_time;
 
                    // Calculate the FPS as the number of frames divided by the interval in seconds
                    current_fps =  (double)frame_count / report_interval;
                    average_fps  = 0.8*average_fps + 0.2*current_fps;
 
                    // Reset the frame counter to 1 (and not zero - which would make our FPS values off)
                    frame_count = 1;
                }
                else // FPS calculation time interval hasn't elapsed yet? Simply increment the FPS frame counter
                {
                    ++frame_count;
                }

            } // End of if we specified a report interval section

            // Calculate how long we should sleep for to stick to our target frame rate
            sleep_duration = target_frame_duration - frame_duration;

            // If we're running faster than our target duration, sleep until we catch up!
            if (sleep_duration > 0.0) {
                glfwSleep(target_frame_duration - frame_duration);
            }

            // Reset the frame start time to be now - this means we only need put a single call into the main loop
            frame_start_time = glfwGetTime();

            // Pass back our total frame duration (including any sleep) to be used as our deltaTime value in the game loop (to allow for framerate independant speeds)
            return glfwGetTime() - frame_end_time;
        } // End of our enforceFPS method
};
