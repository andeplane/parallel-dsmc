#include <camera.h>

const double Camera::TO_RADS = 3.141592654 / 180.0; // The value of 1 degree in radians
 
Camera::Camera(float window_width_, float window_height_, double camera_speed)
{
	init_camera(camera_speed);
	window_width = window_width_;
	window_height = window_height_;
	
	// Calculate the middle of the window
	window_mid_x = window_width  / 2.0f;
	window_mid_y = window_height / 2.0f;
 
	glfwSetMousePos(window_mid_x, window_mid_y);
}
 
Camera::~Camera()
{
	// Nothing to do here - we don't need to free memory as all member variables
	// were declared on the stack.
}
 
void Camera::init_camera(double speed)
{
	// How fast we move (higher values mean we move and strafe faster)
	movement_speed_factor = speed;
 
	pitch_sensitivity = 0.2; // How sensitive mouse movements affect looking up and down
	yaw_sensitivity   = 0.2; // How sensitive mouse movements affect looking left and right
 
	// To begin with, we aren't holding down any keys
	holding_shift 		= false;
	holding_forward     = false;
	holding_backward    = false;
	holding_left_strafe  = false;
	holding_right_strafe = false;
}
 
// Function to convert degrees to radians
double Camera::to_rads(const double &angle_in_degrees) const
{
	return angle_in_degrees * TO_RADS;
}
 
// Function to deal with mouse position changes
void Camera::handle_mouse_move(int mouse_x, int mouse_y)
{
	// Calculate our horizontal and vertical mouse movement from middle of the window
	double movement_x = (mouse_x - window_mid_x) * yaw_sensitivity;
	double movement_y  = (mouse_y - window_mid_y) * pitch_sensitivity;
	
	// Apply the mouse movement to our rotation vector. The vertical (look up and down)
	// movement is applied on the X axis, and the horizontal (look left and right)
	// movement is applied on the Y Axis
	rotation.x += movement_y;
	rotation.y += movement_x;
	
	// Limit loking up to vertically up
	if (rotation.x < -90)
	{
		rotation.x = -90;
	}
 
	// Limit looking down to vertically down
	if (rotation.x > 90)
	{
		rotation.x = 90;
	}
 
	// If you prefer to keep the angles in the range -180 to +180 use this code
	// and comment out the 0 to 360 code below.
	//
	// Looking left and right. Keep the angles in the range -180.0f (anticlockwise turn looking behind) to 180.0f (clockwise turn looking behind)
	/*if (rotation.getY() < -180.0f)
	{
	    rotation.addY(360.0f);
	}
 
	if (rotation.getY() > 180.0f)
	{
	    rotation.addY(-360.0f);
	}*/
 
	// Looking left and right - keep angles in the range 0.0 to 360.0
	// 0 degrees is looking directly down the negative Z axis "North", 90 degrees is "East", 180 degrees is "South", 270 degrees is "West"
	// We can also do this so that our 360 degrees goes -180 through +180 and it works the same, but it's probably best to keep our
	// range to 0 through 360 instead of -180 through +180.
	if (rotation.y < 0)
	{
		rotation.y += 360;
	}
	if (rotation.y > 360)
	{
		rotation.y += -360;
	}
 
	// Reset the mouse position to the centre of the window each frame
	glfwSetMousePos(window_mid_x, window_mid_y);
}
 
// Function to calculate which direction we need to move the camera and by what amount
void Camera::move()
{
	// Vector to break up our movement into components along the X, Y and Z axis
	CVector movement;
 
	// Get the sine and cosine of our X and Y axis rotation
	double sinXRot = sin( to_rads( rotation.x ) );
	double cosXRot = cos( to_rads( rotation.x ) );
 
	double sinYRot = sin( to_rads( rotation.y ) );
	double cosYRot = cos( to_rads( rotation.y ) );
 
	double pitchLimitFactor = cosXRot; // This cancels out moving on the Z axis when we're looking up or down
 
	if (holding_forward)
	{
		movement.x += sinYRot * pitchLimitFactor;
		movement.y += -sinXRot;
		movement.z += -cosYRot * pitchLimitFactor;
	}
 
	if (holding_backward)
	{
		movement.x += -sinYRot * pitchLimitFactor;
		movement.y += sinXRot;
		movement.z += cosYRot * pitchLimitFactor;
	}
 
	if (holding_left_strafe)
	{
		movement.x += -cosYRot;
		movement.z += -sinYRot;
	}
 
	if (holding_right_strafe)
	{
		movement.x += cosYRot;
		movement.z += sinYRot;
	}
 
 	float speed = movement_speed_factor*(1+3*holding_shift);
	// Normalise our movement vector
	movement.normalize();
 
	// Finally, apply the movement to our position
	position = position+movement*speed;
	target.x = sinYRot*cosXRot;
	target.y = -sinXRot;
	target.z = -cosYRot*cosXRot;
}
