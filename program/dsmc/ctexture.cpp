#include <copengl.h>
#include <ctexture.h>
#include <camera.h>
#include <mpi.h>

CTexture::CTexture(COpenGL *ogl)
{
    opengl = ogl;
    v_absolute_value_average = 0;
    average_velocity.resize(3,1);
}


void CTexture::create_sphere1(string name, int size) {
    CBitMap bmp;
    COpenGLTexture texture;

    bmp.width = size;
    bmp.height= size;
    bmp.data = new unsigned char[4*size*size];
    for (int i=0;i<size; i++) {
        for (int j=0;j<size; j++) {
            double x = (i-size/2)/(double)size;
            double y = (j-size/2)/(double)size;
            double dr = sqrt(x*x + y*y)*2;
            dr = min(dr,1.0);

            if(dr>=0.99) {
                bmp.data[4*i + 4*size*j  +0] = 0;
                bmp.data[4*i + 4*size*j  +1] = 0;
                bmp.data[4*i + 4*size*j  +2] = 0;
                bmp.data[4*i + 4*size*j  +3] = 0;
            } else {
                bmp.data[4*i + 4*size*j  +0] = (unsigned char)(255.0*(1 - 0.7*dr));
                bmp.data[4*i + 4*size*j  +1] = (unsigned char)(255.0*(1 - 0.7*dr));
                bmp.data[4*i + 4*size*j  +2] = (unsigned char)(255.0*(1 - 0.7*dr));
                bmp.data[4*i + 4*size*j  +3] = 255;
            }
        }
    }

    load_texture(&bmp, &texture, true);
    textures.push_back(texture);
    names.push_back(name);
}

void CTexture::load_texture(CBitMap* bmp, COpenGLTexture* texture, bool has_alpha) {
    glGenTextures(1, &texture->id);

    texture->has_alpha = true;
    glBindTexture(GL_TEXTURE_2D, texture->id);

    glTexEnvf( GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE );

    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_NEAREST);

    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR );
    glTexParameterf( GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR );

    if(has_alpha) gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGBA, bmp->width, bmp->height, GL_RGBA, GL_UNSIGNED_BYTE, bmp->data  );
    else gluBuild2DMipmaps( GL_TEXTURE_2D, GL_RGB, bmp->width, bmp->height, GL_RGB, GL_UNSIGNED_BYTE, bmp->data  );
}

void CTexture::render_billboards(double *positions, double *velocities, vector<int> &steps_since_collision, int num_particles, float position_scale) {
    Camera *camera = opengl->camera;
    glEnable(GL_LIGHT0);
    glDisable(GL_LIGHTING);
    glDisable(GL_CULL_FACE);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    // glBlendFunc(GL_ONE, GL_ONE);
    COpenGLTexture texture = textures[0];
    GLuint texture_id = texture.id;

    glEnable( GL_TEXTURE_2D );
    glBindTexture(GL_TEXTURE_2D, texture_id);
    glDepthMask(GL_TRUE);

    glEnable(GL_ALPHA_TEST);
    glAlphaFunc(GL_GREATER,0.9);

    glBegin(GL_QUADS);

    CVector left, up, right, direction, v0, v1, v2, v3;
    double cam_x = camera->position.x; double cam_y = camera->position.y; double cam_z = camera->position.z;

    CVector up_on_screen = opengl->coord_to_ray(0,opengl->window_height/2.0);
    direction = camera->target;

    left = direction.cross(up_on_screen);
    up = (direction.cross(left)).normalize();
    right = (direction.cross(up)).normalize();

    v0 = (right + up);
    v1 = (right*-1 + up);
    v2 = (right*-1 + up*-1);
    v3 = (right + up*-1);
    float one_over_color_cutoff = 1.0/1000;
    float scale = 0.05;

    glNormal3f(direction.x, direction.y, direction.z);
    double average_this_time_step = 0;
    vector<float> new_average_values(3,0);

    for(int index=0; index<num_particles; index++) {
        float x = positions[3*index+0]*position_scale;
        float y = positions[3*index+1]*position_scale;
        float z = positions[3*index+2]*position_scale;
        
        double delta_x = x - cam_x;
        double delta_y = y - cam_y;
        double delta_z = z - cam_z;

        double v_abs = velocities[3*index+0]*velocities[3*index+0] + velocities[3*index+1]*velocities[3*index+1] + velocities[3*index+2]*velocities[3*index+2];
        average_this_time_step += v_abs;

        double dr2 = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z;
        if(dr2 < 0.1) continue;

        double cam_target_times_dr = delta_x*direction.x + delta_y*direction.y + delta_z*direction.z;
        if(cam_target_times_dr < 0) continue;
        double velocity_divided_by_average = v_abs / v_absolute_value_average*0.5;
        double collision_factor = 4.0/steps_since_collision[index];
        float color = 0.3 + 0.7*collision_factor;

        new_average_values[0] += abs(velocities[3*index + 0]);
        new_average_values[1] += abs(velocities[3*index + 1]);
        new_average_values[2] += abs(velocities[3*index + 2]);

        float red = abs(velocities[3*index + 0]) / average_velocity[0] * 0.1;
        float green = abs(velocities[3*index + 1]) / average_velocity[1] * 0.1;
        float blue = abs(velocities[3*index + 2]) / average_velocity[2] * 0.1;
        color = 0.3 + 0.7*velocity_divided_by_average + collision_factor;

        glColor4f(color+red, color+green, color+blue, 1.0);
        // glColor4f(red, green, blue, 1.0);

        glTexCoord2f(0,0);
        glVertex3f(v0.x*scale + x,v0.y*scale + y,v0.z*scale + z);
        glTexCoord2f(1,0);
        glVertex3f(v1.x*scale + x,v1.y*scale + y,v1.z*scale + z);
        glTexCoord2f(1,1);
        glVertex3f(v2.x*scale + x,v2.y*scale + y,v2.z*scale + z);
        glTexCoord2f(0,1);
        glVertex3f(v3.x*scale + x,v3.y*scale + y,v3.z*scale + z);
    }
    average_this_time_step /= num_particles;
    v_absolute_value_average = average_this_time_step;
    average_velocity[0] = new_average_values[0]/num_particles;
    average_velocity[1] = new_average_values[1]/num_particles;
    average_velocity[2] = new_average_values[2]/num_particles;
    new_average_values.clear();

    glEnd();
    glDisable(GL_BLEND);
    // glEnable(GL_CULL_FACE);
    glDisable(GL_ALPHA_TEST);

    glDepthMask(GL_TRUE);
    glDisable( GL_TEXTURE_2D );
    glColor4f(1.0,1.0,1.0,1.0);
}
