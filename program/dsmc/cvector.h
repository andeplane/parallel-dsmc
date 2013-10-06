#pragma once

using namespace std;

#include <math.h>
#include <iostream>
#include <vector>
#include <cutil.h>

class CVector  {
 public:
  double x,y,z;

  // Constructors
  inline CVector() {
    x=y=z=0;         
  }

  inline CVector(const double& px, const double& py, const double& pz) {
    x = px;
    y = py;
    z = pz;
  }

  inline CVector(const CVector& v)
  {
    x=v.x;
    y=v.y;
    z=v.z;
  }


  // Operators 

  inline CVector operator+(const CVector& a) const {
     return CVector(a.x +x, a.y+y, a.z+z);
   }

  inline CVector operator-(const CVector& a) const {
     return CVector(x - a.x, y - a.y, z - a.z);
   }

  inline void subtract(const CVector& a)   {
    x-=a.x;
    y-=a.y;
    z-=a.z;
  }
  inline void add(const CVector& a)   {
    x+=a.x;
    y+=a.y;
    z+=a.z;
  }


  inline double operator[](const int i) const {
    if (i==0) return x;
    if (i==1) return y;
    return z;
  }
 
  inline CVector operator*(double scale) const {
     return CVector(x * scale, y * scale, z * scale);
   }

  inline CVector operator/(const double scale) const {
     return CVector(x / scale, y / scale, z / scale);
   }
 


  inline CVector operator/(const CVector& o) const {
     return CVector(x / o.x, y / o.y, z / o.z);
   }

  friend ostream& operator<<(ostream& os, const CVector& o) {
    os <<"[ " << o.x << " " << o.y << " " << o.z << " ";
    os << "]" << endl;
    return os;
  }

  inline void operator=(const double& v) {
     x = v; y= v; z=v;
  }
  inline bool operator==(const CVector& v) const {
       if (x==v.x && y==v.y && z==v.z) {
          return true;
       }

       return false;     
  }

  inline void operator=(const int& v) {
     x = v; y= v; z=v;
  }

  // doubleransformations

  inline void inverse() {
     x=-x;y=-y;z=-z;
   }

  inline CVector mul(const CVector& v) const {
     return CVector(x * v.x, y * v.y, z * v.z);
   }
 
  inline CVector rotate_2d(double t) const {
      return CVector(cos(t)*x - sin(t)*y, sin(t)*x+cos(t)*y,0);
  }

  inline CVector rotate_y(double t) const {
   return   CVector(x * cos(t) - z * sin(t) ,y, x*sin(t) + z*cos(t));
 }

  inline CVector rotate_x(double t)  const {
   return CVector(x, y * cos(t) - z * sin(t), y*sin(t) + z*cos(t));
 }

  inline CVector rotate_z(double t) const {
   return CVector(x * cos(t) - y * sin(t), x*sin(t) + y*cos(t),z);
  }
  
  inline double length() const {
    return sqrt(x*x + y*y + z*z);
  }
  
  inline double length_2() const {
    return (x*x + y*y + z*z);
  }
  
  inline CVector normalize() const {
    double l = length();
    if (l!=0) {
        return *this/l;
    }

    return *this;
  }
  
  inline double dot( const CVector& o) const {
    return (double)(x*o.x + y*o.y + z*o.z);
  }
  
  
  inline CVector cross(const CVector& o) const {
    return CVector(o.y*z - o.z*y, o.z*x - o.x*z, o.x*y - o.y*x);
  }
  
  // Return statements
  inline void min_max(CVector& min, CVector& max) {
    if (x>max.x) max.x=x;
    if (y>max.y) max.y=y;
    if (z>max.z) max.z=z;
    
    if (x<min.x) min.x=x;
    if (y<min.y) min.y=y;
    if (z<min.z) min.z=z;
  }

  inline CVector xz() const {
    return CVector(x,0,z);
  }
  
  inline CVector xy() const {
    return CVector(x,y,0);
  }
  
  inline CVector yz() const {
    return CVector(0,y,z);
  }
  
  // Utilities
  
  inline void set(const double px, const double py, const double pz) {
     x = px; y=py; z = pz;
   }

  inline void set(const CVector& v) {
     x = v.x; y=v.y; z = v.z;
   }

  inline void to_double(double* a) {
    if (a==0) {
        throw string("CVector::todouble error: array not allocated");
    }
    a[0] = x;
    a[1] = y;
    a[2] = z;
  }

  inline void to_float(float* a) {
    if (a==0) {
        throw string("CVector::todouble error: array not allocated");
    }
    a[0] = x;
    a[1] = y;
    a[2] = z;
  }

  inline double distance_from_plane(CVector& plane_normal, CVector& V) const {
       return plane_normal.dot(*this - V);
  }

 inline CVector from_spherical() const {
   return CVector( cos(y)*sin(x), sin(y)*sin(x), cos(x)).normalize() * z;
 }

  static bool get_plane_equation(CVector p0, CVector p1, CVector p2, double* eq);
  static CVector interpolate(CVector& v1, CVector& v2, CVector& v3, const double& val);
  double distance_from_plane(CVector& plane_normal, CVector& V);
};
