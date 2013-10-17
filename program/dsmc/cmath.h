#pragma once

using namespace std;


#include <math.h>
#include <CMatrix.h>
#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

class CMath {
 public:
             
  static double* MatrixZYZ(double phi, double theta, double psi);
  static double Minmax(double, const double&, const double&);    
  static double* Normal2Euler(double x, double y, double z);

  const static double pi;
  //  const static double pi = 3.14159265;


  static double RandomUniform();
  static double NormalDistribution(double x, double sigma, double mu);
  static double SmoothStep(double a, double b, double x);
  static double RandomGauss();

};

/*
class Quaternion
{
 public:
  
  Quaternion();
  ~Quaternion();
  
  void CreateMatrix(double *pMatrix);
  void CreateFromAxisAngle(const float &in_x,
			   const float &in_y,
			   const float &in_z,
			   const float &in_degrees);
  
  Quaternion operator *(const Quaternion &q);
  float x,y,z,w;
};
*/
