#pragma once

#include <vector>
using namespace std;

class CMatrix
{
 public:  
  double** M;
    CMatrix() 
    {
      M = 0;
      M = new double*[3];
      for (int i=0;i<3;i++)
	M[i] = new double[3];

    }

    void free() {
      if (M) {
	for (int i=0;i<3;i++)
	  delete[] M[i];
	delete[] M;
	
      } 
    }
    void Mul(CMatrix o)
    {
	CMatrix temp;
	int i,j;
	for (j=0;j<3;j++)
		for (i=0;i<3;i++)
		  temp.M[i][j] = M[i][0] * o.M[0][j] + 
		                 M[i][1] * o.M[1][j] +
		                 M[i][2] * o.M[2][j] ;

	Copy(temp);
	temp.free();
    }


    void Zyz(double phi, double theta, double psi) {
         
        double  sphi, cphi, sth, cth, spsi, cpsi;

        sphi = sin(phi);
        cphi = cos(phi);

        sth  = sin(theta);
        cth  = cos(theta);

        spsi = sin(psi);
        cpsi = cos(psi);
        

        M[0][0] = -sphi * spsi + cth * cphi * cpsi;
        M[1][0]= -sphi * cpsi - cth * cphi * spsi;
        M[2][0] =                sth * cphi;

        M[0][1] =  cphi * spsi + cth * sphi * cpsi;
        M[1][1] =  cphi * cpsi - cth * sphi * spsi;
        M[2][1] =                sth * sphi;

        M[0][2] =              - sth * cpsi;
        M[1][2] =                sth * spsi;
        M[2][2] =                cth;

        }
      

    void Copy(CMatrix src)
    {
	int i,j;
	for (i=0;i<3;i++)
		for (j=0;j<3;j++)
                  M[i][j] = src.M[i][j];
    }
    
    void Identity()
    {
	M[0][0] = 1; M[0][1] = 0; M[0][2] = 0;  
	M[1][0] = 0; M[1][1] = 1; M[1][2] = 0;  
	M[2][0] = 0; M[2][1] = 0; M[2][2] = 1; 
    }

    void RotateXY(double angle)
    {
	double c, s;
	c = cos(angle);
	s = sin(angle);
	M[0][0] = c; M[0][1] = s;  M[0][2] = 0; 
	M[1][0] =-s; M[1][1] = c;  M[1][2] = 0; 
	M[2][0] = 0; M[2][1] = 0;  M[2][2] = 1; 
    }

    void RotateYZ(double angle)
    {
	double c, s;
	c = cos(angle);
	s = sin(angle);
	M[0][0] = 1; M[0][1] = 0;  M[0][2] = 0; 
	M[1][0] = 0; M[1][1] = c;  M[1][2] = s; 
	M[2][0] = 0; M[2][1] =-s;  M[2][2] = c; 
    }
    
    void RotateXZ(double angle)
    {
	double c, s;
	c = cos(angle);
	s = sin(angle);
	M[0][0] = c; M[0][1] = 0;  M[0][2] =-s; 
	M[1][0] = 0; M[1][1] = 1;  M[1][2] = 0; 
	M[2][0] = s; M[2][1] = 0;  M[2][2] = c; 
    }

};


