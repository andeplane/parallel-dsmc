#include <cmath.h>

const double CMath::pi = 3.14159265;


  
double* CMath::MatrixZYZ(double phi, double theta, double psi) {

        double  sphi, cphi, sth, cth, spsi, cpsi;

        sphi = sin(phi);
        cphi = cos(phi);

        sth  = sin(theta);
        cth  = cos(theta);

        spsi = sin(psi);
        cpsi = cos(psi);
        
        double* mat = new double[4*4];

        mat[0 + 0*4] = -sphi * spsi + cth * cphi * cpsi;
        mat[1 + 0*4]= -sphi * cpsi - cth * cphi * spsi;
        mat[2 + 0*4] =                sth * cphi;
        mat[3 + 0*4] = 0;

        mat[0 + 1*4] =  cphi * spsi + cth * sphi * cpsi;
        mat[1 + 1*4] =  cphi * cpsi - cth * sphi * spsi;
        mat[2 + 1*4] =                sth * sphi;
        mat[3 + 1*4] = 0;

        mat[0 + 2*4] =              - sth * cpsi;
        mat[1 + 2*4] =                sth * spsi;
        mat[2 + 2*4] =                cth;
        mat[3 + 2*4] = 0;

        mat[0 + 3*4] = 0;
        mat[1 + 3*4] = 0;
        mat[2 + 3*4] = 0;
        mat[3 + 3*4] = 1.0;

        return mat;
}


double CMath::Minmax(double val, const double& min, const double& max) {
   if (val>max) val = max;
   if (val<min) val = min;
   return val;       
} 

double* CMath::Normal2Euler(double x, double y, double z) {
   double phi = acos(y);// - 3.14159/2;
   double theta = atan(z / x) ;//- 3.14159/2;       
   return MatrixZYZ(-phi, -theta, -phi);
}


double CMath::RandomUniform() {
  return (double)rand()/(RAND_MAX+1);
}

double CMath::NormalDistribution(double x, double sigma, double mu) {
  return (1.0/(sigma*2.0*pi))*exp(-pow(x-mu,2)/(2.0*sigma*sigma));
}


double CMath::SmoothStep(double a, double b, double x) {
  if (x<a) return 0;
  if (x>b) return 1;
  x = (x-a)/(double)(b-a);
    return (x*x)*(3-2*x);      
}

double CMath::RandomGauss() {
  double r,val, x =0;
  do
    {
      double sigma=0.2;
      r = (rand()%10000)/(double)10000.0 - 0.5;
      x = (rand()%10000)/(double)1000.0;
      val = 1.0/(sqrt(2*3.14156*sigma*sigma)) * exp(-r*r/(2*sigma*sigma));
      //cout << "r:"<<r<<"     val:" << val <<"    x:"<<x<<endl;
    }
  while (val<x);
  return r + 0.5;
}


/*
Quaternion::Quaternion()
        : x(0.0f), y(0.0f), z(0.0f), w(1.0f)
{

}

Quaternion::~Quaternion()
{

}

void Quaternion::CreateFromAxisAngle(const double &in_x, const double &in_y,
const double &in_z, const double &in_degrees)
{

  double angle = double((in_degrees / 180.0f) * CMath::pi);
        double result = double(sin(angle/2.0f));
        w = double(cos(angle/2.0f));

        // Calculate the x, y and z of the quaternion
        x = double(in_x * result);
        y = double(in_y * result);
        z = double(in_z * result);
}

void Quaternion::CreateMatrix(double *pMatrix)
{
        if(pMatrix)
        {
                // First row
            pMatrix[ 0] = 1.0f - 2.0f * ( y * y + z * z );
            pMatrix[ 1] = 2.0f * ( x * y - w * z );
            pMatrix[ 2] = 2.0f * ( x * z + w * y );
            pMatrix[ 3] = 0.0f;

            // Second row
            pMatrix[ 4] = 2.0f * ( x * y + w * z );
            pMatrix[ 5] = 1.0f - 2.0f * ( x * x + z * z );
            pMatrix[ 6] = 2.0f * ( y * z - w * x );
            pMatrix[ 7] = 0.0f;

            // Third row
            pMatrix[ 8] = 2.0f * ( x * z - w * y );
            pMatrix[ 9] = 2.0f * ( y * z + w * x );
            pMatrix[10] = 1.0f - 2.0f * ( x * x + y * y );
            pMatrix[11] = 0.0f;

            // Fourth row
            pMatrix[12] = 0;
            pMatrix[13] = 0;
            pMatrix[14] = 0;
            pMatrix[15] = 1.0f;
        }
}

Quaternion Quaternion::operator *(const Quaternion &q)
{
        Quaternion r;

        r.w = w*q.w - x*q.x - y*q.y - z*q.z;
        r.x = w*q.x + x*q.w + y*q.z - z*q.y;
        r.y = w*q.y + y*q.w + z*q.x - x*q.z;
        r.z = w*q.z + z*q.w + x*q.y - y*q.x;

        return r;
}
*/
