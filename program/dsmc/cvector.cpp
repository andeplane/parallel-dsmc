#include <cvector.h>

bool CVector::get_plane_equation(CVector p0, CVector p1, CVector p2, double* eq)
{
	// following two expressions order of arguments is very important 
	// still we need to keep plane's positive normal in plane's equation	
    CVector v0 = p0-p1;
    CVector v1 = p2-p1;
    CVector n;
    n = v0.cross(v1).normalize();
	
    double &A = eq[0], &B = eq[1], &C = eq[2], &D = eq[3];
	
	// positive normal
    A = n.x; 
    B = n.y; 
    C = n.z;
	// Here instead of p0, p1 and p2 can be used as well since they all belong to this plane
    D = - (A*p0.x + B*p0.y + C*p0.z);
    return true;
}



CVector CVector::interpolate(CVector& v1, CVector& v2, CVector& v3, const double& val) {
    // val is between 0 and 1, 0.5 = v2
    double sc1,sc2,sc0;
    if (val<=0.5) {
        sc2 = 0.0;
        sc1 = 2.0*val;
        sc0 = 1.0-sc1;
    }
    else
    {
        double v = val - 0.5;
        sc1 = 1.0 - (2*v);
        sc2 = 1.0 - sc1;
        sc0 = 0.0;
    }

    return v1 * sc0 +v2 * sc1 + v3 * sc2;

}  


