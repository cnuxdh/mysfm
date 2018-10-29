
#ifndef ROTATION_HPP
#define ROTATION_HPP

#include"defines.hpp"

// rotation to axis angle vector format
DLL_EXPORT void  rot2aa(double *R, double *aa);

DLL_EXPORT int  eular2rot(double* R, double pitch, double roll, double yaw);


// go from a vector representing a rotation in 
// axis-angle format to a 3x3 rotation matrix 
// in column major form
template<typename T>
void aa2rot(const T* x, T * R)
{
	const T epsilon = T(1e-6);

	T theta = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	
	T wx = x[0]/(theta+epsilon);
	T wy = x[1]/(theta+epsilon);
	T wz = x[2]/(theta+epsilon);

	T costheta = cos(theta);
	T sintheta = sin(theta);

	R[0] = costheta + wx*wx*(T(1.0)-costheta);
	R[1] = wz*sintheta + wx*wy*(T(1.0)-costheta);
	R[2] = -wy*sintheta + wx*wz*(T(1.0)-costheta);
	R[3] = wx*wy*(T(1.0)-costheta) - wz*sintheta;
	R[4] = costheta + wy*wy*(T(1.0)-costheta);
	R[5] = wx*sintheta + wy*wz*(T(1.0)-costheta);
	R[6] = wy*sintheta + wx*wz*(T(1.0)-costheta);
	R[7] = -wx*sintheta + wy*wz*(T(1.0)-costheta);
	R[8] = costheta + wz*wz*(T(1.0)-costheta);
};



//calculate the Eular angle (degree) from rotation matrix
template<typename T>
int  rot2eular(T* R, T* ea)
{
	ea[0] = atan2( R[5], R[8] ) / PI * T(180.0);  //x:pitch
	ea[1] = asin( -R[2] ) / PI * T(180.0);        //y:roll
	ea[2] = atan2( R[1], R[0]) / PI * T(180.0);   //z:yaw

	return 0;
}

/*
//calculate rotation matrix from the Eular angle (degree) 
template<typename T>
int  eular2rot(T* R, T* ea)
{
	T pt = ea[0] / T(180.0) * PI; //pitch
	T ro = ea[1] / T(180.0) * PI; //roll
	T he = ea[2] / T(180.0) * PI; //yaw
	
	ro = -ro;
	pt = -pt;
	
	rr = [cos(ro), 0, sin(ro);
		  0, 1, 0;
		  -sin(ro), 0, cos(ro)];

	rp = [1, 0, 0;
	      0, cos(pt), -sin(pt);
	      0, sin(pt), cos(pt)];
	
    rh = [cos(he), -sin(he), 0;
	      sin(he), cos(he), 0;
	      0, 0, 1];
	
	r = rr * rp * rh;


	return 0;
}
*/





#endif

