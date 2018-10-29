
#include "math.h"

#include"matrix.h"

#include "rotation.hpp"

//baselib dll 
#include "baselib.h"


#include "defines.hpp"


#ifndef PI
#define PI 3.1415926
#endif

/*
// go from a vector representing a rotation in 
// axis-angle format to a 3x3 rotation matrix 
// in column major form
void aa2rot(double const * x, double * R)
{
	const double epsilon = 1e-18;
	double theta = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	// dblas::norm2(3,const_cast<double *>(x));
	double wx = x[0]/(theta+epsilon);
	double wy = x[1]/(theta+epsilon);
	double wz = x[2]/(theta+epsilon);

	double costheta = cos(theta);
	double sintheta = sin(theta);

	R[0] = costheta + wx*wx*(1-costheta);
	R[1] = wz*sintheta + wx*wy*(1-costheta);
	R[2] = -wy*sintheta + wx*wz*(1-costheta);
	R[3] = wx*wy*(1-costheta) - wz*sintheta;
	R[4] = costheta + wy*wy*(1-costheta);
	R[5] = wx*sintheta + wy*wz*(1-costheta);
	R[6] = wy*sintheta + wx*wz*(1-costheta);
	R[7] = -wx*sintheta + wy*wz*(1-costheta);
	R[8] = costheta + wz*wz*(1-costheta);
};
*/


// rotation to axis angle vector format
void rot2aa(double *R, double *aa)
{
	double tr = R[0] + R[4] + R[8];
	double costheta = CLAMP(0.5 * (tr - 1.0), -1, 1);

	double RT[9], RRT[9];
	dll_matrix_transpose(3, 3, R, RT);
	dll_matrix_diff(3, 3, 3, 3, R, RT, RRT);
	double sintheta = dll_matrix_norm(3, 3, RRT)/sqrt(8.0);

	double theta = atan2(sintheta,costheta);
	double factor = theta / (2.0 * sintheta + 1.0e-10);

	aa[0] = factor * (R[7]-R[5]);
	aa[1] = factor * (R[2]-R[6]);
	aa[2] = factor * (R[3]-R[1]);
}


int  eular2rot(double* R, double pitch, double roll, double yaw)
{
	
	/*
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
	*/


	roll  = roll  / 180.0 * PI;     //y
	pitch = pitch / 180.0 * PI;     //x
	yaw   = yaw   / 180.0 * PI;     //z

	double roty[9];
	double rotx[9];
	double rotz[9];

	//y
	roty[0] = cos(roll);      roty[1] = 0;    roty[2] = sin(roll);
	roty[3] = 0;              roty[4] = 1;    roty[5] = 0;
	roty[6] = -sin(roll);     roty[7] = 0;    roty[8] = cos(roll);

	//x
	rotx[0] = 1;     rotx[1] = 0;             rotx[2] = 0;
	rotx[3] = 0;     rotx[4] = cos(pitch);    rotx[5] = -sin(pitch);
	rotx[6] = 0;     rotx[7] = sin(pitch);    rotx[8] = cos(pitch);

	//z
	rotz[0] = cos(yaw);     rotz[1] = -sin(yaw);    rotz[2] = 0;
	rotz[3] = sin(yaw);     rotz[4] =  cos(yaw);    rotz[5] = 0;
	rotz[6] = 0;            rotz[7] =  0;           rotz[8] = 1;
	
	double R1[9];
	mult(rotx, rotz, R1, 3, 3, 3);
	mult(roty, R1, R, 3, 3, 3);

	//double inRoll  = roll;
	//double inPitch = pitch;
	//double inYaw   = yaw;

	//R[0] = cos(inRoll)*cos(inYaw) + sin(inPitch)*sin(inRoll)*sin(inYaw);
	//R[1] = cos(inPitch)*sin(inYaw);
	//R[2] = -sin(inRoll)*cos(inYaw) + sin(inPitch)*cos(inRoll)*sin(inYaw);
	//R[3] = -cos(inRoll)*sin(inYaw) + sin(inPitch)*sin(inRoll)*cos(inYaw);
	//R[4] = cos(inPitch)*cos(inYaw);
	//R[5] = sin(inRoll)*sin(inYaw) + sin(inPitch)*cos(inRoll)*cos(inYaw);
	//R[6] = cos(inPitch)*sin(inRoll);
	//R[7] = -sin(inPitch); // + or - ?
	//R[8] = cos(inRoll)*cos(inPitch);


	return 0;
}
