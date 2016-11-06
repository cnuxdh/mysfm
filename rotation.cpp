
#include "math.h"

#include "rotation.hpp"

//baselib dll 
#include "baselib.h"

//imagelib
#include "defines.h"


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


