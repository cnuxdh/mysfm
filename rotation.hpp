
#ifndef ROTATION_HPP
#define ROTATION_HPP



//void aa2rot(double const * x, double * R);
void rot2aa(double *R, double *aa);


template<typename T>
void aa2rot(const T* x, T * R)
{
	const T epsilon = T(1e-18);
	T theta = sqrt(x[0]*x[0]+x[1]*x[1]+x[2]*x[2]);
	// dblas::norm2(3,const_cast<double *>(x));
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



//calculate the Eular angle from rotation matrix
int  rot2eular(double* R, double* ea);




#endif

