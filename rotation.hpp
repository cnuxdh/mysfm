
#ifndef ROTATION_HPP
#define ROTATION_HPP


// rotation to axis angle vector format
void rot2aa(double *R, double *aa);


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



//calculate the Eular angle from rotation matrix
template<typename T>
int  rot2eular(T* R, T* ea)
{
	ea[0] = atan2( R[5], R[8] ) / PI * T(180.0);  //x:pitch
	ea[1] = asin( -R[2] ) / PI * T(180.0);        //y:roll
	ea[2] = atan2( R[1], R[0]) / PI * T(180.0);   //z:yaw

	return 0;
}




#endif

