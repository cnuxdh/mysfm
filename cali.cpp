

#include "cali.hpp"
#include "baselib.h"

//#include "matrix/matrix.h"


/*
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#ifdef OPENCV_1X 
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#else
//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
using namespace cv;
#endif
*/



void InitializeCameraParams(CameraPara camPara, camera_params_t &camera)
{
	dll_matrix_ident(3, camera.R);
	camera.t[0] = camera.t[1] = camera.t[2] = 0.0;
	camera.f = 0.0;
	camera.k[0] = camera.k[1] = 0.0;

	camera.k_inv[0] = camera.k_inv[2] = camera.k_inv[3] = 0.0;
	camera.k_inv[4] = camera.k_inv[5] = 0.0;
	camera.k_inv[1] = 1.0;

	camera.f_scale = 1.0;
	camera.k_scale = 1.0;

	for (int i = 0; i < NUM_CAMERA_PARAMS; i++) 
	{
		camera.constrained[i] = 0;
		camera.constraints[i] = 0.0;
		camera.weights[i] = 0.0;
	}
	camera.known_intrinsics = 0;
    
	camera.f    = camPara.focalLen;
	camera.k[0] = camPara.k1;
	camera.k[1] = camPara.k2;
	memcpy(camera.R, camPara.R, sizeof(double)*9);
	memcpy(camera.t, camPara.T, sizeof(double)*3);
}


void InitializeCameraParams(camera_params_t &camera)
{
	dll_matrix_ident(3, camera.R);
	camera.t[0] = camera.t[1] = camera.t[2] = 0.0;
	camera.f = 0.0;
	camera.k[0] = camera.k[1] = 0.0;

	camera.k_inv[0] = camera.k_inv[2] = camera.k_inv[3] = 0.0;
	camera.k_inv[4] = camera.k_inv[5] = 0.0;
	camera.k_inv[1] = 1.0;

	camera.f_scale = 1.0;
	camera.k_scale = 1.0;

	for (int i = 0; i < NUM_CAMERA_PARAMS; i++) 
	{
		camera.constrained[i] = 0;
		camera.constraints[i] = 0.0;
		camera.weights[i] = 0.0;
	}

	/*
	camera.fisheye = data.m_fisheye;
	camera.f_cx = data.m_fCx;
	camera.f_cy = data.m_fCy;
	camera.f_rad = data.m_fRad;
	camera.f_angle = data.m_fAngle;
	camera.f_focal = data.m_fFocal;
    */

	camera.known_intrinsics = 0;

	/*
	if (data.m_known_intrinsics) {
		camera.known_intrinsics = 1;
		memcpy(camera.K_known, data.m_K, 9 * sizeof(double));
		memcpy(camera.k_known, data.m_k, 5 * sizeof(double));
	} else {
		camera.known_intrinsics = 0;
	}
	*/
}


DLL_EXPORT void CreateChessBoard()
{

	/*
	int blockSize=75;
	int imageSize=blockSize*8;
	Mat chessBoard(imageSize,imageSize,CV_8UC3,Scalar::all(0));
	unsigned char color=0;

	for(int i=0;i<imageSize;i=i+blockSize)
	{
		color=~color;
		for(int j=0;j<imageSize;j=j+blockSize)
		{
			Mat ROI=chessBoard(Rect(i,j,blockSize,blockSize));
			ROI.setTo(Scalar::all(color));
			color=~color;
		}
	}

	imshow("Chess board", chessBoard);
	*/
}



