
#include "assert.h"
#include "math.h"



#ifdef OPENCV_1X 
	#include "cv.h"
	#include "highgui.h"
	#include "cxcore.h"
#else
	#include "opencv2/core/core.hpp"
	#include "opencv2/highgui/highgui.hpp"
	#include "opencv2/calib3d/calib3d.hpp"
	using namespace cv;
#endif




#include "sift.hpp"
#include "relativepose.hpp"
#include "cali.hpp"
#include "distortion.hpp"
//#include "ba.hpp"
#include "panorama.hpp"
#include "rotation.hpp"
#include "triangulate.hpp"


//corelib
#include "Matrix.h"
#include "commonfile.h"

#include "CalcAngle.h"

//sfm driver lib
//#include "sfm.h"

//5point lib 
//#include "5point.h"

//matrix lib
//#include "matrix/matrix.h"


//imagelib
//#include "defines.h"
//#include "triangulate.h"
//#include "fmatrix.h"

#include "baselib.h"


#include <algorithm>
using namespace std;


//convert from RX+T to R(X - (-inv(R)*T) )
int CalculateExplicitT(double* R, double* T, double* explicitT)
{
	double iR[9];
	memset(iR, 0, sizeof(double)*9);
	dll_matrix_invert(3, R, iR);

	double eT[3];
	dll_matrix_product(3,3,3,1, iR, T, eT);
	
	explicitT[0] = -eT[0];
	explicitT[1] = -eT[1];
	explicitT[2] = -eT[2];
	
	return 0;
}


/* Use a 180 rotation to fix up the intrinsic matrix */
void FixIntrinsics(double *P, double *K, double *R, double *t) 
{
	/* Check the parity along the diagonal */
	int neg = (K[0] < 0.0) + (K[4] < 0.0) + (K[8] < 0.0);

	/* If odd parity, negate the instrinsic matrix */
	if ((neg % 2) == 1) 
	{
		dll_matrix_scale(3, 3, K, -1.0, K);
		dll_matrix_scale(3, 4, P, -1.0, P);
	}

	/* Now deal with case of even parity */
	double fix[9];
	dll_matrix_ident(3, fix);
	double tmp[9], tmp2[12];

	if (K[0] < 0.0 && K[4] < 0.0) 
	{
		fix[0] = -1.0;
		fix[4] = -1.0;
	} 
	else if (K[0] < 0.0) 
	{
		fix[0] = -1.0;
		fix[8] = -1.0;
	} 
	else if (K[4] < 0.0) 
	{
		fix[4] = -1.0;
		fix[8] = -1.0;
	} else {
		/* No change needed */
	}

	dll_matrix_product(3, 3, 3, 3, K, fix, tmp);
	memcpy(K, tmp, sizeof(double) * 3 * 3);

	double Kinv[9];
	dll_matrix_invert(3, K, Kinv);

	dll_matrix_product(3, 3, 3, 4, Kinv, P, tmp2);

	memcpy(R + 0, tmp2 + 0, sizeof(double) * 3);
	memcpy(R + 3, tmp2 + 4, sizeof(double) * 3);
	memcpy(R + 6, tmp2 + 8, sizeof(double) * 3);

	t[0] = tmp2[3];
	t[1] = tmp2[7];
	t[2] = tmp2[11];
}


/* Estimate relative pose from a given set of point matches for panorama
   Model: x = RX+T
   threshold: 0.25*9 default
*/
int EstimatePose5Point_Pano( vector<Point3DDouble>& p1, 
							 vector<Point3DDouble>& p2,
							 double radius,
						     int num_trials, double threshold, 
						     double *R, double *t, vector<double>& residual)
{
	int num_pts = (int) p1.size();

	v3_t *k1_pts = new v3_t[num_pts];
	v3_t *k2_pts = new v3_t[num_pts];

	for (int i = 0; i < num_pts; i++) 
	{
		k1_pts[i] = dll_v3_new(p1[i].p[0], p1[i].p[1], p1[i].p[2]);
		k2_pts[i] = dll_v3_new(p2[i].p[0], p2[i].p[1], p2[i].p[2]);
	}

	double em[9];
	int num_inliers = dll_compute_pose_ransac_pano(num_pts, k1_pts, k2_pts, 
		radius, threshold, num_trials, em, R, t);
    
	//printf("residuals ... \n");
	for (int i = 0; i < num_pts; i++) 
	{
		double dis = dll_fmatrix_compute_residual_pano(em, k1_pts[i], k2_pts[i], radius);
		//printf("%d  %lf \n", i, dis);
		residual[i] = dis;
	}
	
	delete [] k1_pts;
	delete [] k2_pts;

	return num_inliers;
}

/* Estimate relative pose from a given set of point matches 
   Model: x = RX+T
   threshold: 0.25*9 default
*/
int EstimatePose5Point( vector<Point2DDouble>& p1, 
						vector<Point2DDouble>& p2,
						int num_trials, double threshold, 
						double *K1, double *K2, 
						double *R, double *t)

{
	int num_pts = (int) p1.size();

	v2_t *k1_pts = new v2_t[num_pts];
	v2_t *k2_pts = new v2_t[num_pts];

	for (int i = 0; i < num_pts; i++) 
	{
		k1_pts[i] = dll_v2_new(p1[i].p[0], p1[i].p[1]);
		k2_pts[i] = dll_v2_new(p2[i].p[0], p2[i].p[1]);
	}

	int num_inliers = dll_compute_pose_ransac(num_pts, k1_pts, k2_pts, 
		K1, K2, threshold, num_trials, R, t);

	delete [] k1_pts;
	delete [] k2_pts;

	return num_inliers;
}



/* Check cheirality for a camera and a point */
bool CheckCheirality(v3_t p, const camera_params_t &camera) 
{
	double pt[3] = { Vx(p), Vy(p), Vz(p) };
	double cam[3];

	pt[0] -= camera.t[0];
	pt[1] -= camera.t[1];
	pt[2] -= camera.t[2];
	dll_matrix_product(3, 3, 3, 1, (double *) camera.R, pt, cam);
		
	// EDIT!!!
	if (cam[2] > 0.0)
		return false;
	else
		return true;
}


/*
v2_t UndistortNormalizedPoint(v2_t p, camera_params_t c) 
{
	double r = sqrt(Vx(p) * Vx(p) + Vy(p) * Vy(p));
	if (r == 0.0)
		return p;

	double t = 1.0;
	double a = 0.0;

	for (int i = 0; i < POLY_INVERSE_DEGREE; i++) 
	{
		a += t * c.k_inv[i];
		t = t * r;
	}

	double factor = a / r;

	return v2_scale(factor, p);
}
*/


/* Compute the angle between two rays */
double ComputeRayAngle(v2_t p, v2_t q, 
					   camera_params_t &cam1, 
					   camera_params_t &cam2)
{
	double K1[9], K2[9];
	GetIntrinsics(cam1, K1);
	GetIntrinsics(cam2, K2);

	double K1_inv[9], K2_inv[9];
	dll_matrix_invert(3, K1, K1_inv);
	dll_matrix_invert(3, K2, K2_inv);

	double p3[3] = { Vx(p), Vy(p), 1.0 };
	double q3[3] = { Vx(q), Vy(q), 1.0 };

	double p3_norm[3], q3_norm[3];
	dll_matrix_product331(K1_inv, p3, p3_norm);
	dll_matrix_product331(K2_inv, q3, q3_norm);

	v2_t p_norm = dll_v2_new(p3_norm[0] / p3_norm[2], p3_norm[1] / p3_norm[2]);
	v2_t q_norm = dll_v2_new(q3_norm[0] / q3_norm[2], q3_norm[1] / q3_norm[2]);

	double R1_inv[9], R2_inv[9];
	dll_matrix_transpose(3, 3, (double *) cam1.R, R1_inv);
	dll_matrix_transpose(3, 3, (double *) cam2.R, R2_inv);

	double p_w[3], q_w[3];

	double pv[3] = { Vx(p_norm), Vy(p_norm), -1.0 };
	double qv[3] = { Vx(q_norm), Vy(q_norm), -1.0 };

	double Rpv[3], Rqv[3];

	dll_matrix_product331(R1_inv, pv, Rpv);
	dll_matrix_product331(R2_inv, qv, Rqv);

	dll_matrix_sum(3, 1, 3, 1, Rpv, (double *) cam1.t, p_w);
	dll_matrix_sum(3, 1, 3, 1, Rqv, (double *) cam2.t, q_w);

	/* Subtract out the camera center */
	double p_vec[3], q_vec[3];
	dll_matrix_diff(3, 1, 3, 1, p_w, (double *) cam1.t, p_vec);
	dll_matrix_diff(3, 1, 3, 1, q_w, (double *) cam2.t, q_vec);

	/* Compute the angle between the rays */
	double dot;
	dll_matrix_product(1, 3, 3, 1, p_vec, q_vec, &dot);

	double mag = dll_matrix_norm(3, 1, p_vec) * dll_matrix_norm(3, 1, q_vec);

	return acos(CLAMP(dot / mag, -1.0 + 1.0e-8, 1.0 - 1.0e-8));
}



//////////////////////////////////////////////////////////////////////////
CEstimatePose5Point::CEstimatePose5Point()
{
	m_bIsExplicit = false;
}
CEstimatePose5Point::~CEstimatePose5Point()
{

}
/* 
inputs:
  pairMatches:         match pair
  lImgFeat,rImageFeat: feature point 
outputs:
  cam1,cam2:           camera structure
*/
int CEstimatePose5Point::EstimatePose( PairMatchRes pairMatches, ImgFeature& lImageFeat, ImgFeature& rImageFeat, 
									   CameraPara& cam1, CameraPara& cam2 )
{
	//sfm interface
	camera_params_t camera1;
	camera_params_t camera2;
	camera1.f = cam1.focalLen;
	camera2.f = cam2.focalLen;
	camera1.known_intrinsics = 0;
	camera2.known_intrinsics = 0;

	/* Put first camera at origin */
	camera1.R[0] = 1.0;  camera1.R[1] = 0.0;  camera1.R[2] = 0.0;
	camera1.R[3] = 0.0;  camera1.R[4] = 1.0;  camera1.R[5] = 0.0;
	camera1.R[6] = 0.0;  camera1.R[7] = 0.0;  camera1.R[8] = 1.0;
	camera1.t[0] = 0.0;  camera1.t[1] = 0.0;  camera1.t[2] = 0.0;


	double K1[9], K2[9];
	GetIntrinsics(camera1, K1);
	GetIntrinsics(camera2, K2);
    
	vector<Point2DDouble> p1;
	vector<Point2DDouble> p2;

	int numPt = pairMatches.matchs.size();
	for(int i=0; i<numPt; i++)
	{
		Point2DDouble tp1,tp2;
		int li = pairMatches.matchs[i].l;
		int ri = pairMatches.matchs[i].r;
		
		//normalize the coordinates
		tp1.p[0] = lImageFeat.featPts[li].cx;
		tp1.p[1] = lImageFeat.featPts[li].cy;
		tp2.p[0] = lImageFeat.featPts[li].cx;
		tp2.p[1] = lImageFeat.featPts[li].cy;

		p1.push_back(tp1);
		p2.push_back(tp2);
	}

	double R0[9], t0[3];
	int num_inliers = 0;
    double fmatrix_threshold = 9;
	num_inliers = EstimatePose5Point(p1, p2, 
		512,/* m_fmatrix_rounds, 8 * m_fmatrix_rounds */ 
		0.25*fmatrix_threshold, // 0.003, // 0.004 /*0.001,*/ // /*0.5 **/ fmatrix_threshold, 
		K1, K2, R0, t0 );
	
	printf("number of inliers : %d \n", num_inliers);
    

	//bool initialized = false;
	//if (!initialized) 
	{
		memcpy(camera2.R, R0, sizeof(double) * 9);
		dll_matrix_transpose_product(3, 3, 3, 1, R0, t0, camera2.t);
		dll_matrix_scale(3, 1, camera2.t, -1.0, camera2.t);
	}
	m_bIsExplicit = true;

	memcpy(cam1.R, camera1.R, sizeof(double)*9);
	memcpy(cam1.T, camera1.t, sizeof(double)*3);
	memcpy(cam2.R, camera2.R, sizeof(double)*9);
	memcpy(cam2.T, camera2.t, sizeof(double)*3);
   
	return 1;
}

/*
  Output model: x = R(X-T)
  inputs:
       lpts,rpts: normalized point coordinates of image, the origin is the center of image
  output:
	   cam2: the results are output into cam2
*/
int CEstimatePose5Point::EstimatePose( vector<Point2DDouble>& lPts, vector<Point2DDouble>& rPts, 
									  CameraPara& cam1, CameraPara& cam2 )
{
	//printf(" \n Relative Pose Estimation ... \n");

	double K1[9],K2[9],R[9],T[3],R0[9],t0[3];
	memset(K1, 0, sizeof(double)*9);
	memset(K2, 0, sizeof(double)*9);
	memset(R0, 0, sizeof(double)*9);
	memset(t0, 0, sizeof(double)*3);
	memset(R, 0, sizeof(double)*9);
	memset(T, 0, sizeof(double)*3);

	double focus1 = cam1.focalLen; //947.2;
	double focus2 = cam2.focalLen; //947.2;
	K1[0] = focus1;  K1[4] = focus1; K1[8] = 1;
	K2[0] = focus2;  K2[4] = focus2; K2[8] = 1;
	R0[0] = R0[4] = R0[8] = 1;    

	memcpy(cam1.R, R0, sizeof(double)*9);
	memcpy(cam1.T, t0, sizeof(double)*3);

	int num_inliers =  EstimatePose5Point(lPts, rPts, 512, 4,  K1, K2, R, T);
	printf("number of inliers : %d \n", num_inliers);

	printf("RX+T - translation: %lf %lf %lf \n", T[0], T[1], T[2]);

	//decompose the three angles around x,y,z axis
	cam2.ax = atan( R[5]/R[8] )/PI*180; 
	cam2.ay = asin( -R[2] )/PI*180;
	cam2.az = atan( R[1]/R[0])/PI*180;

	double nt[3];
	//bool initialized = false;
	//if (!initialized) 
	{
		memcpy(cam2.R, R, sizeof(double) * 9);
		
		dll_matrix_transpose_product(3, 3, 3, 1, R, T, nt);
		dll_matrix_scale(3, 1, nt, -1.0, nt);
		memcpy(cam2.T,nt,sizeof(double)*3);
	}	
   
	m_bIsExplicit = true;

	/*
	//print results
#ifdef _DEBUG
	
	FILE* fp = fopen("c:\\temp\\myRelativePose.txt", "w");
	fprintf(fp, "%d %d \n", 0, 1);
	fprintf(fp, "focalLen: %lf \n", cam2.focalLen);
	fprintf(fp, "Rotation Matrix: \n");
	for(int j=0; j<3; j++)
	{
		for(int i=0; i<3; i++)
		{
			fprintf(fp,"%lf ", cam2.R[j*3+i]);
		}
		fprintf(fp,"\n");
	}
	fprintf(fp, "Translation: \n");
	for(int j=0; j<3; j++)
		fprintf(fp,"%lf ", cam2.T[j]);
	fprintf(fp, "\n");
	fclose(fp);
#endif
	*/

	////output
	//printf("focalLen: %lf \n", cam2.focalLen);
	//printf("Rotation Angle x:%lf y:%lf z:%lf \n", cam2.ax, cam2.ay, cam2.az);
	//for(int j=0; j<3; j++)
	//{
	//	for(int i=0; i<3; i++)
	//	{
	//		printf("%lf ", cam2.R[j*3+i]);
	//	}
	//	printf("\n");
	//}
	//printf("Translation: \n");
	//for(int j=0; j<3; j++)
	//	printf("%lf ", cam2.T[j]);
	//printf("\n\n");

  	return 1;
}


CEstimatePose5PointPano::CEstimatePose5PointPano()
{
	m_bIsExplicit = false;
}

CEstimatePose5PointPano::~CEstimatePose5PointPano()
{

}

/*
   Output model: x = R(X-T)   
*/
int CEstimatePose5PointPano::EstimatePose( vector<Point2DDouble>& lPts, vector<Point2DDouble>& rPts, 
	CameraPara& cam1, CameraPara& cam2 )
{
	int ht = cam1.rows;
	int wd = cam1.cols;
	double radius = (double)(wd) / (2*PI);

	vector<Point3DDouble> pl;
	vector<Point3DDouble> pr;

	//from pano 2d to 3d 
	for(int i=0; i<lPts.size(); i++)
	{
		double x = lPts[i].p[0];
		double y = lPts[i].p[1];
		//from image center to top-left
		//x = wd*0.5 + x;
		//y = ht*0.5 - y;
		double gx,gy,gz;
		SphereTo3D_center(x, y, radius, gx, gy, gz);
		Point3DDouble p3;
		p3.p[0] = gx / radius;
		p3.p[1] = gy / radius;
		p3.p[2] = gz / radius;
		pl.push_back(p3);


		x = rPts[i].p[0];
		y = rPts[i].p[1];
		//from image center to top-left
		//x = wd*0.5 + x;
		//y = ht*0.5 - y;
		SphereTo3D_center(x, y, radius, gx, gy, gz);
		p3.p[0] = gx / radius;
		p3.p[1] = gy / radius;
		p3.p[2] = gz / radius;
		pr.push_back(p3);
	}
	
	int num_trials = 512;
	double threshold = 2.5;
	double R[9];
	double T[3];

	vector<double> residual;
	residual.resize(pl.size());
	EstimatePose5Point_Pano(pl, pr, radius, num_trials, threshold, R, T, residual);

	/*
	//remove the outliers
	double epipolarThreshold = wd*0.005;
	vector<Point2DDouble> inlierLeftPts, inlierRightPts;
	for(int i=0; i<lPts.size(); i++)
	{
		if(residual[i]<epipolarThreshold)
		{
			inlierLeftPts.push_back(lPts[i]);
			inlierRightPts.push_back(rPts[i]);
		}
	}
	lPts = inlierLeftPts;
	rPts = inlierRightPts;
	*/

	//rotation matrix
	for(int i=0; i<9; i++)
	{
		cam1.R[i] = 0;
		cam2.R[i] = R[i];
	}
	cam1.R[0] = cam1.R[4] = cam1.R[8] = 1;


	//from RX+T to R(X-T')
	double Rt[9];
	transpose(R, Rt, 3, 3);
	double gt[3];
	mult(Rt, T, gt, 3, 3, 1);
	for(int i=0; i<3; i++)
	{
		cam1.T[i] = 0;
		cam2.T[i] = -gt[i];
	}

	m_bIsExplicit = true;
	
	//save eular angle
	double ea[3];
	rot2eular(R, ea);
	cam1.ax = cam1.ay = cam1.az = 0;
	cam1.bIsExplicit = m_bIsExplicit;

	cam2.ax = ea[0];
	cam2.ay = ea[1];
	cam2.az = ea[2];
	cam2.bIsExplicit = m_bIsExplicit;

	printf("relative pose eular angle: %lf %lf %lf \n", ea[0], ea[1], ea[2]);
	printf("relative pose translation: %lf %lf %lf \n", cam2.T[0], cam2.T[1], cam2.T[2]);
	
	return 0;
}

int CEstimatePose5PointPano::EstimatePose( vector<Point2DDouble>& lPts, vector<Point2DDouble>& rPts, 
	CameraPara& cam1, CameraPara& cam2, vector<int>& inliers )
{
	int ht = cam1.rows;
	int wd = cam1.cols;
	double radius = (double)(wd) / (2*PI);

	vector<Point3DDouble> pl;
	vector<Point3DDouble> pr;

	//from pano 2d to 3d 
	for(int i=0; i<lPts.size(); i++)
	{
		double x = lPts[i].p[0];
		double y = lPts[i].p[1];
		//from image center to top-left
		//x = wd*0.5 + x;
		//y = ht*0.5 - y;
		double gx,gy,gz;
		SphereTo3D_center(x, y, radius, gx, gy, gz);
		Point3DDouble p3;
		p3.p[0] = gx / radius;
		p3.p[1] = gy / radius;
		p3.p[2] = gz / radius;
		pl.push_back(p3);


		x = rPts[i].p[0];
		y = rPts[i].p[1];
		//from image center to top-left
		//x = wd*0.5 + x;
		//y = ht*0.5 - y;
		SphereTo3D_center(x, y, radius, gx, gy, gz);
		p3.p[0] = gx / radius;
		p3.p[1] = gy / radius;
		p3.p[2] = gz / radius;
		pr.push_back(p3);
	}
	
	int num_trials = 512;
	double threshold = 2.5;
	double R[9];
	double T[3];

	vector<double> residual;
	residual.resize(pl.size());
	EstimatePose5Point_Pano(pl, pr, radius, num_trials, threshold, R, T, residual);

	//save the inliers
	double epipolarThreshold = wd*0.005;
	inliers.clear();
	for(int i=0; i<lPts.size(); i++)
	{
		if(residual[i]<epipolarThreshold)
		{
			inliers.push_back(i);
		}
	}

	//rotation matrix
	for(int i=0; i<9; i++)
	{
		cam1.R[i] = 0;
		cam2.R[i] = R[i];
	}
	cam1.R[0] = cam1.R[4] = cam1.R[8] = 1;


	//from RX+T to R(X-T')
	double Rt[9];
	transpose(R, Rt, 3, 3);
	double gt[3];
	mult(Rt, T, gt, 3, 3, 1);
	for(int i=0; i<3; i++)
	{
		cam1.T[i] = 0;
		cam2.T[i] = -gt[i];
	}

	m_bIsExplicit = true;
	
	//save eular angle
	double ea[3];
	rot2eular(R, ea);
	cam1.ax = cam1.ay = cam1.az = 0;
	cam1.bIsExplicit = m_bIsExplicit;

	cam2.ax = ea[0];
	cam2.ay = ea[1];
	cam2.az = ea[2];
	cam2.bIsExplicit = m_bIsExplicit;

	printf("relative pose eular angle: %lf %lf %lf \n", ea[0], ea[1], ea[2]);
	printf("relative pose translation: %lf %lf %lf \n", cam2.T[0], cam2.T[1], cam2.T[2]);
	

	return 0;
}


//////////////////////////////////////////////////////////////////////////
CDLTPose::CDLTPose()
{

}

CDLTPose::~CDLTPose()
{

}

// the projection matrix: P=K[R|T]
int CDLTPose::EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, double* K, double* R, double* t)
{	
	double P[12];
	int r = -1;
	int num_points = pt3.size();

	//too few points
	if (num_points < 8)
	{
		printf("too few dlt point number: %d \n", num_points);
		return -1;
	}

	v3_t *points_solve = (v3_t*)malloc( num_points*sizeof(v3_t) );
	v2_t *projs_solve  = (v2_t*)malloc( num_points*sizeof(v2_t) );

	for(int i=0; i<num_points; i++)
	{
		points_solve[i].p[0] = pt3[i].p[0];
		points_solve[i].p[1] = pt3[i].p[1];
		points_solve[i].p[2] = pt3[i].p[2];

		projs_solve[i].p[0] = pt2[i].p[0];
		projs_solve[i].p[1] = pt2[i].p[1];
	}

	double proj_estimation_threshold = 4.0;
	double proj_estimation_threshold_weak = 16.0*proj_estimation_threshold;

#ifdef _DEBUG
	FILE* fp = fopen("c:\\temp\\my_dlt.txt", "w");
	for(int i=0; i<num_points; i++)
	{
		fprintf(fp, "%lf %lf %lf   %lf %lf \n", 
			points_solve[i].p[0], points_solve[i].p[1], points_solve[i].p[2],
			projs_solve[i].p[0], projs_solve[i].p[1]);
	}
	fclose(fp);
#endif

	//find the 3*4 projection matrix
	r = dll_find_projection_3x4_ransac(num_points, 
		points_solve, projs_solve, 
		P, /* 2048 */ 4096 /* 100000 */, 
		proj_estimation_threshold);

	if (r == -1 ) 
	{
		printf("[FindAndVerifyCamera] Couldn't find projection matrix\n");
		free(points_solve);
		free(projs_solve);

		return -1;
	}
	
	/* If number of inliers is too low, fail */
	if ( r <= MIN_INLIERS_EST_PROJECTION) 
	{
		printf("[FindAndVerifyCamera] Too few inliers to use projection matrix\n");
		free(points_solve);
		free(projs_solve);

		return -1;
	}

	//decompose the P to K[R|t]
	double KRinit[9], Kinit[9], Rinit[9], tinit[3];
	memcpy(KRinit + 0, P + 0, 3 * sizeof(double));
	memcpy(KRinit + 3, P + 4, 3 * sizeof(double));
	memcpy(KRinit + 6, P + 8, 3 * sizeof(double));

	// RQ factorization of KR
	dll_dgerqf_driver(3, 3, KRinit, Kinit, Rinit);	    

	/* We want our intrinsics to have a certain form */
	FixIntrinsics(P, Kinit, Rinit, tinit);
	dll_matrix_scale(3, 3, Kinit, 1.0 / Kinit[8], Kinit);

	printf("[FindAndVerifyCamera] Estimated intrinsics:\n");
	dll_matrix_print(3, 3, Kinit);
	printf("[FindAndVerifyCamera] Estimated extrinsics:\n");
	dll_matrix_print(3, 3, Rinit);
	dll_matrix_print(1, 3, tinit);
	fflush(stdout);

	/* Check cheirality constraint */
	printf("[FindAndVerifyCamera] Checking consistency...\n");

	double Rigid[12] = 
	{ Rinit[0], Rinit[1], Rinit[2], tinit[0],
	  Rinit[3], Rinit[4], Rinit[5], tinit[1],
	  Rinit[6], Rinit[7], Rinit[8], tinit[2] };

	//std::vector<int> inliers;
	//std::vector<int> inliers_weak;
	std::vector<int> outliers;

	//int *idxs_solve = (int*)malloc( sizeof(int)*num_points );

	inliers.clear();
	inliers_weak.clear();
	outliers.clear();

	int num_behind = 0;
	for (int j = 0; j < num_points; j++) 
	{
		double p[4] = { Vx(points_solve[j]), 
						Vy(points_solve[j]),
						Vz(points_solve[j]), 1.0 };
		double q[3], q2[3];

		dll_matrix_product(3, 4, 4, 1, Rigid, p, q);
		dll_matrix_product331(Kinit, q, q2);

		double pimg[2] = { -q2[0] / q2[2], -q2[1] / q2[2] };
		double diff = 
			(pimg[0] - Vx(projs_solve[j])) * 
			(pimg[0] - Vx(projs_solve[j])) + 
			(pimg[1] - Vy(projs_solve[j])) * 
			(pimg[1] - Vy(projs_solve[j]));

		diff = sqrt(diff);

		if (diff < proj_estimation_threshold)
			inliers.push_back(j);

		if (diff < proj_estimation_threshold_weak) 
		{
			inliers_weak.push_back(j);
		}
		else 
		{
			printf("[FindAndVerifyCamera] Removing point [%d] "
				"(reproj. error = %0.3f)\n", j, diff);
			outliers.push_back(j);
		}

		// EDIT!!!
		if (q[2] > 0.0)
			num_behind++;  /* Cheirality constraint violated */
	}

	//free(idxs_solve);
	free(points_solve);
	free(projs_solve);

	if (num_behind >= 0.9 * num_points) 
	{
		printf("[FindAndVerifyCamera] Error: camera is pointing "
			"away from scene\n");

		return -1;
	}

	memcpy(K, Kinit, sizeof(double) * 9);
	memcpy(R, Rinit, sizeof(double) * 9);
	memcpy(t, tinit, sizeof(double) * 3);	

	//conversion from RX+t to R(X-t')
	printf(" Conversion from RX+t to R(X-t') \n");
	double nt[3];
	dll_matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, nt);
	dll_matrix_scale(3, 1, nt, -1.0, nt);
	dll_matrix_print(1, 3, nt);
	t[0]=nt[0]; t[1]=nt[1]; t[2]=nt[2];

	/*
	//calculate the rotation angle
	double ax = atan( R[5]/R[8] )/PI*180; 
	double ay = asin( -R[2] )/PI*180;
	double az = atan( R[1]/R[0])/PI*180;
	printf("Rotation Angle about x,y,z: %lf %lf %lf \n", ax, ay, az);
	*/
	
	return 0;
}


//Output model: R(X-T)
int CDLTPose::EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, CameraPara& cam)
{
	printf("\n  DLT Pose Estimation ... \n");

	int res = EstimatePose(pt3, pt2, cam.K, cam.R, cam.T);
	
	//double ratio = 
	cam.focalLen = (cam.K[0]+cam.K[4])*0.5;
	
	//calculate the rotation angle
	cam.ax = atan2( cam.R[5], cam.R[8] )/PI*180; //atan( cam.R[5]/cam.R[8] )/PI*180; 
	cam.ay = asin( -cam.R[2] )/PI*180;
	cam.az = atan2(cam.R[1], cam.R[0]) /PI*180;  //atan( cam.R[1]/cam.R[0])/PI*180;

	printf("Rotation Angle about x,y,z: %lf %lf %lf \n", cam.ax, cam.ay, cam.az);


	return res;
}

vector<int> CDLTPose::GetInliers()
{
	return inliers;
}
vector<int> CDLTPose::GetWeakInliers()
{
	return inliers_weak;
}



CPanoDLTPose::CPanoDLTPose()
{

}
CPanoDLTPose::~CPanoDLTPose()
{

}

/* Output Model: R(X-T)
   pt3: 3D control points,
   pt2: 2D centered image points
*/
int CPanoDLTPose::EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, CameraPara& cam)
{
	int ht = cam.rows;
	int wd = cam.cols;

	double radius = (double)(wd) / (2*PI);

	//from spherical image to 3D
	vector<Point3DDouble> panoPt3;
	for(int i=0; i<pt2.size(); i++)
	{
		Point3DDouble p3;
		
		double x,y;

		//center point
		x = pt2[i].p[0]; //pt2[i].p[0] + wd*0.5;
		y = pt2[i].p[1]; //ht*0.5 - pt2[i].p[1];

		//calculate the 3D ray light 
		double gx,gy,gz;
		SphereTo3D_center(x,y,radius,gx,gy,gz);
		p3.p[0] = gx;
		p3.p[1] = gy;
		p3.p[2] = gz;
		
		panoPt3.push_back(p3);
	}
	
	//panorama pose estimation
	int res = EstimatePose(pt3, panoPt3, cam);
	
	return res;
}


//pt3 and pt2 are 3D vectors
int CPanoDLTPose::EstimatePose(vector<Point3DDouble> pt3, vector<Point3DDouble> pt2, CameraPara& cam)
{
	double P[12];
	int r = -1;
	int num_points = pt3.size();

	//too few points
	if (num_points <= 6)
	{
		printf("too few points, exit! \n");
		return -1;
	}

	v3_t *points_solve = (v3_t*)malloc( num_points*sizeof(v3_t) );
	v3_t *projs_solve  = (v3_t*)malloc( num_points*sizeof(v3_t) );

	for(int i=0; i<num_points; i++)
	{
		points_solve[i].p[0] = pt3[i].p[0];
		points_solve[i].p[1] = pt3[i].p[1];
		points_solve[i].p[2] = pt3[i].p[2];

		projs_solve[i].p[0] = pt2[i].p[0];
		projs_solve[i].p[1] = pt2[i].p[1];
		projs_solve[i].p[2] = pt2[i].p[2];
	}

	double proj_estimation_threshold = 2; //cam.cols*0.002;
	double proj_estimation_threshold_weak = 4.0*proj_estimation_threshold;

	//find the 3*4 projection matrix
	r = dll_find_projection_3x4_ransac_pano(num_points, 
		points_solve, projs_solve, 
		P, /* 2048 */ 4096 /* 100000 */, 
		proj_estimation_threshold);

	if (r == -1 ) 
	{
		printf("[FindAndVerifyCamera] Couldn't find projection matrix\n");
		free(points_solve);
		free(projs_solve);

		return -1;
	}
	
	/* If number of inliers is too low, fail */
	if ( r <= MIN_INLIERS_EST_PROJECTION) 
	{
		printf("[FindAndVerifyCamera] Too few inliers to use projection matrix\n");
		free(points_solve);
		free(projs_solve);

		return -1;
	}

	//print projective matrix
	if(1)
	{
		FILE* fp = fopen("c:\\temp\\p.txt", "w");
		for(int j=0; j<3; j++)
		{
			for(int i=0; i<4; i++)
			{
				fprintf(fp, "%lf ", P[j*4+i]);
			}
			fprintf(fp, "\n");
		}
		fclose(fp);
	}
	
	//decompose the P to K[R|t]
	double KRinit[9], Kinit[9], Rinit[9], tinit[3];
	memcpy(KRinit + 0, P + 0, 3 * sizeof(double));
	memcpy(KRinit + 3, P + 4, 3 * sizeof(double));
	memcpy(KRinit + 6, P + 8, 3 * sizeof(double));

	// RQ factorization of KR
	dll_dgerqf_driver(3, 3, KRinit, Kinit, Rinit);	    
	printf("[FindAndVerifyCamera] K from RQ factoriztion :\n");
	dll_matrix_print(3, 3, Kinit);
	printf("[FindAndVerifyCamera] R from RQ factoriztion: \n");
	dll_matrix_print(3, 3, Rinit);


	/* We want our intrinsics to have a certain form */
	FixIntrinsics(P, Kinit, Rinit, tinit);
	dll_matrix_scale(3, 3, Kinit, 1.0 / Kinit[8], Kinit);

	printf("[FindAndVerifyCamera] Estimated intrinsics:\n");
	dll_matrix_print(3, 3, Kinit);
	printf("[FindAndVerifyCamera] Estimated extrinsics:\n");
	dll_matrix_print(3, 3, Rinit);
	dll_matrix_print(1, 3, tinit);
	fflush(stdout);

	double Rigid[12] = 
	{ Rinit[0], Rinit[1], Rinit[2], tinit[0],
	  Rinit[3], Rinit[4], Rinit[5], tinit[1],
	  Rinit[6], Rinit[7], Rinit[8], tinit[2] };

	std::vector<int> outliers;

	inliers.clear();
	inliers_weak.clear();
	outliers.clear();

	int num_behind = 0;
	for (int j = 0; j < num_points; j++) 
	{
		double p[4] = { Vx(points_solve[j]), 
						Vy(points_solve[j]),
						Vz(points_solve[j]), 
						1.0 };

		double q[3], q2[3];

		dll_matrix_product(3, 4, 4, 1, Rigid, p, q);
		
		//calculate the angle	
		v3_t projPt;
		projPt.p[0] = q[0];
		projPt.p[1] = q[1];
		projPt.p[2] = q[2];
		double diff = dll_cross_angle(projs_solve[j], projPt);
		
		if (diff < proj_estimation_threshold)
			inliers.push_back(j);

		if (diff < proj_estimation_threshold_weak) 
		{
			inliers_weak.push_back(j);
		}
		else 
		{
			printf("[FindAndVerifyCamera] Removing point [%d] "
				"(reproj. error = %0.3f)\n", j, diff);
			outliers.push_back(j);
		}	
	}

	//free(idxs_solve);
	free(points_solve);
	free(projs_solve);

	double K[9];
	double R[9];
	double t[3];
	memcpy(K, Kinit, sizeof(double) * 9);
	memcpy(R, Rinit, sizeof(double) * 9);
	memcpy(t, tinit, sizeof(double) * 3);	

	//conversion from RX+t to R(X-t')
	printf(" Conversion from RX+t to R(X-t') \n");
	double nt[3];
	dll_matrix_transpose_product(3, 3, 3, 1, Rinit, tinit, nt);
	dll_matrix_scale(3, 1, nt, -1.0, nt);
	dll_matrix_print(1, 3, nt);
	t[0]=nt[0]; t[1]=nt[1]; t[2]=nt[2];

	memcpy(cam.K, K, sizeof(double)*9);
	memcpy(cam.R, R, sizeof(double)*9);
	memcpy(cam.T, t, sizeof(double)*3);	 

	return 0;
}

//////////////////////////////////////////////////////////////////////////
COpenCVFM::COpenCVFM()
{
	//m_FM.create(3, 3, CV_32F);
}

COpenCVFM::~COpenCVFM()
{

}

int COpenCVFM::Estimate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, vector<double>& fm)
{
	assert(lPts.size()==rPts.size());

	fm.resize(9);

	int point_count = lPts.size();
	CvMat* points1;
	CvMat* points2;
	CvMat* status;
	CvMat* fundamental_matrix;
	
	points1 = cvCreateMat(1,point_count,CV_32FC2);
	points2 = cvCreateMat(1,point_count,CV_32FC2);
	status  = cvCreateMat(1,point_count,CV_8UC1);

	/* Fill the points here ... */
	for( int i = 0; i < point_count; i++ )
	{
		points1->data.fl[i*2]   = lPts[i].p[0]; 
		points1->data.fl[i*2+1] = lPts[i].p[1]; 
		points2->data.fl[i*2]   = rPts[i].p[0]; 
		points2->data.fl[i*2+1] = rPts[i].p[1];
	}

	fundamental_matrix = cvCreateMat(3,3,CV_32FC1);
	int fm_count = cvFindFundamentalMat( points1, points2, fundamental_matrix,
										 CV_FM_RANSAC,1.0,0.99,status );
	
	//print result:
	for(int i=0; i<9; i++)
	{
		fm[i]   = fundamental_matrix->data.fl[i];
		m_FM[i] = fundamental_matrix->data.fl[i];
		printf(" %lf ", fm[i]);		
	}

	cvReleaseMat(&points1);
	cvReleaseMat(&points2);
	cvReleaseMat(&status);
	cvReleaseMat(&fundamental_matrix);
	return 1;
}

// x2.FM.X1 = 0
void COpenCVFM::CalculateEpipolarLine(Point2DDouble pt, int flag, double& a, double& b, double& c)
{
	double ip[3];
	ip[0] = pt.p[0];
	ip[1] = pt.p[1];
	ip[2] = 1;

	double line[3];

	//point from left image
	if(flag==0)
	{
		dll_matrix_product(3,3,3,1,m_FM, ip, line);
		//mult(m_FM, ip, line, 3, 3, 1);
	}

	////point from right image
	if(flag==1)
	{
		dll_matrix_product(1,3,3,3, ip, m_FM, line);
		//mult(ip, m_FM, line, 1, 3, 3);
	}

	a = line[0];
	b = line[1];
	c = line[2];
}

double COpenCVFM::CalculateError(Point2DDouble lpt, Point2DDouble rpt)
{
	double dis = 0;

	double lp[3];
	lp[0] = lpt.p[0];
	lp[1] = lpt.p[1];
	lp[2] = 1;

	double rp[3];
	rp[0] = rpt.p[0];
	rp[1] = rpt.p[1];
	rp[2] = 1;

	double lline[3];
	double rline[3];

	//point from left image
	dll_matrix_product(3,3,3,1,m_FM, lp, rline);

	//point from right image
	dll_matrix_product(1,3,3,3, rp, m_FM, lline);
	
	//calculate the distance from image point to epipolar line, average of left and right 	
	double dotValue = rp[0]*rline[0] + rp[1]*rline[1] + rp[2]*rline[2];
	dis = 0.5* fabs(dotValue)  / 
		( sqrt(rline[0]*rline[0]+rline[1]*rline[1]) + sqrt(lline[0]*lline[0]+lline[1]*lline[1]) ) ;
    
	return dis;
}




int RelativePoseEstimation(char* filename1, char* filename2, double focalLen1, double focalLen2, double& relativeHei)
{    		
	IplImage* pLeft  = cvLoadImage(filename1);
	IplImage* pRight = cvLoadImage(filename1);

	//1. detect feature points 
	CFeatureBase* pFeatDetect = new CSIFT();
	ImgFeature lImageFeature, rImageFeature;
	pFeatDetect->Detect(filename1, lImageFeature);
	pFeatDetect->Detect(filename2, rImageFeature);

	//2. matching 
	CMatchBase* pMatch = new CKNNMatch();
	PairMatchRes mp;
	pMatch->Match(lImageFeature, rImageFeature, mp);
	int nMatch = mp.matchs.size();
	printf("Matching Point Number: %d \n", nMatch);

	//3. relative pose 
	vector<Point2DDouble> lpts,rpts;
	nMatch = mp.matchs.size();
	for(int i=0; i<nMatch; i++)
	{
		int li = mp.matchs[i].l;
		int ri = mp.matchs[i].r;
		Point2DDouble pl,pr;
		pl.p[0] = lImageFeature.featPts[li].cx;
		pl.p[1] = lImageFeature.featPts[li].cy;
		lpts.push_back(pl);
		pr.p[0] = rImageFeature.featPts[ri].cx;
		pr.p[1] = rImageFeature.featPts[ri].cy;
		rpts.push_back(pr);
	}

	CRelativePoseBase* pRP = new CEstimatePose5Point();
	CameraPara cam1, cam2;
	//intrinsic parameters
	cam1.focalLen = focalLen1; //pLeft->width*0.5;
	cam2.focalLen = focalLen2; //pRight->width*0.5;
	pRP->EstimatePose(lpts, rpts, cam1, cam2 );


	//output
	printf("focalLen: %lf \n", cam2.focalLen);
	printf("Rotation Angle x:%lf y:%lf z:%lf \n", cam2.ax, cam2.ay, cam2.az);
	for(int j=0; j<3; j++)
	{
		for(int i=0; i<3; i++)
		{
			printf("%lf ", cam2.R[j*3+i]);
		}
		printf("\n");
	}
	printf("Translation: \n");
	for(int j=0; j<3; j++)
		printf("%lf ", cam2.T[j]);
	printf("\n\n");
    
	//4. triangulation
    CTriangulateBase* pTriangulate = new CTriangulateCV();
	vector<Point3DDouble> gpts;
	pTriangulate->Triangulate(lpts, rpts, cam1, cam2, gpts);
  
	//bundle adjustment

	vector<double> gz;
	for(int i=0; i<gpts.size(); i++)
	{
		gz.push_back(gpts[i].p[2]);
	}
	sort(gz.begin(),gz.end());	
	int half = gz.size() / 2;
	relativeHei = fabs(gz[half]);
	


	/*
	int nImage = 2;
	char** filenames = f2c(nImage, 256);
	strcpy(filenames[0], filename1);
	strcpy(filenames[1], filename2);

	vector<double> focalLen;
	focalLen.push_back(focalLen1);
	focalLen.push_back(focalLen2);

	//1. feature detection
	//vector<IplImage*> pImages;
	vector<CImageDataBase*> vecImageDataPointer;
	for(int i=0; i<nImage; i++)
	{
		CImageDataBase* pImageData = new CImageFeature();
		//char file[256];
		//strcpy(file, fileVector[i]);
		pImageData->Load(filenames[i]);
		pImageData->DetectPtFeature(SIFT_FEATURE);
		pImageData->SetInitFocus(focalLen[i]);
		vecImageDataPointer.push_back(pImageData);
	}

	//2. matching
	vector<PairMatchRes> vecMatch;
	//CMatchBase* pMatch = new CKNNMatch();
	CMatchBase* pMatch = new CSiftMatch();
	for(int i=0; i<nImage; i++)
		for(int j=i+1; j<nImage; j++)
		{
			printf("Matching %d-%d \n", i,j);
			PairMatchRes mt;
			mt.lId = i;
			mt.rId = j;
			pMatch->Match( vecImageDataPointer[i]->GetImageFeature(), vecImageDataPointer[j]->GetImageFeature(), mt);

			//if( mt.inlierRatio>0.5 && mt.matchs.size()>16 )
			if( mt.matchs.size()>16 )
			{
				vecMatch.push_back(mt);
				printf(" inlier: %lf  match: %d \n", mt.inlierRatio, mt.matchs.size());
			}
		}

	//3. generating tracks
	printf("\n Generating Tracks.... \n");
	CGenerateTracksBase* pGenerateTracks = new CMyGenerateTrack();
	vector<TrackInfo> vecTrack;
	pGenerateTracks->GenerateTracks(vecMatch, vecImageDataPointer, vecTrack);
	PrintTracks(vecTrack, "d:\\tracks.txt");

	if(vecTrack.size()<8)
	{
		printf("Tracks are less than the threshold ! \n");
		return -1;
	}


	//4. bundle adjustment
    CBABase* pBA = new CSBA();
	vector<CameraPara> cameras;
    cameras.resize(nImage);
	//double focalLen = vecImageDataPointer[0]->GetHt();
	for(int i=0; i<cameras.size(); i++)
	{
		if( vecImageDataPointer[i]->IsHasInitFocus() )
		{
			cameras[i].focalLen = vecImageDataPointer[i]->GetInitFocus();
		}
		else
		{
			double focalLen = vecImageDataPointer[i]->GetHt();
			cameras[i].focalLen = focalLen; //vecImageDataPointer[i]->GetHt();
			vecImageDataPointer[i]->SetInitFocus(focalLen);
		}
	}
	pBA->BundleAdjust( cameras.size(), cameras, vecImageDataPointer, vecMatch, vecTrack, "d:\\");

	
	delete pGenerateTracks;
	delete pMatch;
	delete pBA;
	*/

	return 0;
}


int CalculateEssentialMatrix1(double* R, double* T, double* F)
{
	//generate the asymmetric matrix
	double aT[3][3];
    aT[0][0] = 0;		aT[0][1]=-T[2];  aT[0][2]=T[1];
	aT[1][0] = T[2];	aT[1][1]=0;		 aT[1][2]=-T[0];
	aT[2][0] = -T[1];	aT[2][1]=T[0];	 aT[2][2]=0;

	//mult(*aT, R, F, 3, 3, 3);
	dll_matrix_product(3,3,3,3,*aT, R, F);
    
	return 0;
}


/*
int GrdToImg(Point3DDouble grdpt, Point2DDouble& imgpt, CameraPara cam)
{
	//int ht, wd;
	double focal = cam.focalLen;
	double k1    = cam.k1;
	double k2    = cam.k2;
	double omiga = cam.ax;
	double phi   = cam.ay;
	double kapa  = cam.az;

	double R[9];
	GenerateRMatrixDirect(omiga, phi, kapa, R);

	double t[3];
	t[0] = cam.t[0];
	t[1] = cam.t[1];
	t[2] = cam.t[2];

	double ix1,iy1;
	double gx = grdpt.p[0];
	double gy = grdpt.p[1];
	double gz = grdpt.p[2];

	GrdToImgWithDistort(gx, gy, gz, &ix1, &iy1, R, t, focal, 0.0, 0.0, k1, k2);

	imgpt.p[0] = ix1;
	imgpt.p[1] = iy1;

	return 0;
}*/