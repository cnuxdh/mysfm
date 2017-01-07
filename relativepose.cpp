
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
#include "defines.h"
//#include "triangulate.h"
//#include "fmatrix.h"

#include "baselib.h"


#include <algorithm>
using namespace std;



// model: R(X-T), from ground point to image
int GrdToImg(Point3DDouble gp, Point2DDouble& ip, CameraPara cam, bool explicit_camera_center)
{
	const double epsilon = 1e-18;
	double p[3];
	double res[3];

	if( cam.camtype == PerspectiveCam )
	{
		double f = cam.focus;
		double k1 = cam.k1;
		double k2 = cam.k2;

		p[0] = gp.p[0] - cam.t[0];
		p[1] = gp.p[1] - cam.t[1];
		p[2] = gp.p[2] - cam.t[2];
		mult(cam.R, p, res, 3, 3, 1);

		//res[0] = R[0]*p[0]+R[1]*p[1]+R[2]*p[2];
		//res[1] = R[3]*p[0]+R[4]*p[1]+R[5]*p[2];
		//res[2] = R[6]*p[0]+R[7]*p[1]+R[8]*p[2];

		double x0 = 0;
		double y0 = 0;
		double cx =  res[0]/(res[2]+epsilon)*(-f) + x0;
		double cy =  res[1]/(res[2]+epsilon)*(-f) + y0;

		double dx = cx;
		double dy = cy;
		Undistort(cx, cy, dx, dy, k1, k2);

		ip.p[0] = dx;
		ip.p[1] = dy;
	}
	else if(cam.camtype == PanoramCam)
	{
		p[0] = gp.p[0] - cam.t[0];
		p[1] = gp.p[1] - cam.t[1];
		p[2] = gp.p[2] - cam.t[2];
		mult(cam.R, p, res, 3, 3, 1);

		double radius = double(cam.cols) / (2*PI);
		double ix,iy;

		GrdToSphere(res[0], res[1], res[2], radius, ix, iy);

		ip.p[0] = ix;
		ip.p[1] = iy;
	}

	return 0;
}




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

void GetIntrinsics(const camera_params_t &camera, double *K) 
{
	if (!camera.known_intrinsics) 
	{
		K[0] = camera.f;  K[1] = 0.0;       K[2] = 0.0;
		K[3] = 0.0;       K[4] = camera.f;  K[5] = 0.0;
		K[6] = 0.0;       K[7] = 0.0;       K[8] = 1.0;    
	} 
	else
	{
		memcpy(K, camera.K_known, 9 * sizeof(double));
	}
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


/* Triangulate two points
   explicit_camera_centers:  true: transformation is R(X-T); false: transformation is RX+T
   in_front: ?
   angle:  ?
*/ 
v3_t PtTriangulate(v2_t p, v2_t q, 
				 camera_params_t c1, camera_params_t c2, 
				 double &proj_error, bool &in_front, double &angle,
				 bool explicit_camera_centers)
{
	double K1[9], K2[9];
	double K1inv[9], K2inv[9];

	GetIntrinsics(c1, K1);
	GetIntrinsics(c2, K2);

	dll_matrix_invert(3, K1, K1inv);
	dll_matrix_invert(3, K2, K2inv);

	/* Set up the 3D point */
	// EDIT!!!
	double proj1[3] = { Vx(p), Vy(p), -1.0 };
	double proj2[3] = { Vx(q), Vy(q), -1.0 };

	double proj1_norm[3], proj2_norm[3];

	dll_matrix_product(3, 3, 3, 1, K1inv, proj1, proj1_norm);
	dll_matrix_product(3, 3, 3, 1, K2inv, proj2, proj2_norm);

	v2_t p_norm = dll_v2_new(proj1_norm[0] / proj1_norm[2],
		proj1_norm[1] / proj1_norm[2]);

	v2_t q_norm = dll_v2_new(proj2_norm[0] / proj2_norm[2],
		proj2_norm[1] / proj2_norm[2]);

	/* Undo radial distortion */
	p_norm = UndistortNormalizedPoint(p_norm, c1);
	q_norm = UndistortNormalizedPoint(q_norm, c2);

	/* Compute the angle between the rays */
	angle = ComputeRayAngle(p, q, c1, c2);

	/* Triangulate the point */
	v3_t pt;
	if (!explicit_camera_centers) 
	{
		pt = dll_triangulate(p_norm, q_norm, c1.R, c1.t, c2.R, c2.t, &proj_error);
	}
	else 
	{
		double t1[3];
		double t2[3];

		/* Put the translation in standard form */
		dll_matrix_product(3, 3, 3, 1, c1.R, c1.t, t1);
		dll_matrix_scale(3, 1, t1, -1.0, t1);
		dll_matrix_product(3, 3, 3, 1, c2.R, c2.t, t2);
		dll_matrix_scale(3, 1, t2, -1.0, t2);

		pt = dll_triangulate(p_norm, q_norm, c1.R, t1, c2.R, t2, &proj_error);
	}

	proj_error = (c1.f + c2.f) * 0.5 * sqrt(proj_error * 0.5);

	/* Check cheirality */
	bool cc1 = CheckCheirality(pt, c1);
	bool cc2 = CheckCheirality(pt, c2);

	in_front = (cc1 && cc2);

	return pt;
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
	camera1.f = cam1.focus;
	camera2.f = cam2.focus;
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
	memcpy(cam1.t, camera1.t, sizeof(double)*3);
	memcpy(cam2.R, camera2.R, sizeof(double)*9);
	memcpy(cam2.t, camera2.t, sizeof(double)*3);
   
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

	double K1[9],K2[9],R[9],t[3],R0[9],t0[3];
	memset(K1, 0, sizeof(double)*9);
	memset(K2, 0, sizeof(double)*9);
	memset(R0, 0, sizeof(double)*9);
	memset(t0, 0, sizeof(double)*3);
	memset(R, 0, sizeof(double)*9);
	memset(t, 0, sizeof(double)*3);

	double focus1 = cam1.focus; //947.2;
	double focus2 = cam2.focus; //947.2;
	K1[0] = focus1;  K1[4] = focus1; K1[8] = 1;
	K2[0] = focus2;  K2[4] = focus2; K2[8] = 1;
	R0[0] = R0[4] = R0[8] = 1;    

	memcpy(cam1.R, R0, sizeof(double)*9);
	memcpy(cam1.t, t0, sizeof(double)*3);

	int num_inliers =  EstimatePose5Point(lPts, rPts, 512, 4,  K1, K2, R, t);
	printf("number of inliers : %d \n", num_inliers);

	//decompose the three angles around x,y,z axis
	cam2.ax = atan( R[5]/R[8] )/PI*180; 
	cam2.ay = asin( -R[2] )/PI*180;
	cam2.az = atan( R[1]/R[0])/PI*180;

	double nt[3];
	//bool initialized = false;
	//if (!initialized) 
	{
		memcpy(cam2.R, R, sizeof(double) * 9);
		
		dll_matrix_transpose_product(3, 3, 3, 1, R, t, nt);
		dll_matrix_scale(3, 1, nt, -1.0, nt);
		memcpy(cam2.t,nt,sizeof(double)*3);
	}	
   
	m_bIsExplicit = true;

	/*
	//print results
#ifdef _DEBUG
	
	FILE* fp = fopen("c:\\temp\\myRelativePose.txt", "w");
	fprintf(fp, "%d %d \n", 0, 1);
	fprintf(fp, "focus: %lf \n", cam2.focus);
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
		fprintf(fp,"%lf ", cam2.t[j]);
	fprintf(fp, "\n");
	fclose(fp);
#endif
	*/

	////output
	//printf("Focus: %lf \n", cam2.focus);
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
	//	printf("%lf ", cam2.t[j]);
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
		x = wd*0.5 + x;
		y = ht*0.5 - y;
		double gx,gy,gz;
		SphereTo3D(x, y, radius, gx, gy, gz);
		Point3DDouble p3;
		p3.p[0] = gx / radius;
		p3.p[1] = gy / radius;
		p3.p[2] = gz / radius;
		pl.push_back(p3);


		x = rPts[i].p[0];
		y = rPts[i].p[1];
		//from image center to top-left
		x = wd*0.5 + x;
		y = ht*0.5 - y;
		SphereTo3D(x, y, radius, gx, gy, gz);
		p3.p[0] = gx / radius;
		p3.p[1] = gy / radius;
		p3.p[2] = gz / radius;
		pr.push_back(p3);
	}
	
	int num_trials = 512;
	double threshold = 2.5;
	double R[9];
	double t[3];

	vector<double> residual;
	residual.resize(pl.size());
	EstimatePose5Point_Pano(pl, pr, radius, num_trials, threshold, R, t, residual);

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
	mult(Rt, t, gt, 3, 3, 1);
	for(int i=0; i<3; i++)
	{
		cam1.t[i] = 0;
		cam2.t[i] = -gt[i];
	}

	m_bIsExplicit = true;
	
	//save eular angle
	double ea[3];
	rot2eular(R, ea);
	cam1.ax = cam1.ay = cam1.az = 0;
	cam2.ax = ea[0];
	cam2.ay = ea[1];
	cam2.az = ea[2];

	printf("relative pose eular angle: %lf %lf %lf \n", ea[0], ea[1], ea[2]);
	printf("relative pose translation: %lf %lf %lf \n", cam2.t[0], cam2.t[1], cam2.t[2]);
	
	return 0;
}

//////////////////////////////////////////////////////////////////////////
CTriangulateCV::CTriangulateCV()
{

}
CTriangulateCV::~CTriangulateCV()
{

}
void CTriangulateCV::Triangulate(std::vector<Point2DDouble> lPts, std::vector<Point2DDouble> rPts, 
								 CameraPara cam1, CameraPara cam2, std::vector<Point3DDouble> &gps)
{
	//printf(" \n Triangulation .... \n");
	//assert( lPts.size() == rPts.size() );
	
	camera_params_t camera1,camera2;
	InitializeCameraParams(camera1);
	InitializeCameraParams(camera2);
    
	camera1.f = cam1.focus;
	camera1.known_intrinsics = false;
	memcpy(camera1.R, cam1.R, sizeof(double)*9);
	memcpy(camera1.t, cam1.t, sizeof(double)*3);

	camera2.f = cam2.focus;
	camera2.known_intrinsics = false;
	memcpy(camera2.R, cam2.R, sizeof(double)*9);
	memcpy(camera2.t, cam2.t, sizeof(double)*3);

	for(int i=0; i<lPts.size(); i++)
	{
		v2_t p = dll_v2_new(lPts[i].p[0], lPts[i].p[1]);
		v2_t q = dll_v2_new(rPts[i].p[0], rPts[i].p[1]);
        
		double error = 0.0;
		bool   in_front = true;
		double angle = 0.0;
         
		v3_t gp = PtTriangulate( p, q, camera1, camera2, error, in_front, angle, true);
		
		//the projection error is beyond the threshold
		if (error > 4) 
		{
			//printf(" skipping point\n");
			continue;
		}		

		Point3DDouble p3;
		p3.p[0] = gp.p[0]; p3.p[1] = gp.p[1]; p3.p[2] = gp.p[2];
        p3.extra = i;

		gps.push_back(p3);

		//printf("%lf %lf ")
	}
}

void CTriangulateCV::Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
								 CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, 
								 vector<double>& errorArray)
{
	//printf(" \n Triangulation .... \n");
	//assert( lPts.size() == rPts.size() );
	/*
	camera_params_t camera1,camera2;
	InitializeCameraParams(camera1);
	InitializeCameraParams(camera2);

	camera1.f = cam1.focus;
	camera1.known_intrinsics = false;
	memcpy(camera1.R, cam1.R, sizeof(double)*9);
	memcpy(camera1.t, cam1.t, sizeof(double)*3);

	camera2.f = cam2.focus;
	camera2.known_intrinsics = false;
	memcpy(camera2.R, cam2.R, sizeof(double)*9);
	memcpy(camera2.t, cam2.t, sizeof(double)*3);
	*/

	for(int i=0; i<lPts.size(); i++)
	{
		vector<CameraPara> cams;
		vector<Point2DDouble> pts;

		cams.push_back(cam1);
		cams.push_back(cam2);
		pts.push_back( lPts[i] );
		pts.push_back( rPts[i] );
		
		/*
		v2_t p = dll_v2_new(lPts[i].p[0], lPts[i].p[1]);
		v2_t q = dll_v2_new(rPts[i].p[0], rPts[i].p[1]);
		double error = 0.0;
		bool   in_front = true;
		double angle = 0.0;
		v3_t gp = PtTriangulate( p, q, camera1, camera2, error, in_front, angle, true);
		
		Point3DDouble p3;
		p3.p[0] = gp.p[0]; p3.p[1] = gp.p[1]; p3.p[2] = gp.p[2];
		p3.extra = i;
		gps.push_back(p3);
		errorArray.push_back(error);
		*/

		Point3DDouble p3;
		double error;
		Triangulate(pts, cams, p3, true, error);
		p3.extra = i;

		gps.push_back(p3);
		errorArray.push_back(error); 
	}
}


//default model: RX+T, single point
void CTriangulateCV::Triangulate(vector<Point2DDouble> pts,  
								 vector<CameraPara> cams, 
							 	 Point3DDouble& gps,
								 bool explicit_camera_centers,
								 double& ferror)
{
	
	bool estimate_distortion = false;

	int num_views = (int) cams.size();

	camera_params_t *cameras = new camera_params_t[num_views];

	for(int i=0; i<num_views; i++)
	{
		InitializeCameraParams(cameras[i]);
		cameras[i].f = cams[i].focus;
		cameras[i].known_intrinsics = false; //we do not use the input intrisics!!!!
		memcpy(cameras[i].R, cams[i].R, sizeof(double)*9);
		memcpy(cameras[i].t, cams[i].t, sizeof(double)*3);
	}
	
	v2_t *pv = new v2_t[num_views];
	double *Rs = new double[9 * num_views];
	double *ts = new double[3 * num_views];

	for (int i = 0; i < num_views; i++) 
	{
		camera_params_t *cam = NULL;

		double cx = pts[i].p[0];
		double cy = pts[i].p[1];

		if(cams[i].camtype == PerspectiveCam)
		{
			//undistort the image points, added by xiedonghai, 2016.11.18
			double dx,dy;
			double k1 = cams[i].k1;
			double k2 = cams[i].k2;
			Undistort(cx, cy, dx, dy, k1, k2);
			double p3[3] = { dx, dy, 1.0 };
			double K[9], Kinv[9];
			GetIntrinsics(cameras[i], K);
			dll_matrix_invert(3, K, Kinv);
			double p_n[3];
			dll_matrix_product(3, 3, 3, 1, Kinv, p3, p_n);
			// EDIT!!!
			pv[i] = dll_v2_new(-p_n[0], -p_n[1]);
			//pv[i] = UndistortNormalizedPoint(pv[i], cameras[i]);
		}
		else if( cams[i].camtype == PanoramCam )
		{
			int ht = cams[i].rows;
			int wd = cams[i].cols;
			double radius = double(wd) / (2*PI);

			//from spherical panorama image point to 3D point
			double gx,gy,gz;
			double x = cx + wd*0.5;
			double y = ht*0.5 - cy;
			SphereTo3D(x,y,radius,gx,gy,gz);
			//normalized the 3D point
			pv[i].p[0] = gx/(gz+MINIMAL_VALUE);
			pv[i].p[1] = gy/(gz+MINIMAL_VALUE);
		}

		//
		cam = cameras + i;
		memcpy(Rs + 9 * i, cam->R, 9 * sizeof(double));
		if (!explicit_camera_centers) 
		{
			memcpy(ts + 3 * i, cam->t, 3 * sizeof(double));
		}
		else 
		{
			dll_matrix_product(3, 3, 3, 1, cam->R, cam->t, ts + 3 * i);
			dll_matrix_scale(3, 1, ts + 3 * i, -1.0, ts + 3 * i);
		}
	}

	ferror = 0;

	//model: RX+T
	v3_t pt = dll_triangulate_n(num_views, pv, Rs, ts, &ferror);


	//calculate the average error 
	ferror = 0.0;
	for (int i = 0; i < num_views; i++)
	{
		double cx = pts[i].p[0];
		double cy = pts[i].p[1];

		FeatPoint key; // = imageData[image_idx]->GetKeyPoint(key_idx);
		key.cx = cx;
		key.cy = cy;

		/*
		v2_t pr = dll_sfm_project_final(cameras + i, pt, 
			explicit_camera_centers ? 1 : 0,
			estimate_distortion ? 1 : 0);
			*/

		Point3DDouble gpt;
		Point2DDouble ipt;
		gpt.p[0] = pt.p[0];
		gpt.p[1] = pt.p[1];
		gpt.p[2] = pt.p[2];
		
		GrdToImg(gpt, ipt, cams[i]);
		
		double dx = ipt.p[0] - key.cx; //key.m_x;
		double dy = ipt.p[1] - key.cy; //key.m_y;

		ferror += dx * dx + dy * dy;
	}
	ferror = sqrt(ferror / num_views);

	gps.p[0] = pt.p[0];
	gps.p[1] = pt.p[1];
	gps.p[2] = pt.p[2];

	delete [] pv;
	delete [] Rs;
	delete [] ts;

	delete [] cameras;
}


CTriangulatePano::CTriangulatePano()
{

}
CTriangulatePano::~CTriangulatePano()
{

}

void CTriangulatePano::Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
	CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps)
{

}

void CTriangulatePano::Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
	CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, vector<double>& errorArray)
{
	int ht = cam1.rows;
	int wd = cam1.cols;
	double radius = (double)(wd) / (2*PI);

	bool explicit_camera_centers = true;

	double R1[9];
	double t1[3];
	double R2[9];
	double t2[3];

	memcpy(R1, cam1.R, sizeof(double)*9);
	memcpy(t1, cam1.t, sizeof(double)*3);
	memcpy(R2, cam2.R, sizeof(double)*9);
	memcpy(t2, cam2.t, sizeof(double)*3);

	//from R(X-T) to RX+T
	if (explicit_camera_centers) 
	{
		/* Put the translation in standard form */
		mult(R1, cam1.t, t1, 3, 3, 1);
		t1[0] = -t1[0];
		t1[1] = -t1[1];
		t1[2] = -t1[2];

		mult(R2, cam2.t, t2, 3, 3, 1);
		t2[0] = -t2[0];
		t2[1] = -t2[1];
		t2[2] = -t2[2];
	}

	for(int i=0; i<lPts.size(); i++)
	{
		//from spherical panorama image point to 3D point
		double gx,gy,gz;
		double x = lPts[i].p[0] + wd*0.5;
		double y = ht*0.5 - lPts[i].p[1];
		SphereTo3D(x,y,radius,gx,gy,gz);
		//normalized the 3D point
		Point2DDouble p;
		p.p[0] = gx/(gz+MINIMAL_VALUE);
		p.p[1] = gy/(gz+MINIMAL_VALUE);


		x = rPts[i].p[0] + wd*0.5;
		y = ht*0.5 - rPts[i].p[1];
		SphereTo3D(x,y,radius,gx,gy,gz);
		Point2DDouble q;
		q.p[0] = gx/(gz+MINIMAL_VALUE);
		q.p[1] = gy/(gz+MINIMAL_VALUE);
		
		double error = 0.0;
		Point3DDouble gp = TriangulatePt(p, q, R1, t1, R2, t2, &error);
				
		gps.push_back(gp);
		errorArray.push_back(error);
	}
}

void CTriangulatePano::Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, 
	Point3DDouble& gps,bool explicit_camera_centers,double& ferror)
{

}

//////////////////////////////////////////////////////////////////////////
CDLTPose::CDLTPose()
{

}

CDLTPose::~CDLTPose()
{

}

// the projection matrix: P=K[R|t]
int CDLTPose::EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, double* K, double* R, double* t)
{	
	double P[12];
	int r = -1;
	int num_points = pt3.size();

	//too few points
	if (num_points < 9)
		return -1;

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

	int res = EstimatePose(pt3, pt2, cam.K, cam.R, cam.t);
	
	//double ratio = 
	cam.focus = (cam.K[0]+cam.K[4])*0.5;
	
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

		//from center to top-left
		x = pt2[i].p[0] + wd*0.5;
		y = ht*0.5 - pt2[i].p[1];

		//calculate the 3D ray light 
		double gx,gy,gz;
		SphereTo3D(x,y,radius,gx,gy,gz);
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
	if (num_points < 9)
		return -1;

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

	double proj_estimation_threshold = 2.0;
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
	memcpy(cam.t, t, sizeof(double)*3);	 

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


Point3DDouble TriangulatePt(Point2DDouble p, Point2DDouble q, 
				  double *R0, double *t0, 
				  double *R1, double *t1, double *error)
{

	v2_t lp,rp;
    
	lp.p[0] = p.p[0];
	lp.p[1] = p.p[1];

	rp.p[0] = q.p[0];
	rp.p[1] = q.p[1];

	v3_t gp = dll_triangulate(lp, rp, R0, t0, R1, t1, error );

	Point3DDouble grdPt;
	grdPt.p[0] = gp.p[0];
	grdPt.p[1] = gp.p[1];
	grdPt.p[2] = gp.p[2];

	return grdPt;
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
	cam1.focus = focalLen1; //pLeft->width*0.5;
	cam2.focus = focalLen2; //pRight->width*0.5;
	pRP->EstimatePose(lpts, rpts, cam1, cam2 );


	//output
	printf("Focus: %lf \n", cam2.focus);
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
		printf("%lf ", cam2.t[j]);
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
	//double focus = vecImageDataPointer[0]->GetHt();
	for(int i=0; i<cameras.size(); i++)
	{
		if( vecImageDataPointer[i]->IsHasInitFocus() )
		{
			cameras[i].focus = vecImageDataPointer[i]->GetInitFocus();
		}
		else
		{
			double focus = vecImageDataPointer[i]->GetHt();
			cameras[i].focus = focus; //vecImageDataPointer[i]->GetHt();
			vecImageDataPointer[i]->SetInitFocus(focus);
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
	double focal = cam.focus;
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