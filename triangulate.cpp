
#include "triangulate.hpp"
#include "cali.hpp"
#include "baselib.h"
#include "distortion.hpp"
#include "CalcAngle.h"
#include "panorama.hpp"

//sfm lib
#include "sfm.h"

//corelib
#include "matrix.h"


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
	//angle = ComputeRayAngle(p, q, c1, c2);

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

	in_front = true;

	/* Check cheirality */
	/*
	bool cc1 = CheckCheirality(pt, c1);
	bool cc2 = CheckCheirality(pt, c2);
	in_front = (cc1 && cc2);
	*/

	return pt;
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
    
	camera1.f = cam1.focalLen;
	camera1.known_intrinsics = false;
	memcpy(camera1.R, cam1.R, sizeof(double)*9);
	memcpy(camera1.t, cam1.T, sizeof(double)*3);

	camera2.f = cam2.focalLen;
	camera2.known_intrinsics = false;
	memcpy(camera2.R, cam2.R, sizeof(double)*9);
	memcpy(camera2.t, cam2.T, sizeof(double)*3);

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

	int num_views = pts.size();// (int)cams.size();

	camera_params_t *cameras = new camera_params_t[num_views];

	for(int i=0; i<num_views; i++)
	{
		InitializeCameraParams(cameras[i]);
		cameras[i].f = cams[i].focalLen;
		cameras[i].known_intrinsics = false; //we do not use the input intrisics!!!!
		memcpy(cameras[i].R, cams[i].R, sizeof(double)*9);
		memcpy(cameras[i].t, cams[i].T, sizeof(double)*3);
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
			double x = cx; //;cx + wd*0.5;
			double y = cy; //ht*0.5 - cy;
			SphereTo3D_center(x, y, radius, gx, gy, gz);
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
		//printf("%lf %lf \n", dx, dy);

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
	memcpy(t1, cam1.T, sizeof(double)*3);
	memcpy(R2, cam2.R, sizeof(double)*9);
	memcpy(t2, cam2.T, sizeof(double)*3);

	//from R(X-T) to RX+T
	if (explicit_camera_centers) 
	{
		/* Put the translation in standard form */
		mult(R1, cam1.T, t1, 3, 3, 1);
		t1[0] = -t1[0];
		t1[1] = -t1[1];
		t1[2] = -t1[2];

		mult(R2, cam2.T, t2, 3, 3, 1);
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