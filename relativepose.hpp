
#ifndef CV_RELATIVE_POSE_H
#define CV_RELATIVE_POSE_H


#include "register.hpp"
#include "defines.hpp"


//sfm-driver
#include "sfm.h"


void FixIntrinsics(double *P, double *K, double *R, double *t); 
void GetIntrinsics(const camera_params_t &camera, double *K); 
bool CheckCheirality(v3_t p, const camera_params_t &camera);



//////////////////////////////////////////////////////////////////////////
//relative pose 
class DLL_EXPORT CRelativePoseBase
{
public:
	CRelativePoseBase(){}
	virtual ~CRelativePoseBase(){}
	virtual int EstimatePose( PairMatchRes pairMatches, ImgFeature lImageFeat, ImgFeature rImageFeat, CameraPara& cam1, CameraPara& cam2 ){return 0;}
	virtual int EstimatePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, CameraPara& cam1, CameraPara& cam2 ){return 0;}
};

class DLL_EXPORT CEstimatePose5Point: public CRelativePoseBase
{
public:
	CEstimatePose5Point();
	~CEstimatePose5Point();
	int EstimatePose( PairMatchRes pairMatches, ImgFeature lImageFeat, ImgFeature rImageFeat, CameraPara& cam1, CameraPara& cam2 );
	int EstimatePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, CameraPara& cam1, CameraPara& cam2 );
};




//////////////////////////////////////////////////////////////////////////
//recover the pose from 3D to 2D
class DLL_EXPORT CPoseEstimationBase
{
public:
	CPoseEstimationBase(){}
	virtual ~CPoseEstimationBase(){}
	virtual int EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, double* K, double* R, double* t){return 0;}
	virtual int EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, CameraPara& cam){return 0;}
};


class DLL_EXPORT CDLTPose: public CPoseEstimationBase
{
public:
	CDLTPose();
	~CDLTPose();
	int EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, double* K, double* R, double* t);
	int EstimatePose(vector<Point3DDouble> pt3, vector<Point2DDouble> pt2, CameraPara& cam);
};
/////////////////////////////////////////////////////////////////////////


/////////////////////////////////////////////////////////////////////////
/* calculate the projection of 3D point
   input: 
		grdpt: ground point
		cam:   camera parameters
   output:
		imgpt: normalized point, origin is the center of image
*/
int GrdToImg(Point3DDouble grdpt, Point2DDouble& imgpt, CameraPara cam);






//////////////////////////////////////////////////////////////////////////
//recover 3D coordinates of point
class DLL_EXPORT CTriangulateBase
{
public:
	CTriangulateBase(){}
	virtual ~CTriangulateBase(){}
	
	//for track sequence with two projections
	virtual void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps){} 
	virtual void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, vector<double>& errorArray){} 
		
	//for one track with multiple projections
	virtual void Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, Point3DDouble& gps,bool explicit_camera_centers,double& ferror){} 

};

class DLL_EXPORT CTriangulateCV: public CTriangulateBase
{
public:
	CTriangulateCV();
    ~CTriangulateCV();
	void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps);
	void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, vector<double>& errorArray);
	void Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, Point3DDouble& gps,bool explicit_camera_centers,double& ferror);
};
//////////////////////////////////////////////////////////////////////////



/**********************  Fundamental Matrix *****************************/
//////////////////////////////////////////////////////////////////////////
class DLL_EXPORT CFMBase
{
public:
	CFMBase(){}
	virtual ~CFMBase(){}
	
	//
	virtual int    Estimate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, vector<double>& fm){return 0;}
	
	//flag: 0-left image, 1-right image
	virtual void   CalculateEpipolarLine(Point2DDouble pt, int flag, double& a, double& b, double& c){}
	
	//calculate the distance from the image point to epipolar lines
	virtual double CalculateError(Point2DDouble lpt, Point2DDouble rpt){return 0;}
};


class DLL_EXPORT COpenCVFM: public CFMBase
{
public:
	COpenCVFM();
	~COpenCVFM();
	int    Estimate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, vector<double>& fm);
	void   CalculateEpipolarLine(Point2DDouble pt, int flag, double& a, double& b, double& c);
	double CalculateError(Point2DDouble lpt, Point2DDouble rpt);

private:
	//Mat m_FM;
	double m_FM[9];
};
//////////////////////////////////////////////////////////////////////////



/**********************  Homography Matrix *****************************/







DLL_EXPORT int CalculateExplicitT(double* R, double* T, double* explicitT);

/*
   relative pose estimation for panorama
*/
DLL_EXPORT int EstimatePose5Point_Pano( vector<Point3DDouble>& p1, 
							 vector<Point3DDouble>& p2,
							 double radius,
							 int num_trials, double threshold, 
							 double *R, double *t, vector<double>& residual);


DLL_EXPORT Point3DDouble TriangulatePt(Point2DDouble p, Point2DDouble q, 
							 double *R0, double *t0, 
							 double *R1, double *t1, double *error);


DLL_EXPORT int RelativePoseEstimation(char* filename1, char* filename2, double focalLen1, double focalLen2, double& relativeHei);


//calculate the fundamental matrix using: F=[T]x.R
DLL_EXPORT int CalculateEssentialMatrix1(double* R, double* T, double* F);




#endif