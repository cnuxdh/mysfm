
#ifndef CV_INVOKE_H
#define CV_INVOKE_H

#include <vector>
using namespace std;



#include"dataBase.hpp"
#include "defines.hpp"


int DLL_EXPORT DetectFileFeaturePts(char** filenames, int nFile, char* outpath);

int DLL_EXPORT DetectFileFeaturePts(char** filenames, int nFile, vector<ImgFeature>& imgFeatures, int maxHt);

//for single image
int DLL_EXPORT DetectFileFeaturePts(char* filenames, ImgFeature& imgFeatures, int maxHt);


int DLL_EXPORT MatchImageFiles(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& matchRes);


DLL_EXPORT int dll_EstimatePose5Point_Pano( vector<Point3DDouble>& pl, 
									vector<Point3DDouble>& pr,
									double radius,
									int num_trials, double threshold, 
									double *R, double *t, vector<double>& residual);



//interface for relative pose estimation
DLL_EXPORT int dll_EstimatePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts,
				CameraPara& cam1, CameraPara& cam2 );


DLL_EXPORT Point3DDouble dll_TriangulatePt(Point2DDouble p, Point2DDouble q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error);

DLL_EXPORT void dll_Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, Point3DDouble& gps,
				bool explicit_camera_centers,double& ferror);


DLL_EXPORT void dll_GenerateRMatrix(double omiga, double phi, double kapa, double* R);


DLL_EXPORT void dll_GrdToImg(double gx, double gy, double gz, double* ix, double* iy, 
	double* R, double* Ts, double f, double x0, double y0,
	int ht, int wd);

DLL_EXPORT int dll_EstimatePose5Point_Pano( vector<Point3DDouble>& p1, 
	vector<Point3DDouble>& p2,
	double radius,
	int num_trials, double threshold, 
	double *R, double *t, vector<double>& residual);

#endif





