
#ifndef CV_INVOKE_H
#define CV_INVOKE_H

#include <vector>
using namespace std;

#include "dataBase.hpp"
#include "defines.hpp"

//
//#include <opencv2/core/core.hpp>
//#include <opencv2/highgui/highgui.hpp>
//#include <opencv2/features2d/features2d.hpp>
//#include <opencv2/imgproc/imgproc.hpp>
//using namespace cv;



int DLL_EXPORT DetectFileFeaturePts(char** filenames, int nFile, char* outpath);

//multiple images
int DLL_EXPORT DetectFileFeaturePts(char** filenames, int nFile, 
	vector<ImgFeature>& imgFeatures, int maxHt, double& dProgress);

//for single image
int DLL_EXPORT DetectFileFeaturePts(char* filenames, ImgFeature& imgFeatures, int maxHt);


int DLL_EXPORT MatchImageFiles(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& matchRes,
	CameraType camtype, int matchSteps);


DLL_EXPORT int dll_EstimatePose5Point_Pano( vector<Point3DDouble>& pl, 
									vector<Point3DDouble>& pr,
									double radius,
									int num_trials, double threshold, 
									double *R, double *t, vector<double>& residual);


DLL_EXPORT double dll_CalculatePanoEpipolarError(double* em, Point3DDouble lp, Point3DDouble rp, double radius);



//interface for relative pose estimation
DLL_EXPORT int dll_EstimatePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts,
				CameraPara& cam1, CameraPara& cam2, CameraType camtype );


//interface for DLT algorithm
DLL_EXPORT int dll_DLT(vector<Point3DDouble>& grds, vector<Point2DDouble>& projs,
						CameraPara& cam, CameraType camType);


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


DLL_EXPORT int dll_GenerateRainbowMapping(vector<int>& r, vector<int>& g, vector<int>& b);



//DLL_EXPORT int dll_DetectFeatPts(char** filenames, int nfile,
//					vector<CImageDataBase*>& vecImageDataPointer, double& dProgress);

//ORB feature detection
DLL_EXPORT int dll_GenerateORBFeature(string filepath, vector<KeyPoint>& Keys, Mat& Descriptors);




#endif





