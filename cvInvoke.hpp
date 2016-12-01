
#ifndef CV_INVOKE_H
#define CV_INVOKE_H

#include <vector>
using namespace std;


#include"sift.hpp"
#include"dataBase.hpp"
#include"register.hpp"



int DLL_EXPORT DetectFileFeaturePts(char** filenames, int nFile, char* outpath);

int DLL_EXPORT DetectFileFeaturePts(char** filenames, int nFile, vector<ImgFeature>& imgFeatures, int maxHt);

//for single image
int DLL_EXPORT DetectFileFeaturePts(char* filenames, ImgFeature& imgFeatures, int maxHt);


int DLL_EXPORT MatchImageFiles(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& matchRes);


DLL_EXPORT int dll_EstimatePose5Point_Pano( vector<Point3DDouble>& p1, 
									vector<Point3DDouble>& p2,
									double radius,
									int num_trials, double threshold, 
									double *R, double *t, vector<double>& residual);


DLL_EXPORT Point3DDouble dll_TriangulatePt(Point2DDouble p, Point2DDouble q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error);


#endif





