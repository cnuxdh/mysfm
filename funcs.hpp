
#ifndef FUNCS_HPP
#define FUNCS_HPP

#ifdef OPENCV_1X 
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#else
//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#endif

#include "defines.hpp"
#include "dataBase.hpp"
#include "pos.hpp"



double norm(Point3DDouble p);

double dot(Point3DDouble p1, Point3DDouble p2);

double distanceVec(Point2DDouble p1, Point2DDouble p2);

//return value: degree
double angleOfVector(Point3DDouble p1, Point3DDouble p2);


DLL_EXPORT int	   DrawMatches(char* filename, IplImage* pLeft, IplImage* pRight, 
							vector<Point2DDouble> lpts, vector<Point2DDouble> rpts);

int    CalculateColorHist(IplImage* pImage, int grayStep, vector<double>& hist);

DLL_EXPORT int    InitCamera(CameraPara& cam, POSInfo pos);

DLL_EXPORT int	  CalculateViewField(Point3DDouble T, Point3DDouble leftLight, 
	double maxDistance,  double minDistance,
	Point3DDouble& rightLightBottom, Point3DDouble& rightLightTop );

IplImage* ResizeImage(IplImage* src, int maxHt);

DLL_EXPORT int DrawFeatPt(ImgFeature& featpts, IplImage* pImage);



#endif
