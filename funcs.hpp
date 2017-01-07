
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
#endif



#include "defines.hpp"

double distanceVec(Point2DDouble p1, Point2DDouble p2);

//return value: degree
double angleOfVector(Point3DDouble p1, Point3DDouble p2);


int DrawMatches(char* filename, IplImage* pLeft, IplImage* pRight, 
	vector<Point2DDouble> lpts, vector<Point2DDouble> rpts);


#endif