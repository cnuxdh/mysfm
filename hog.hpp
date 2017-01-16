
#ifndef HOG_HPP
#define HOG_HPP

#include "defines.hpp"



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



DLL_EXPORT int CalculateHOG(IplImage* pImage, vector<double>& hog);



#endif