

#ifndef  CV_FEATURES_H
#define  CV_FEATURES_H


#include "defines.hpp"
#include "badata.hpp"


#ifdef OPENCV_1X 
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#else
#include "opencv2/core/core.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/highgui/highgui.hpp"
using namespace cv;
#endif



//the maximal image height for feature detection
//#define DETECT_IMAGE_HT 420

//get the image dimension after scale
void GetResizeDimension(int srcHt, int srcWd, int& dstHt, int& dstWd);

//////////////////////////////////////////////////////////////////////////
// feature detection interface 
class DLL_EXPORT CFeatureBase
{
public:
	CFeatureBase(){}
	virtual ~CFeatureBase(){}

	virtual int Detect(char* filePath, char* featurePath){return 0;}
	virtual int Detect(char* filePath, vector<double>& px, vector<double>& py){return 0;}
	virtual int Detect(char* filePath, vector<PtFeature>& featPts){return 0;};
	
	virtual int Detect(char* filePath, ImgFeature& imgFeat){return 0;};
	virtual int Detect(char* filePath, ImgFeature& imgFeat, int maxHt){return 0;};

	virtual int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat){return 0;};

	virtual int Detect(IplImage* pImage, ImgFeature& imgFeat){return 0;};
};



//feature I/O base class
class DLL_EXPORT CFeatureIO
{
public:
	CFeatureIO(){}
	virtual ~CFeatureIO(){}
	
	virtual int WriteBin(ImgFeature& imgFeat){return 0;}
	virtual int WriteTxt(ImgFeature& imgFeat){return 0;}
};




#endif