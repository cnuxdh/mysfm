
#ifndef  CV_SIFT_H
#define  CV_SIFT_H

#include "export.hpp"
#include "feature.hpp"
#include "defines.hpp"



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


class DLL_EXPORT CSIFT: public CFeatureBase
{
public:
	CSIFT();
	~CSIFT();
	int Detect(char* filePath, char* featurePath);
	int Detect(char* filePath, ImgFeature& imgFeat);
	int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat);
	int Detect(IplImage* pImage, ImgFeature& imgFeat);
private:

};



class DLL_EXPORT CSIFTFloat: public CFeatureBase
{
public:
	CSIFTFloat();
	~CSIFTFloat();
	int Detect(char* filePath, char* featurePath);
    int Detect(char* filePath, ImgFeature& imgFeat);
	int Detect(char* filePath, ImgFeature& imgFeat, int maxHt);
	int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat);
	int Detect(IplImage* pImage, ImgFeature& imgFeat);
private:

};

class DLL_EXPORT CSIFTPano: public CFeatureBase
{
public:
	CSIFTPano();
	~CSIFTPano();

	int Detect(IplImage* pImage, ImgFeature& imgFeat, int maxHt);
	int Detect(char* filePath, ImgFeature& imgFeat, int maxHt);
};





#endif