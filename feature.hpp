

#ifndef  CV_FEATURES_H
#define  CV_FEATURES_H

#include "defines.hpp"


//the maximal image height for feature detection
#define DETECT_IMAGE_HT 640

//get the image dimension after scale
void GetResizeDimension(int srcHt, int srcWd, int& dstHt, int& dstWd);


//////////////////////////////////////////////////////////////////////////
// feature detection interface 
class CFeatureBase
{
public:
	CFeatureBase(){}
	virtual ~CFeatureBase(){}

	virtual int Detect(char* filePath, char* featurePath){return 0;}
	virtual int Detect(char* filePath, vector<double>& px, vector<double>& py){return 0;}
	virtual int Detect(char* filePath, vector<PtFeature>& featPts){return 0;};
	virtual int Detect(char* filePath, ImgFeature& imgFeat){return 0;};
	virtual int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat){return 0;};
};



//feature I/O
class CFeatureIO
{
public:
	CFeatureIO(){}
	virtual ~CFeatureIO(){}
	
	virtual int WriteBin(ImgFeature& imgFeat){return 0;}
	virtual int WriteTxt(ImgFeature& imgFeat){return 0;}
};




#endif