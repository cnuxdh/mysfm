
#ifndef  CV_SIFT_H
#define  CV_SIFT_H

#include "export.hpp"
#include "feature.hpp"
#include "defines.hpp"


class DLL_EXPORT CSIFT: public CFeatureBase
{
public:
	CSIFT();
	~CSIFT();
	int Detect(char* filePath, char* featurePath);
	int Detect(char* filePath, ImgFeature& imgFeat);
	int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat);
private:
};



class DLL_EXPORT CSIFTFloat: public CFeatureBase
{
public:
	CSIFTFloat();
	~CSIFTFloat();
	int Detect(char* filePath, char* featurePath);
    int Detect(char* filePath, ImgFeature& imgFeat);
	int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat);
private:

};





#endif