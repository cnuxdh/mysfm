


#ifndef  CV_SURF_H
#define  CV_SURF_H

#include "export.hpp"
#include "feature.hpp"
#include "defines.hpp"



class DLL_EXPORT CSURF: public CFeatureBase
{
public:
	CSURF();
	~CSURF();
	int Detect(char* filePath, char* featurePath);
	int Detect(char* filePath, vector<double>& px, vector<double>& py);
	int Detect(char* filePath, vector<PtFeature>& featPts);
	int Detect(char* filePath, ImgFeature& imgFeat);
private:
};

class DLL_EXPORT CSURFKeyPoint: public CFeatureBase
{
public:
	CSURFKeyPoint();
	~CSURFKeyPoint();
	//int Detect(char* filePath, char* featurePath);
	//int Detect(char* filePath, vector<double>& px, vector<double>& py);
	int Detect(char* filePath, vector<PtFeature>& featPts);
	//int Detect(char* filePath, ImgFeature& imgFeat);
private:
};









#endif