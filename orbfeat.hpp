
#ifndef ORB_FEAT_HPP
#define ORB_FEAT_HPP


#include "export.hpp"
#include "feature.hpp"




class DLL_EXPORT CORBFeat : public CFeatureBase
{
public:
	CORBFeat();
	~CORBFeat();
	//int Detect(char* filePath, char* featurePath);
	//int Detect(char* filePath, ImgFeature& imgFeat);
	//int Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat);
	//int Detect(IplImage* pImage, ImgFeature& imgFeat);
	int Detect(char* filePath, ImgFeature& imgFeat, int maxHt);
	int Detect(Mat& image, ImgFeature& imgFeat, int maxHt);
	
	//vector<KeyPoint> mKeys;
	//Mat mDescriptors;

private:

};


class DLL_EXPORT CORBFrame
{
public:
	CORBFrame();
	~CORBFrame();

	int Detect(string filePath, int maxHt);

	//save the feature
	vector<KeyPoint> mKeys;
	Mat mDescriptors;
};




#endif


