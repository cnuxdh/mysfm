
#ifndef PANOBA_HPP
#define PANOBA_HPP

#include "ba.hpp"



class DLL_EXPORT CPanoBA: public CBABase
{
public:
	CPanoBA();
	~CPanoBA();

	int RunSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
		vector<ImgFeature> imageFeatures,  vector<int> cameraIDOrder,
		vector<CameraPara>& cameras){return 0;}

	//
	int BundleAdjust( int numCameras, vector<CameraPara>& cameras,vector<ImgFeature> imageFeatures, 
		vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir){return 0;}


	//new interface, including more input parameters
	int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
		vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir){return 0;}


	int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<ImgFeature>& imageFeatures, 
		vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks);

private:
};





#endif