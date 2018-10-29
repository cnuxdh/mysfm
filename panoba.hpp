
#ifndef PANOBA_HPP
#define PANOBA_HPP

#include "ba.hpp"

#include "ceresba.hpp"


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


class DLL_EXPORT CPanoBAWithPos 
{
public:
	CPanoBAWithPos();
	~CPanoBAWithPos();

	int SetCenter(double ax, double ay, double az);

	int BundleAdjust(int numCameras, vector<CameraPara>& cameras,
					vector<ImgFeature>& imageFeatures,
					vector<TrackInfo>& tracks);

private:
	double mAx, mAy, mAz;
};





//new api which can add control points
class DLL_EXPORT CPanoGlobalBA
{
public:
	CPanoGlobalBA();
	~CPanoGlobalBA();

	int Init(vector<CameraPara>& cameras,
			 vector<ImgFeature>& imageFeatures,
		     vector<TrackInfo>&  tracks);
	
	int AddFreePtBlock();

	int AddCtrlPtBlock(vector<stCtrlPt>& ctrlPts);

	int Run();

	int Output(vector<CameraPara>& cameras, vector<TrackInfo>&  tracks);


private:
	ceres::Problem mProblem;
	vector<CameraPara> mCameras;

	//projections
	int     mnProjNum;
	double* mpProjPt;

	//ground points
	int     mnGrdNum;
	double* mpGrdPt;

	//camera exterior parameters
	int     mnCamNum;
	double* mpCamPoseParas;

	//camera interior parameters
	double  mCamCaliParas[3]; //focal len, k1, k2

	vector<int> mVecCamIndex;    // camera index for each projection point
	vector<int> mVecTrackIndex;  // track point index for each projection point

};




#endif