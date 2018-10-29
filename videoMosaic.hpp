
#ifndef  VIDEO_MOSAIC_HPP
#define  VIDEO_MOSAIC_HPP

#include"commondata.h"

//cvlib
#include"badata.hpp"
#include"videoMosaic.hpp"
#include"defines.hpp"


//parameter class for fast video mosaic, written by xdh, 2018.1.13
class DLL_EXPORT CPlaneMosaicPara
{
public:
	CPlaneMosaicPara();
	~CPlaneMosaicPara();
	
	int SetSimilarityParas(double a, double b, double c, double d);
	int GetSimilarityParas(double& a, double& b, double& c, double& d);

	double ma, mb, mc, md; //similarity paramters  x'=ax-by+c, y'=bx+ay+d
};

class CVideoKeyFrame
{
public:
	CVideoKeyFrame(){};
	~CVideoKeyFrame(){};

	void Clear() { mFeature.Clear(); }

	ImgFeature& GetFeat(){ return mFeature; }
	CPlaneMosaicPara& GetCamPara(){ return mCamera; }

//private:
	ImgFeature mFeature;
	CPlaneMosaicPara mCamera;
};


//system for video mosaic
class DLL_EXPORT CVideoMosaicSystem
{
public:
	CVideoMosaicSystem();
	~CVideoMosaicSystem();
	
	//int Detect(Mat& frame);
	//int MatchWithClosestFrame();

	bool Initialize();
	bool IsKeyFrame();

	int  MatchingWithClosetFrame();
	int  SimTransWithLastKeyFrame(
		double& sa, double& sb, double& sc, double& sd, double& inlierRatio);

	//int  AddingKeyFrame();
	int  CalculateNewCamPara();

	int  WarpImage(Mat& mosaic);

	int  LocalBA();

	int  Run(Mat& frameImage);
	
private:
	
	vector<TrackInfo>  mTracks;
	vector<CVideoKeyFrame> mKeyFrames;

	//for current video frame feature and match
	CVideoKeyFrame mCurrentFrame;
	PairMatchRes   mNewFrameMatch;

	//current frame image
	Mat mCurrentFrameImage;

	bool mbIsHavingKeyFrame;

	bool mbIsInitialized;
};


class CVideoFrameBA
{
public:
	CVideoFrameBA();
	~CVideoFrameBA();

	int Run();
	
private:
	
};



DLL_EXPORT int SimilarityRANSAC(MyPointF* lImagePts, MyPointF* rImagePts, int nPt,
	double& a, double& b, double& c, double& d, vector<int>& vInliers);

DLL_EXPORT int VideoMosaicInitialize(vector<TrackInfo>&  tracks, vector<CVideoKeyFrame>& keyFrames);

DLL_EXPORT int VisibleTracksInitialize(vector<CVideoKeyFrame>& keyFrames,
	vector<TrackInfo>&  tracks, vector<int> visibleIndex);

//get the visible tracks before inserting the new image
DLL_EXPORT int GetCurrentFrameVisibleTracks(vector<TrackInfo>& tracks,
	vector<CVideoKeyFrame>& keyFrames,
	vector<PairMatchRes> pm,
	int currentFrameIndex,
	vector<int>& visibleTrackIndex,
	vector<int>& currentFramePtIndex);

//get the visible tracks after inserting the new image
DLL_EXPORT int GetCurrentVisibleTracksNew(vector<TrackInfo>& tracks,
	ImgFeature& currentFrame, vector<int>& visibleTrackIndex);



#endif
