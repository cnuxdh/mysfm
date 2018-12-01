
#ifndef MY_SFM_HPP
#define MY_SFM_HPP

#include"defines.hpp"
#include"badata.hpp"

//thread
#include<thread>

//#include<unordered_map>

//#include"pthread.h"


void *runthread(void* sfmsys);




typedef struct stSFMParameters
{
	string     projectPath;
	string     imagePath;
	string     outfilePath;
	int        maxHt;
	double     outResolution;
	CameraType camtype;
}SFMParameters;



//structure from motion by adding image one by one, written by xiedonghai, 2018.1.4
class DLL_EXPORT CSFMSystem
{
public:
	CSFMSystem();
	CSFMSystem(SFMParameters paras);


	int LoadParams(SFMParameters paras);


	~CSFMSystem();

	int SetProjectPath(string projectPath);
	int LoadImages(string path);

	//////////  sparse reconstruction ///////////
	int DetectFeatures(int maxHt);
	int ImageMatch(CameraType camtype);
	int GenerateTracks();
	int BA(CameraType camType);
	//int BAWithPos(CameraType camType);


	//////////  mosaic & fusion  //////////
	int AbsOrientation();
	int PointCloudSmooth(int numNeibor, double stddev_mult);
	int DEMInterpolation();
	int GenerateOrthoImages();
	int Fusion(double outResolution, string mosaicFile);
    int WeightFusion(double outResolution, string mosaicFile);

	//batch processing
	int Run();

	//stop thread
	int Stop();

	//io
	int SaveDEM(string filepath);
	int SaveCamsAndGrdpts(string filepath);
	int ReadCamsAndGrdpts(string filepath);


	//call back function
	void RegCallBack(sfmCall c) { mpCallback = c; }
	
	SFMParameters GetParameters() { return mSFMParas; }

private:
	
	sfmCall mpCallback;

	vector<string>       mImageFiles;   //the image file names
	string				 mProjectPath;  //the path saving temporary files


	//for match
	double mMedianCamDistance;

	//for BA
	vector<ImgFeature>   mImgFeatures;  //the point features in each image
	vector<PairMatchRes> mMatchRes;     //match index for each image pair
	vector<TrackInfo>    mTracks;       //tracks for BA
	vector<TrackInfo>    mGoodTracks;   //good tracks
	vector<CameraPara>   mCameras;      //camera parameters for BA
	
	//for generation of DOM
	double mdRawResolution;            //the resolution based on POS
	double mdMeanFlyHeight;            //the flying height from the ground
	double mdMeanSurfElevation;        //the mean height of surface 
	double mMinX, mMaxX, mMinY, mMaxY; //the mosaic area
	
	//for DEM
	vector<vector<double>> mdDEM;       //interpolated DEM
	double mdDEMResolution;

	//for mosaic
	vector<string>    mValidDOMFiles;
	vector<stGeoInfo> mValidDOMGeoData;

	//for progress bar
	double mdProgress;   //for progress bar 

	SFMParameters mSFMParas;

	//for thread
	//pthread_t mpThreadID;

	//to stop the thread, addedy by xiedonghai, 2018.7.6
	bool mbIsStop; 


	//handle
	//HANDLE mhCloseEvent;

	thread* mpThread;
	//typedef std::unordered_map<std::string, pthread_t> ThreadMap;
	//ThreadMap mTmap;
};


/*
//////////////////////////////// run in thread ///////////////////////////////
class DLL_EXPORT CSFMThread
{
public:
	CSFMThread();
	~CSFMThread();

	void    RunSFM(string imagePath, string projectPath, string outFilePath);
	double  GetProgressNumber() { return mdProgress; }
	string  GetProgressInfo() { return msProcessInfo; }

private:
	double mdProgress;
	string msProcessInfo;
};
*/


#endif