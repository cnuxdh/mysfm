
//#ifdef WIN32
//#include"windows.h"
//#include "atlimage.h"
//#endif

#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/stitching/stitcher.hpp>

#include<fstream>

#include"sfm.hpp"

#include"panorama.hpp"   //must included, because ceresba.hpp need it
#include"ceresba.hpp"    //must locate in the front, otherwise errors happen

////corelib
#include"commonfile.h"
#include"CommonFuncs.h"
#include"LatLong-UTMconversion.h"


////cvlibmdProgress
#include"sift.hpp"
#include"register.hpp"
#include"bundlerio.hpp"

//
#include"panoba.hpp"
//#include"smooth.h"
#include"absOri.hpp"
#include"OrthoImage.h"
#include"readJpegHeader.hpp"
#include"blend.hpp"
#include"mosaic.hpp"
#include"orbfeat.hpp"

//gdal
#include"geotiff.h"

#include"gdal_priv.h"
#include"ogr_spatialref.h"


////mba dll
#include"mbaExports.h"

#include<thread>

//#include<pthread.h>

using namespace cv;
using namespace std;


CSFMSystem::CSFMSystem()
{
	mpCallback = NULL;
	mbIsStop = false;
}

CSFMSystem::CSFMSystem(SFMParameters paras)
{
	mpCallback = NULL;
	mSFMParas = paras;
	mbIsStop = false;
}


CSFMSystem::~CSFMSystem()
{

}

int CSFMSystem::LoadParams(SFMParameters paras) {
	
	mSFMParas = paras;

	return 0;
}

int CSFMSystem::SaveCamsAndGrdpts(string filepath)
{
	char *sp = const_cast<char*>(filepath.c_str());

	FILE* fp = fopen(sp, "wb");

	//save camera data
	int camnum = mCameras.size();
	fwrite(&camnum,sizeof(int),1,fp);
	for (int i = 0; i < camnum; i++)
	{
		int camtype = mCameras[i].camtype;
		fwrite(&camtype, sizeof(int), 1, fp);
		fwrite(&(mCameras[i].focalLen), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].rows), sizeof(int), 1, fp);
		fwrite(&(mCameras[i].cols), sizeof(int), 1, fp);
		fwrite(mCameras[i].R, sizeof(double), 9, fp);
		fwrite(mCameras[i].T, sizeof(double), 3, fp);
		fwrite(&(mCameras[i].xs), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].ys), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].zs), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].gx), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].gy), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].gz), sizeof(double), 1, fp);
		fwrite(&(mCameras[i].bIsExplicit), sizeof(bool), 1, fp);
		fwrite(&(mCameras[i].bHavingPOS), sizeof(bool), 1, fp);
		fwrite(&(mCameras[i].bIsAddedIntoNet), sizeof(bool), 1, fp);
	}
	//save ground point data
	int nTracks = mGoodTracks.size();
	fwrite(&nTracks, sizeof(int), 1, fp);
	for (int i = 0; i < nTracks; i++)
	{
		double* p = mGoodTracks[i].GetGround().p;
		fwrite(p, sizeof(double), 3, fp);
	}

	fclose(fp);

	return 0;
}

int CSFMSystem::ReadCamsAndGrdpts(string filepath)
{
	char *sp = const_cast<char*>(filepath.c_str());

	FILE* fp = fopen(sp, "rb");

	//save camera data
	int camnum = 0;
	fread(&camnum, sizeof(int), 1, fp);
	mCameras.resize(camnum);
	for (int i = 0; i < camnum; i++)
	{
		int camtype = 0;
		fread(&camtype, sizeof(int), 1, fp);
		mCameras[i].camtype = CameraType(camtype);

		fread(&(mCameras[i].focalLen), sizeof(double), 1, fp);
		fread(&(mCameras[i].rows), sizeof(int), 1, fp);
		fread(&(mCameras[i].cols), sizeof(int), 1, fp);
		fread(mCameras[i].R, sizeof(double), 9, fp);
		fread(mCameras[i].T, sizeof(double), 3, fp);
		fread(&(mCameras[i].xs), sizeof(double), 1, fp);
		fread(&(mCameras[i].ys), sizeof(double), 1, fp);
		fread(&(mCameras[i].zs), sizeof(double), 1, fp);
		fread(&(mCameras[i].gx), sizeof(double), 1, fp);
		fread(&(mCameras[i].gy), sizeof(double), 1, fp);
		fread(&(mCameras[i].gz), sizeof(double), 1, fp);
		fread(&(mCameras[i].bIsExplicit), sizeof(bool), 1, fp);
		fread(&(mCameras[i].bHavingPOS), sizeof(bool), 1, fp);
		fread(&(mCameras[i].bIsAddedIntoNet), sizeof(bool), 1, fp);
	}

	//save ground point data
	int nTracks  = 0;
	fread(&nTracks, sizeof(int), 1, fp);
	mGoodTracks.resize(nTracks);
	//printf("track number: %d \n", )
	for (int i = 0; i < nTracks; i++)
	{
		double p[3]; 
		fread(p, sizeof(double), 3, fp);
		Point3DDouble gp;
		gp.p[0] = p[0];
		gp.p[1] = p[1];
		gp.p[2] = p[2];
		mGoodTracks[i].SetGround(gp);
	}

	fclose(fp);

	return 0;
}

int CSFMSystem::SetProjectPath(string projectPath)
{

	//judge if the project path exist
	mProjectPath = projectPath;

	return 0;
}

int CSFMSystem::LoadImages(string path)
{
	char *sp = const_cast<char*>(path.c_str());

	char** filenames = NULL;
	int n = 0, nfile = 0;
	GetDirFileName(filenames, sp, &n, &nfile, "jpg", 0);
	if (nfile==0)
	{
		GetDirFileName(filenames, sp, &n, &nfile, "JPG", 0);
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, sp, &n, &nfile, "JPG", 1);
	}
	else
	{
		filenames = f2c(nfile, 256);
		GetDirFileName(filenames, sp, &n, &nfile, "jpg", 1);
	}

	for (int i = 0; i < nfile; i++)
	{
		
		string st = filenames[i];
		mImageFiles.push_back(st);
		cout << st << endl;
	}
	
	//reading pos data
	ReadingPosFromJpegImages(filenames, nfile, mCameras);
	
	if (mbIsStop) {
		FreeArray_char(filenames, nfile, 256);
		return -1;
	}

	if (nfile < 2)
		return -1;

	//if camera has POS system
	if (mCameras[0].bHavingPOS){

		//convert from (lat,lon) to ground coordinates
		int nzone = int((mCameras[0].lon + 180) / 6) + 1;
		for (int i = 0; i < mCameras.size(); i++)
		{

			if (mbIsStop) {
				FreeArray_char(filenames, nfile, 256);
				return -1;
			}

			double lat = mCameras[i].lat;
			double lon = mCameras[i].lon;

			double gx, gy;
			LLtoUTM(23, lat, lon, gy, gx, nzone); //change the position of gx and gy, Xie Donghai, 2018.6.14

			mCameras[i].gx = gx;
			mCameras[i].gy = gy;
			mCameras[i].gz = mCameras[i].altitude;
		}

		//calculate the distance of neighbor image
		vector<double> camDistance;
		for (int i = 0; i < mCameras.size() - 1; i++)
		{

			if (mbIsStop) {
				FreeArray_char(filenames, nfile, 256);
				return -1;
			}

			double gx1 = mCameras[i].gx;
			double gy1 = mCameras[i].gy;
			double gx2 = mCameras[i+1].gx;
			double gy2 = mCameras[i+1].gy;
			double d = sqrt((gx1 - gx2)*(gx1 - gx2) + (gy1 - gy2)*(gy1 - gy2));
			camDistance.push_back(d);
		}

		//sort the distance
		std::sort(camDistance.begin(), camDistance.end());

		//calculate the median distance
		mMedianCamDistance = camDistance[int(camDistance.size()*0.5)];
	}

	
	FreeArray_char(filenames, nfile, 256);

	return 0;
}

int CSFMSystem::DetectFeatures(int maxHt)
{
	printf("[DetectFileFeaturePts] ... \n");

	//*mpsProcessInfo = "Detect Feature points....";

	CFeatureBase* pFeatDetect = new CSIFTFloat();
	//CFeatureBase* pFeatDetect = new CORBFeat(); //for orb
	int nFile = mImageFiles.size();

	double step = 100.0 / double(nFile);

	mpCallback(this, 0, "Detecting...");

	mImgFeatures.clear();
	//*mpdProgress = 0;
	for (int i = 0; i<nFile; i++)
	{
		//printf("%s \n", mImageFiles[i]);
		cout << mImageFiles[i] << endl;

		string title    = GetFileTitle(mImageFiles[i]);
		string featFile = mProjectPath + "\\" + title + ".feat";
		ifstream filein(featFile);

		ImgFeature feats;

		if (filein)
		{
			feats.ReadImgFeaturePts(featFile);
		}
		else{
			char *sp = const_cast<char*>(mImageFiles[i].c_str());
			printf("image: %s \n", sp);
			pFeatDetect->Detect(sp, feats, maxHt);	
			feats.SaveImgFeaturePts(featFile);
		}

		mImgFeatures.push_back(feats);

		mdProgress += step;

		mpCallback(this, mdProgress, "Detecting...");
	}

	delete pFeatDetect;

	return 0;
}

int CSFMSystem::ImageMatch(CameraType camtype)
{
	printf("[MatchImageFiles] ... \n");

	int nImageNum = mImgFeatures.size();
	
	CMatchBase* pMatch = NULL;

	if (camtype == PerspectiveCam)
	{
		pMatch = new CSiftMatch();  //new CKNNMatch();
		//pMatch = new CORBPerspectiveMatch();  //for orb;
	}
	else if (camtype == PanoramCam)
	{
		pMatch = new CPanoMatch();
	}

	mMatchRes.clear();
	mdProgress = 0;
	double step = 100.0 / double(nImageNum);
	for (int i = 0; i < nImageNum; i++)
	{
		printf("match for %d ... \n", i);

		double gx1 = mCameras[i].gx;
		double gy1 = mCameras[i].gy;

		for (int j = i + 1; j < nImageNum; j++)
		{

			if (mbIsStop) {
				return -1;
			}

			double gx2 = mCameras[j].gx;
			double gy2 = mCameras[j].gy;
			double d = sqrt((gx1 - gx2)*(gx1 - gx2) + (gy1 - gy2)*(gy1 - gy2));

			if (mCameras[i].bHavingPOS){
				if (d>mMedianCamDistance * 4){
					continue;
				}
			}

			//generate the match file name
			string title1 = GetFileTitle(mImageFiles[i]);
			string title2 = GetFileTitle(mImageFiles[j]);
			string mfile1 = mProjectPath + "\\" + title1 + "_" + title2 + ".match";
			string mfile2 = mProjectPath + "\\" + title2 + "_" + title1 + ".match";
			ifstream file1(mfile1);
			ifstream file2(mfile2);

			PairMatchRes mr;

			if (file1){
				mr.ReadMatchFile(mfile1);
			}
			else if (file2){
				mr.ReadMatchFile(mfile2);
			}
			else
			{
				pMatch->Match(mImgFeatures[i], mImgFeatures[j], mr);
				printf("%d-%d  %d %lf \n", i, j, mr.matchs.size(), mr.inlierRatio );
				
				//save match result
				mr.SaveMatchFile(mfile1);
			}

			mr.lId = i;
			mr.rId = j;

			printf("match: %d-%d  %d  F-inliers: %lf  H-inliers: %lf  \n", 
				i, j,  mr.matchs.size(), mr.inlierRatio, mr.inlierRatioHomography);

			mMatchRes.push_back(mr);

			mpCallback(this, mdProgress, "Matching...");
		}

		mdProgress += step;
	}

	delete pMatch;

	return 0;
}

int CSFMSystem::GenerateTracks()
{
	CGenerateTracksBase* pGenerateTrack = new CFastGenerateTrack();
	mTracks.clear();

	pGenerateTrack->GenerateTracks(mImgFeatures, mMatchRes, mTracks);

	delete pGenerateTrack;

	return 0;
}

int CSFMSystem::BA(CameraType camType)
{
	//*mpsProcessInfo = "Bundle Adjustment....";

	int numImage = mImgFeatures.size();

	if (numImage < 2)
		return -1;

	bool bHavePos = mCameras[0].bHavingPOS;

	//vector<CameraPara> cameras;
	//cameras.resize(numImage);
	double focalLen = (mImgFeatures[0].ht + mImgFeatures[0].wd) * 0.5;
	for (int i = 0; i<numImage; i++)
	{
		mCameras[i].focalLen = focalLen; //initialize the focal length 
		memset(mCameras[i].R, 0, sizeof(double) * 9);
		mCameras[i].R[0] = 1;
		mCameras[i].R[4] = 1;
		mCameras[i].R[8] = 1;
		mCameras[i].rows = mImgFeatures[i].ht;
		mCameras[i].cols = mImgFeatures[i].wd;
		mCameras[i].camtype = camType;
		mCameras[i].bIsAddedIntoNet = false;
	}

	string baResFile = mProjectPath + "\\" + "ba.out";
	char *sp = const_cast<char*>(baResFile.c_str());

	int res = 0;

	//if (IsFileExist(sp)){
	//	ReadCamsAndGrdpts(baResFile);
	//}
	//else
	{

		CBABase* pBA = NULL;
		if (camType == PerspectiveCam){
			pBA = new CCeresBA();
		}
		else if (camType == PanoramCam){
			pBA = new CPanoBA();
		}

		//pBA->SetProgress(mpdProgress);
		mpCallback(NULL, 0, "Bundle Adjustment...");

		pBA->SetCallBack(mpCallback);
		res = pBA->BundleAdjust(mCameras.size(), mCameras, mImgFeatures, 
            mMatchRes, mTracks, mbIsStop);

		delete pBA;
		

		if (res >= 0){
			GetGoodTracks(mTracks, mGoodTracks);
			//save the sparse clouds
			string plyfile = mProjectPath + "\\" + "bapt.ply";
			char* cplyfile = (char*)(plyfile.c_str());
			SaveTracksToPly(cplyfile, mGoodTracks, mCameras);
			//save the camera and tracks
			SaveCamsAndGrdpts(baResFile);	
		}
	}

	return res;
}


int CSFMSystem::PointCloudSmooth(int numNeibor, double stddev_mult)
{
	//get the good tracks
	vector<TrackInfo> goodTracks;
	GetGoodTracks(mTracks, goodTracks);
	
	int nPt = goodTracks.size();
	double* px = (double*)malloc(nPt*sizeof(double));
	double* py = (double*)malloc(nPt*sizeof(double));
	double* pz = (double*)malloc(nPt*sizeof(double));
	for (int i = 0; i<nPt; i++)
	{
		px[i] = goodTracks[i].grd.p[0];
		py[i] = goodTracks[i].grd.p[1];
		pz[i] = goodTracks[i].grd.p[2];
	}

	//point cloud filter
	//int nSmoothPt = PointCloudFilter(px, py, pz, nPt, numNeibor, stddev_mult);
	
	free(px);
	free(py);
	free(pz);

	return 0;
}


int CSFMSystem::AbsOrientation()
{
	if (mbIsStop) {
		return -1;
	}

	string aorResFile = mProjectPath + "\\" + "aor.out";
	char *sp = const_cast<char*>(aorResFile.c_str());

	/*if (IsFileExist(sp))
	{
		ReadCamsAndGrdpts(aorResFile);
	}
	else*/
	{
		//calculate the valid number of camera
		int nValidCamera = 0;
		for (int i = 0; i < mCameras.size(); i++)
		{
			if (mCameras[i].bIsAddedIntoNet)
				nValidCamera++;
		}
		
		/*if (nValidCamera < (mCameras.size()*0.5)){
			printf("imcomplete mosaic! \n");
			return -1;
		}*/

		if (nValidCamera < 3)
		{
			printf("the number of valid cameras is less than 3, exit! \n");
			return -1;
		}

		bool bIsHavePos = mCameras[0].bHavingPOS;

		//absolute orientation
		stAbsPOS absPosParams;
		absPosParams.scale = 1;
		memset(absPosParams.T, 0, 3 * sizeof(double));
		memset(absPosParams.R, 0, 9 * sizeof(double));
		absPosParams.R[0] = absPosParams.R[4] = absPosParams.R[8] = 1;
		if (bIsHavePos)
		{
			//AbsPosEstimation(absPosParams, camParas, tracks);
			double err = AbsOriOrthogonal(absPosParams, mCameras, mGoodTracks);
			printf("absolute orientation error: %lf \n", err);
		}
		else
		{
			printf("without absolute orientation ... \n");
			//find the ground plane 
			//(just like absolute pose estimation and transform the exterior parameters and tracks to new coordinate )
			FindGroundPlane(absPosParams, mCameras, mGoodTracks);
		}

		SaveCamsAndGrdpts(aorResFile);
	}

	string plyfile = mProjectPath + "\\" + "abspt.ply";
	char* cplyfile = (char*)(plyfile.c_str());
	SaveTracksToPly(cplyfile, mGoodTracks, mCameras);


	//calculate the mean height of ground surface
	double planeHeiMean = 0;
	int nValidCameras = 0;
	for (int i = 0; i<mCameras.size(); i++)
	{
		if (!mCameras[i].bIsAddedIntoNet)
			continue;
		planeHeiMean += mCameras[i].zs;
		nValidCameras++;
	}
	planeHeiMean /= double(nValidCameras);
	printf("plane mean height: %lf \n", planeHeiMean);


	//calculate the area and height of points after absolute pose estimation
	//double minx, maxx, miny, maxy;
	double meanHeight;
	CalculateGroudArea(mGoodTracks, mMinX, mMaxX, mMinY, mMaxY, meanHeight);
	mdMeanSurfElevation = meanHeight;
	mdMeanFlyHeight = fabs(meanHeight - planeHeiMean);
	printf("mean surface dem: %lf \n", mdMeanSurfElevation);
	printf("relative Height: %lf \n", mdMeanFlyHeight);


	//calculate the resolution after absolute pose estimation
	mdRawResolution = CalculateResolution(mdMeanSurfElevation, mCameras);
	mdRawResolution = max(0.05, mdRawResolution);
	printf("raw Resolution: %lf \n", mdRawResolution);
	

	return 0;
}


int CSFMSystem::DEMInterpolation()
{
	if (mbIsStop) {
		return -1;
	}

	string demfile = mProjectPath + "\\" + "dem.tif";
	char* sdemfile = const_cast<char*>(demfile.c_str());


	//get the good tracks
	//vector<TrackInfo> goodTracks;
	//GetGoodTracks(mTracks, goodTracks);

	
	int nPt = mGoodTracks.size();
	printf("good tracks: %d \n", mGoodTracks.size());


	double* px = (double*)malloc(nPt*sizeof(double));
	double* py = (double*)malloc(nPt*sizeof(double));
	double* pz = (double*)malloc(nPt*sizeof(double));
	for (int i = 0; i<nPt; i++)
	{
		px[i] = mGoodTracks[i].grd.p[0];
		py[i] = mGoodTracks[i].grd.p[1];
		pz[i] = mGoodTracks[i].grd.p[2];
	}

	//point cloud filter
	int nSmoothPt = nPt;
	//int numNeibor = 60;
	//double stddev_mult = 0.6;
	//nSmoothPt = PointCloudFilter(px, py, pz, nPt, numNeibor, stddev_mult);

	//dem interpolation
	mdDEMResolution = mdRawResolution * 10;  //define by xiedonghai
	int demHt = (mMaxY - mMinY) / mdDEMResolution;
	int demWd = (mMaxX - mMinX) / mdDEMResolution;
	float* iz = new float[demHt*demWd];
	MBAInterpolation(px, py, pz, nSmoothPt, 10, demHt, demWd, iz);

	//median filter

	//save it
	mdDEM.resize(demHt);
	for (int i = 0; i < demHt; i++)
		mdDEM[i].resize(demWd);
	for (int j = 0; j < demHt; j++)
		for (int i = 0; i < demWd; i++)
			mdDEM[j][i] = iz[j*demWd + i];
	
	printf("saving dem into file ... \n");
	if (1)
	{
		stGeoInfo geoinfo;	
		memset(&geoinfo, 0, sizeof(stGeoInfo));
		geoinfo.left = mMinX;
		geoinfo.top  = mMaxY;
		geoinfo.dx = mdDEMResolution;
		geoinfo.dy = mdDEMResolution;
		
		GdalWriteFloat(sdemfile, iz, demHt, demWd, geoinfo);
	}
	
	delete iz;
	free(px);
	free(py);
	free(pz);

	return 0;
}

int CSFMSystem::GenerateOrthoImages()
{	
	//*mpsProcessInfo = "Generate OrthoImage....";

	int zoneNumber = GetZoneNumber(mCameras);

	mValidDOMFiles.clear();
	mValidDOMGeoData.clear();

	mdProgress = 0;
	double step = 100.0 / double(mImageFiles.size());
	for (int i = 0; i < mImageFiles.size(); i++)
	{
		if (mbIsStop) {
			return -1;
		}

		mdProgress += step;

		if (!mCameras[i].bIsAddedIntoNet) 
			continue;

		string title = GetFileTitle(mImageFiles[i]);
		string domfile = mProjectPath + "\\" + title + ".jpeg";
		string domGeofile = mProjectPath + "\\" + title + ".geo";

		//int npos = domfile.find(".");
		//domfile.replace("jpeg", npos + 1);
		char *sp = const_cast<char*>(domfile.c_str());

		//if (IsFileExist(sp))
		//{
		//	
		//	mValidDOMFiles.push_back(domfile);

		//	//read the geoinfo into memory
		//	stGeoInfo geoinfo;
		//	char *sp = const_cast<char*>(domGeofile.c_str());
		//	FILE* fp = fopen(sp, "r");
		//	if (fp != NULL){
		//		fscanf(fp, "%d %d %d  %lf %lf %lf %lf ", &(geoinfo.zoneNumber),
		//			&(geoinfo.ht), &(geoinfo.wd),
		//			&(geoinfo.dx), &(geoinfo.dy),
		//			&(geoinfo.left), &(geoinfo.top));
		//		fclose(fp);
		//	}

		//	mValidDOMGeoData.push_back(geoinfo);

		//	continue;
		//}

		stGeoInfo geoinfo;
		vector<vector<unsigned char>> rgb;

		geoinfo.zoneNumber = zoneNumber;

		int oht, owd;
		double outResolution = max(0.1, mdRawResolution * 2);
		GenerateRGBDOM(mImageFiles[i], mdRawResolution, outResolution,
			mdDEMResolution, mMinX, mMinY, mMaxX, mMaxY, mCameras[i], 
			mdMeanSurfElevation,
			mdDEM, geoinfo, rgb, oht, owd);

		//write into the disk
		//string domfile = mImageFiles[i]+".jpeg";

		int domHt = oht;
		int domWd = owd; 
		printf("%d %d \n", domHt, domWd);
		printf("generating r,g,b.... \n");
		unsigned char* r = (unsigned char*)malloc(domHt*domWd);
		unsigned char* g = (unsigned char*)malloc(domHt*domWd);
		unsigned char* b = (unsigned char*)malloc(domHt*domWd);
		if (r == NULL || g == NULL || b == NULL)
		{
			printf("malloc failed ! \n");
		}

		int index = 0;
		for (int j = 0; j < domHt; j++)
			for (int i = 0; i < domWd; i++)
			{
				r[index] = rgb[0][index];
				g[index] = rgb[1][index];
				b[index] = rgb[2][index];
				index++;
			}
		
		printf("ortho file: %s \n", sp);
		GdalWriteImageByteColor(sp, r, g, b, domHt, domWd);
		free(r);
		free(g);
		free(b);

		//write the geoinfo into file
		sp = const_cast<char*>(domGeofile.c_str());
		FILE* fp = fopen(sp, "w");
		fprintf(fp, "%d %d %d  %lf %lf %lf %lf ", geoinfo.zoneNumber,
			geoinfo.ht, geoinfo.wd, 
			geoinfo.dx, geoinfo.dy, 
			geoinfo.left, geoinfo.top);
		fclose(fp);

		//save the dom into disk
		mValidDOMFiles.push_back(domfile);
		//save the geoinfo into memory
		mValidDOMGeoData.push_back(geoinfo);

		mpCallback(NULL, mdProgress, "Generating Single Orthoimage....");
	}	

	return 0;
}

int CSFMSystem::WeightFusion(double outResolution, string mosaicFile)
{
    //mosaic and fusion	
    double minx, maxx, miny, maxy;
    double maxrx, maxry;
    double rx, ry;
    minx = 5000000000;	maxx = -5000000000;
    miny = 5000000000;	maxy = -5000000000;
    maxrx = 0;
    maxry = 0;
    //int zoneNumber = mValidDOMGeoData[0].zoneNumber;
    for (int i = 0; i<mValidDOMGeoData.size(); i++)
    {
        //if (zoneNumber != geoArray[i].zoneNumber)
        //	continue;		
        //resolution
        rx = mValidDOMGeoData[i].dx;
        ry = fabs(mValidDOMGeoData[i].dy);
        maxrx = max(rx, maxrx);
        maxry = max(ry, maxry);
        //position
        minx = min(minx, mValidDOMGeoData[i].left);
        maxx = max(maxx, mValidDOMGeoData[i].left + mValidDOMGeoData[i].wd*rx);
        miny = min(miny, mValidDOMGeoData[i].top - mValidDOMGeoData[i].ht*ry);
        maxy = max(maxy, mValidDOMGeoData[i].top);
    }
    printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

    //double outResolution = 0.1;
    double resolution = outResolution;
    int oht = (maxy - miny) / resolution;
    int owd = (maxx - minx) / resolution;
    printf("%d %d \n", oht, owd);
    Mat oimage(oht, owd, CV_16SC3, Scalar(0, 0, 0));

    std::vector<cv::Mat> out_pyr_laplace;
    int nlevel = 4;
    cv::detail::createLaplacePyr(oimage, nlevel, out_pyr_laplace);
    printf("size of pyramid: %d \n", out_pyr_laplace.size());

    std::vector<cv::Mat> out_pyr_weights(nlevel + 1);
    out_pyr_weights[0].create(oht, owd, CV_32FC1);
    out_pyr_weights[0] = Scalar(0);
    for (int i = 0; i < nlevel; ++i)
        cv::pyrDown(out_pyr_weights[i], out_pyr_weights[i + 1]);

    for (int fi = 0; fi < mValidDOMFiles.size(); fi++)
    {
        //convert the image type
        Mat image = imread(mValidDOMFiles[fi]);
        Mat srcimage;
        image.convertTo(srcimage, CV_16SC3);
        
        //create laplacian images
        std::vector<cv::Mat> pyr_laplace;
        cv::detail::createLaplacePyr(srcimage, nlevel, pyr_laplace);
        printf("size of pyramid: %d \n", pyr_laplace.size());

        //create weight image
        Mat weightImage;
        int w = image.cols;
        int h = image.rows;
        weightImage.create(h, w, CV_32FC1);
        float *p = (float*)weightImage.data;
        float x_center = w / 2;
        float y_center = h / 2;
        float dis_max = sqrt(x_center*x_center + y_center * y_center);
        int weightType = 0; // svar.GetInt("Map2D.WeightType", 0);
        
        for (int i = 0; i < h; i++)
        {
            uchar* pbuffer = image.ptr<uchar>(i);
            for (int j = 0; j < w; j++)
            {
                int r = pbuffer[j * 3];
                int g = pbuffer[j * 3 + 1];
                int b = pbuffer[j * 3 + 2];
                if ((r + g + b) == 0)
                {
                    *p = 0;
                }
                else
                {
                    float dis = (i - y_center)*(i - y_center) + (j - x_center)*(j - x_center);
                    dis = 1 - sqrt(dis) / dis_max;
                    if (0 == weightType)
                        *p = dis;
                    else *p = dis * dis;
                    if (*p <= 1e-5) *p = 1e-5;
                }
                p++;
            }
        }

        //create weight pyramid
        std::vector<cv::Mat> pyr_weights(nlevel + 1);
        pyr_weights[0] = weightImage;
        for (int i = 0; i < nlevel; ++i)
            cv::pyrDown(pyr_weights[i], pyr_weights[i + 1]);

        //fill the pyramid
        int l = (mValidDOMGeoData[fi].left - minx) / outResolution;
        int t = (maxy - mValidDOMGeoData[fi].top) / outResolution;
        for (int i = 0; i <= nlevel; i++)
        {
            int pht = pyr_laplace[i].rows;
            int pwd = pyr_laplace[i].cols;
            printf("%d %d \n", pht, pwd);
            for (int m = 0; m < pht; m++)
            {
                //weight data
                int mt = m + t;
                if (mt >= oht) mt = oht - 1;
                short* outp = out_pyr_laplace[i].ptr<short>(mt);
                short* inp = pyr_laplace[i].ptr<short>(m);

                //image data
                float* outw = out_pyr_weights[i].ptr<float>(mt);
                float* inw = pyr_weights[i].ptr<float>(m);

                for (int n = 0; n < pwd; n++)
                {
                    //ignore the background
                    //int sg = inp[n * 3] + inp[n * 3 + 1] + inp[n * 3 + 2];
                    //if (sg == 0)  continue;
                    if (inw[n] == 0)
                        continue;

                    //update the image and weight
                    int nl = n + l;
                    if (nl >= owd) nl = owd - 1;
                    if (inw[n] > outw[nl])
                    {
                        outw[nl] = inw[n];
                        outp[(nl) * 3] = inp[n * 3];
                        outp[(nl) * 3 + 1] = inp[n * 3 + 1];
                        outp[(nl) * 3 + 2] = inp[n * 3 + 2];
                    }
                }
            }
            //zoom
            l = 0.5*l;
            t = 0.5*t;
        }
    }

    cv::detail::restoreImageFromLaplacePyr(out_pyr_laplace);
    Mat result = out_pyr_laplace[0].clone();
    if (result.type() == CV_16SC3)
        result.convertTo(result, CV_8UC3);
    result.setTo(cv::Scalar::all(0), out_pyr_weights[0] == 0);

    
    //save as tiff file    
    int zoneNumber = mValidDOMGeoData[0].zoneNumber; //geoinfo.zoneNumber;
    OGRSpatialReference oSRS;
    oSRS.SetUTM(zoneNumber);
    oSRS.SetWellKnownGeogCS("WGS84");
    char    *pszWKT = NULL;
    oSRS.exportToWkt(&pszWKT);
    stGeoInfo geoinfo;
    geoinfo.left = minx;
    geoinfo.top =  maxy;
    geoinfo.dx = resolution;
    geoinfo.dy = -resolution;
    geoinfo.projectRef = pszWKT;
    
    int rht = result.rows;
    int rwd = result.cols;
    unsigned char* pr = new unsigned char[rht*rwd];
    unsigned char* pg = new unsigned char[rht*rwd];
    unsigned char* pb = new unsigned char[rht*rwd];
    int index = 0;
    for (int m = 0; m < rht; m++)
    {
        uchar *p = result.ptr<uchar>(m);
        for (int n = 0; n < rwd; n++)
        {
            pr[index] = p[n * 3];
            pg[index] = p[n * 3 + 1];
            pb[index] = p[n * 3 + 2];
            index++;
        }
    }





    delete[] pr;
    delete[] pg;
    delete[] pb;

    return 0;
}

int CSFMSystem::Fusion(double outResolution, string mosaicFile)
{	
	if (mbIsStop) {
		return -1;
	}

	//*mpsProcessInfo = "Mosaic and Fusion....";
	mpCallback(NULL, 0, "Mosaic and Fusion...");

	int nFile = mValidDOMFiles.size();

	//mosaic and fusion	
	double minx, maxx, miny, maxy;
	double maxrx, maxry;
	double rx, ry;
	minx =  5000000000;	maxx = -5000000000;
	miny =  5000000000;	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = mValidDOMGeoData[0].zoneNumber;
	for (int i = 0; i<mValidDOMFiles.size(); i++)
	{
		//if (zoneNumber != geoArray[i].zoneNumber)
		//	continue;		
		//resolution
		rx = mValidDOMGeoData[i].dx;
		ry = fabs(mValidDOMGeoData[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);
		//position
		minx = min(minx, mValidDOMGeoData[i].left);
		maxx = max(maxx, mValidDOMGeoData[i].left + mValidDOMGeoData[i].wd*rx);
		miny = min(miny, mValidDOMGeoData[i].top - mValidDOMGeoData[i].ht*ry);
		maxy = max(maxy, mValidDOMGeoData[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);
	
	//get the proper resolution
	double resolution = outResolution;
	int oht = (maxy - miny) / resolution;
	int owd = (maxx - minx) / resolution;
	double memRatio = CalculateRatioFor64bitOS(oht, owd, 1);
	oht *= memRatio;
	owd *= memRatio;
	outResolution = resolution / memRatio;
	printf("%lf %d %d \n", memRatio, oht, owd);
	printf(" restricted resolution: %lf \n", resolution);
	//WriteProgressValueToFile(5.0);

	//////////////  2. generate the original mask image and save them //////////////
	//set the resolution for the mask
	int iht = mValidDOMGeoData[0].ht;
	int iwd = mValidDOMGeoData[0].wd;
	int nMaskSize = max(iht, iwd);
	double maskResolution = fabs(mValidDOMGeoData[0].dx);
	if (nMaskSize>MAX_MASK_SIZE)
	{
		double maskRatio = (double)(MAX_MASK_SIZE) / (double)(nMaskSize);
		maskResolution = maskResolution / maskRatio;
	}
	if (maskResolution<outResolution)
		maskResolution = outResolution;

	vector<MyRect> vecRectMosaic;
	vecRectMosaic.resize(nFile);
	vector<MyRect> vecRectMask;
	vecRectMask.resize(nFile);
	
	unsigned char** pMaskArray = (unsigned char**)malloc(nFile*sizeof(unsigned char*));
	vector<int> masksHt;
	vector<int> masksWd;
	masksHt.resize(nFile);
	masksWd.resize(nFile);
	
	printf("generating original mask... \n");
	for (int i = 0; i<mValidDOMFiles.size(); i++)
	{
		printf("%d \n", i);
		double ratio = fabs(mValidDOMGeoData[i].dx) / maskResolution;		
		stGeoInfo geoInfo;
		char *sp = const_cast<char*>(mValidDOMFiles[i].c_str());
		GetGeoInformation(sp, geoInfo);
		int ht = geoInfo.ht;
		int wd = geoInfo.wd;
		masksHt[i] = ht*ratio;
		masksWd[i] = wd*ratio;
		printf("mask: %d  ht-%d wd-%d \n", i, masksHt[i], masksWd[i]);
		ReadGeoFileByte(sp, 1, ratio, &pMaskArray[i], masksHt[i], masksWd[i]);

		//convert to mask value
		for (int k = 0; k<masksHt[i] * masksWd[i]; k++)
		{
			if (pMaskArray[i][k]>0)
				pMaskArray[i][k] = 10;
		}

		//calculate the translation of each image in mask space
		double mosaicRatio = fabs(mValidDOMGeoData[i].dx) / outResolution;
		int x = (mValidDOMGeoData[i].left - minx) / outResolution;
		int y = (maxy - mValidDOMGeoData[i].top) / outResolution;
		vecRectMosaic[i].top = y;
		vecRectMosaic[i].left = x;
		vecRectMosaic[i].bottom = y + mValidDOMGeoData[i].ht*mosaicRatio + 0.5;
		vecRectMosaic[i].right = x + mValidDOMGeoData[i].wd*mosaicRatio + 0.5;

		double maskRatio = fabs(mValidDOMGeoData[i].dx) / maskResolution;
		x = (mValidDOMGeoData[i].left - minx) / maskResolution;
		y = (maxy - mValidDOMGeoData[i].top) / maskResolution;
		vecRectMask[i].top = y;
		vecRectMask[i].left = x;
		vecRectMask[i].bottom = y + mValidDOMGeoData[i].ht*maskRatio + 0.5;
		vecRectMask[i].right = x + mValidDOMGeoData[i].wd*maskRatio + 0.5;
	}

	//////////////////////   3. gain compensation ////////////////////////////////////
	double* pg = (double*)malloc(nFile*sizeof(double));
	for (int i = 0; i<nFile; i++)
		pg[i] = 1;

	//calculate the gain only for RGB true colr images
	//char* postfix;
	//GetPostfix(filenames[0], &postfix);
	//if (strcmp(postfix, "jpg") == 0 || strcmp(postfix, "JPG") == 0)
	//{
	//	CalculateGain(maskNames, filenames, nFile, vecRectMask, pg);
	//}
	//free(postfix);


	//////////////////////   4. finding seams ////////////////////////////////////
	int res = SeamlinesFindingVoronoi(pMaskArray, nFile, masksHt, masksWd, vecRectMask, mbIsStop);

	if (res < 0)
		return -1;
	

#ifdef _DEBUG
	for (int k = 0; k < nFile; k++)
	{
		char mfile[256];
		sprintf(mfile, "c:\\temp\\m_%d.jpg", k);
		SaveToJpg(pMaskArray[k], masksHt[k], masksWd[k], mfile);
	}
#endif

	int nPyramidLevel = CalculateSeqsPyramidLevel(mValidDOMGeoData, outResolution);
	
	char* outFile = const_cast<char*>(mosaicFile.c_str());
	char** filenames = f2c(mValidDOMFiles.size(), 512); // (char**)malloc(sizeof(unsigned char*));

	printf("dom file number: %d \n", mValidDOMFiles.size());
	for (int i = 0; i < mValidDOMFiles.size(); i++)
	{
		//filenames[i] = const_cast<char*>(mValidDOMFiles[i].c_str());
#ifdef _DEBUG
		char* ps = (char*)(mValidDOMFiles[i].c_str());
		strcpy(filenames[i], ps);
#else
		memset(filenames[i], '\0', 512);
		int strlen = mValidDOMFiles[i].length();
		mValidDOMFiles[i].copy(filenames[i], strlen);
#endif
	}


	WeightedLoGBlendGeneral(filenames, pMaskArray, nFile, masksHt, masksWd,
		mValidDOMGeoData, vecRectMosaic,
		resolution, maskResolution, minx, maxx, miny, maxy, nPyramidLevel, pg, outFile, 
		mdProgress, mpCallback, mbIsStop);

	//release memory
	for (int i = 0; i < nFile; i++)
	{
		free(pMaskArray[i]);
	}
	free(pMaskArray);
	//free(filenames);
	FreeArray_char(filenames, nFile, 512);
	free(pg);
	
	mdProgress = 0;

	mpCallback(NULL, 0, "Finished!");

	return 0;
}
int CSFMSystem::SaveDEM(string filepath)
{

	return 0;
}

//void CSFMSystem::SetProgress(double* pdProgress, string* psInfo)
//{
//	//mpdProgress = pdProgress;
//	//mpsProcessInfo = psInfo;
//}

int CSFMSystem::Run()
{

	//pthread_create(&mpThreadID, NULL, runthread, this);


	//start the thread
	//thread* pt = new thread(runthread, this);
	
	
	//string tname = "mosaic";
	mpThread = new thread(runthread, this);
	//mTmap[tname] = (mpThread->native_handle());
	mpThread->detach();
	//std::cout << "Thread " << tname << " created:" << std::endl;
	

	return 0;
}

int CSFMSystem::Stop(){

	mbIsStop = true;

	//pthread_cancel(mpThreadID);
	//pthread_kill(mpThreadID, 0);

	/*
	string tname = "mosaic";
	//mpThread->detach();
	ThreadMap::const_iterator it = mTmap.find(tname);
	if (it != mTmap.end()) {
		//delete it->second; // thread not killed
		//it->second->std::thread::~thread(); // thread not killed
		pthread_cancel(it->second);
		mTmap.erase(tname);
		std::cout << "Thread " << tname << " killed:" << std::endl;
	}
	*/
	
	return 0;
}


void* runthread(void* sfmsys)
{
	CSFMSystem* pSFM = (CSFMSystem*)(sfmsys);

	string imagepath   = pSFM->GetParameters().imagePath;
	string projectpath = pSFM->GetParameters().projectPath;
	string mosaicfile  = pSFM->GetParameters().outfilePath;
	CameraType camType = PerspectiveCam;
	double outReslution = pSFM->GetParameters().outResolution;
	int maxHt = pSFM->GetParameters().maxHt;

	pSFM->SetProjectPath(projectpath);
	pSFM->LoadImages(imagepath);
	pSFM->DetectFeatures(maxHt);
	pSFM->ImageMatch(camType);
	pSFM->GenerateTracks();

	if (pSFM->BA(camType) < 0)
	{
		printf("BA failed! \n");
		return NULL;
	}

	//generate orthoimage and mosaic
	if (pSFM->AbsOrientation() >= 0) {
		pSFM->DEMInterpolation();
		pSFM->GenerateOrthoImages();
		pSFM->Fusion(outReslution, mosaicfile);
	}
	else {
		printf("absolute orientation failed! \n");
	}	

	return NULL;
}


/*
///////////////////////////// for thread /////////////////////////////////

void runsfm(string imagePath, string projectPath, string outFilePath,
	double* pdProgress, string* psProcessInfo)
{

	CSFMSystem sfm;
	
	//sfm.SetProgress(pdProgress, psProcessInfo);
	sfm.SetProjectPath(projectPath);
	sfm.LoadImages(imagePath);
	sfm.DetectFeatures(640);

	CameraType camType = PerspectiveCam;
	sfm.ImageMatch(camType);
	sfm.GenerateTracks();

	if (sfm.BA(camType) < 0)
	{
		printf("BA failed! \n");
		//return -1;;
	}

	//generate orthoimage and mosaic
	if (sfm.AbsOrientation() >= 0) {
		sfm.DEMInterpolation();
		sfm.GenerateOrthoImages();
		sfm.Fusion(0.1, outFilePath);
	}
	else {
		printf("absolute orientation failed! \n");
	}
}


CSFMThread::CSFMThread()
{

}
CSFMThread::~CSFMThread()
{

}
void CSFMThread::RunSFM(string imagePath, string projectPath, string outFilePath)
{
	thread t;
	t = thread(runsfm, imagePath, projectPath, outFilePath, &mdProgress, &msProcessInfo);
	t.detach();
	//t.join();
}
*/
