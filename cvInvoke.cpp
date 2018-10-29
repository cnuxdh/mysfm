
#include"stdio.h"

#include"cvInvoke.hpp"
#include"sift.hpp"
#include"dataBase.hpp"
#include"register.hpp"
#include"relativepose.hpp"
#include "CalcAngle.h"
#include "triangulate.hpp"
#include "ORBextractor.h"
#include "orbfeat.hpp"


//corelib
#include"commonfile.h"
#include"ImageFunc.h"
#include"commondata.h"


#include "baselib.h"

//rslib
#include"Aerosol.h"
#include"lut.h"
#include"modis.h"

//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgproc/imgproc.hpp"

using namespace cv;
//using namespace ORB_SLAM2;


// function for feature detection of multiple images, for input files
// Input
//	filenames : the image files
//  nFile :
// Output 
//	outpath: image feature files saved into the output path
int DetectFileFeaturePts(char** filenames, int nFile, char* outpath)
{
	//retrive the file title
	
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	CPointFeatureBase* pFeatureData = new CSiftFeatureDataBinary();
	
	for(int i=0; i<nFile; i++)
	{
		printf("image: %s \n", filenames[i]);
		
		ImgFeature feats;
		pFeatDetect->Detect(filenames[i], feats);
		
		//save the feat into the file
		//printf("save the feature points into the file.... \n");
		char* title;
		GetTitleName(filenames[i], &title);
		printf("title: %s \n", title);
		
		char featPath[512];
		sprintf(featPath, "%s/%s.dat", outpath, title);
		printf("%s \n", featPath);
		
		pFeatureData->Write(featPath, feats);		
	}	
	delete pFeatDetect;
	
	return 0;
}


int DLL_EXPORT DetectFileFeaturePts(char* filenames, ImgFeature& imgFeatures, int maxHt)
{
	CFeatureBase* pFeatDetect = new CSIFTFloat();

	pFeatDetect->Detect(filenames, imgFeatures, maxHt);

	delete pFeatDetect;

	return 0;
}



//all features are saved in memory, and progress is also added, 2017.12.31
int DetectFileFeaturePts(char** filenames, int nFile, vector<ImgFeature>& imgFeatures,
	int maxHt, double& dProgress)
{
	printf("[DetectFileFeaturePts] ... \n");

	CFeatureBase* pFeatDetect = new CSIFTFloat();
	//CFeatureBase* pFeatDetect = new CORBFeat(); //for orb
	
	double step = 100.0 / double(nFile);

	for(int i=0; i<nFile; i++)
	{
		printf("image: %s \n", filenames[i]);
		
		ImgFeature feats;
		pFeatDetect->Detect(filenames[i], feats, maxHt);
		imgFeatures.push_back(feats);
	
		dProgress += step;
	}
	
	delete pFeatDetect;
	return 0;
}


//matching between image files
int MatchImageFiles(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& matchRes, 
	CameraType camtype, int matchSteps)
{
	printf("[MatchImageFiles] ... \n");

	int nImageNum = imgFeatures.size();
		

	CMatchBase* pMatch = NULL;
	
	if(camtype==PerspectiveCam)
	{
		pMatch = new CSiftMatch();  //new CKNNMatch();
		//pMatch = new CORBPerspectiveMatch();  //for orb;
	}
	else if(camtype==PanoramCam)
	{
		pMatch = new CPanoMatch();
	}

	for (int i = 0; i < nImageNum; i++)
	{
		int end = min(nImageNum, i + 1 + matchSteps);



		for (int j = i + 1; j < end; j++)
		{

			PairMatchRes mr;
			mr.lId = i;
			mr.rId = j;

			pMatch->Match(imgFeatures[i], imgFeatures[j], mr);

			//printf("%d-%d  %d %lf \n", i, j, mr.matchs.size(), mr.inlierRatio );
			printf("%d-%d  %d  \n", i, j, mr.matchs.size());

			matchRes.push_back(mr);
		}
	}
	
	delete pMatch;

	return 0;
}


DLL_EXPORT int dll_EstimatePose5Point_Pano( vector<Point3DDouble>& pl, 
	vector<Point3DDouble>& pr,
	double radius,
	int num_trials, double threshold, 
	double *R, double *t, vector<double>& residual)
{
	return EstimatePose5Point_Pano(pl, pr, radius, num_trials, threshold, R, t, residual);
}


DLL_EXPORT int dll_EstimatePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts,
	CameraPara& cam1, CameraPara& cam2, CameraType camtype )
{
	if(camtype == PerspectiveCam)
	{
		CRelativePoseBase* pRP = new CEstimatePose5Point(); 
		pRP->EstimatePose(lPts, rPts, cam1, cam2);
		delete pRP;
	}

	if(camtype == PanoramCam)
	{
		CRelativePoseBase* pRP = new CEstimatePose5PointPano(); 
		pRP->EstimatePose(lPts, rPts, cam1, cam2);
		delete pRP;
	}

	return 0;
}

DLL_EXPORT Point3DDouble dll_TriangulatePt(Point2DDouble p, Point2DDouble q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error)
{
	return  TriangulatePt(p, q, R0, t0, R1, t1, error);
}



DLL_EXPORT void dll_Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, Point3DDouble& gps,
	bool explicit_camera_centers,double& ferror)
{

	CTriangulateBase* pTriangulate = new CTriangulateCV();
	pTriangulate->Triangulate(pts, cams, gps, explicit_camera_centers, ferror);
	delete pTriangulate;
}


void dll_GenerateRMatrix(double omiga, double phi, double kapa, double* R)
{
	GenerateRMatrix(omiga, phi, kapa, R);
}

void dll_GrdToImg(double gx, double gy, double gz, double* ix, double* iy, 
	double* R, double* Ts, double f, double x0, double y0,
	int ht, int wd)
{
	GrdToImg(gx, gy, gz, ix, iy, R, Ts, f, x0, y0, ht, wd);
}

int dll_DLT(vector<Point3DDouble>& grds, vector<Point2DDouble>& projs,
			CameraPara& cam, CameraType camType)
{
	int res = 0;

	if(camType == PerspectiveCam)
	{
		//for normal perspective camera
		CPoseEstimationBase* pDLT = new CDLTPose();
		res = pDLT->EstimatePose(grds, projs, cam);
		delete pDLT;
	}

	//for panorama camera
	if(camType == PanoramCam)
	{
		//for normal perspective camera
		CPoseEstimationBase* pDLT = new CPanoDLTPose();
		res = pDLT->EstimatePose(grds, projs, cam);
		delete pDLT;
	}

	return res;
}


DLL_EXPORT double dll_CalculatePanoEpipolarError(double* em, Point3DDouble lp, Point3DDouble rp, double radius)
{
	double error = 0;

	v3_t vl,vr;
	vl.p[0] = lp.p[0];	vl.p[1] = lp.p[1];	vl.p[2] = lp.p[2];
	vr.p[0] = rp.p[0];	vr.p[1] = rp.p[1];	vr.p[2] = rp.p[2];
	
	error = dll_fmatrix_compute_residual_pano(em, vl, vr, radius);

	return error;
}

int dll_GenerateRainbowMapping(vector<int>& r, vector<int>& g, vector<int>& b)
{
	GenerateRainBowMapping(r, g, b);
	
	return 0;
}




int dll_GenerateORBFeature(string filepath, vector<KeyPoint>& Keys, Mat& Descriptors)
{	
	// Read image from file
	Mat im;
	im = imread(filepath, CV_LOAD_IMAGE_UNCHANGED);
	
	if (im.channels() == 3)
	{
		cvtColor(im, im, CV_RGB2GRAY);
	}
	//write the image 
	//imwrite("c:\\temp\\temp.jpg", im);

	//mCurrentFrame = Frame(mImGray, timestamp, 
	//	mpIniORBextractor, mpORBVocabulary, 
	//	mK, mDistCoef, mbf, mThDepth);
		
	//string strSettingPath;
	//cv::FileStorage fSettings(strSettingPath, cv::FileStorage::READ);

	int   nFeatures = 1000;   // fSettings["ORBextractor.nFeatures"];
	float fScaleFactor = 1.2; // fSettings["ORBextractor.scaleFactor"];
	int   nLevels = 8;        // fSettings["ORBextractor.nLevels"];
	int   fIniThFAST = 20;    // fSettings["ORBextractor.iniThFAST"];
	int   fMinThFAST = 7;     // fSettings["ORBextractor.minThFAST"];

	ORBextractor* pIniORBextractor = new ORBextractor(2 * nFeatures,
		fScaleFactor, nLevels, fIniThFAST, fMinThFAST);
	
	std::vector<cv::KeyPoint> tKeys;
	cv::Mat tDescriptors;
	
	//operator() to detect feature points
	(*pIniORBextractor)(im, cv::Mat(), tKeys, tDescriptors);

	Keys = tKeys;
	//tDescriptors.copyTo(Descriptors);
	Descriptors = tDescriptors.clone();

	
	////draw feature points
	//Mat colorImg = imread(filepath, CV_LOAD_IMAGE_COLOR);
	//for (int i = 0; i < Keys.size(); i++)
	//{
	//	//circle(colorImg, Point(Keys[i].pt.x, Keys[i].pt.y), 1, CV_RGB(255, 0, 0), 1);
	//	drawMarker(colorImg, Point(Keys[i].pt.x, Keys[i].pt.y), CV_RGB(255, 0, 0), MARKER_CROSS, 2);
	//}
	//imwrite("c:\\temp\\orb-feat.jpg", colorImg);
	//
	
	//for matching
	//vector<int> matches;
	//OrbMatchGeneral()

	delete pIniORBextractor;


	return 0;
}
