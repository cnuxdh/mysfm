

#include"orbfeat.hpp"
#include"ORBextractor.h"

//opencv
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>


CORBFeat::CORBFeat()
{

}

CORBFeat::~CORBFeat()
{

}

int CORBFeat::Detect(Mat& image, ImgFeature& imgFeat, int maxHt)
{
	//Mat inputImage(image);
	Mat inputImage;
	image.copyTo(inputImage);

	if (inputImage.channels() == 3)
	{
		cvtColor(inputImage, inputImage, CV_RGB2GRAY);
	}
	int sht = inputImage.rows;
	int swd = inputImage.cols;
	if (inputImage.empty())
		return -1;
	
	//printf("%d %d \n", sht, swd);

	//resize the image
	int dHt, dWd;
	double scale = 1.0;
	if (inputImage.rows > maxHt)
	{
		dHt = maxHt;
		scale = double(maxHt) / double(inputImage.rows);
		//printf("scale: %lf \n", scale);

		dWd = scale * double(inputImage.cols);
		resize(inputImage, inputImage, Size(dWd, dHt));
	}

	int   nFeatures = 5000;   // fSettings["ORBextractor.nFeatures"];
	float fScaleFactor = 1.2; // fSettings["ORBextractor.scaleFactor"];
	int   nLevels = 3;        // fSettings["ORBextractor.nLevels"];
	int   fIniThFAST = 20;    // fSettings["ORBextractor.iniThFAST"];
	int   fMinThFAST = 7;     // fSettings["ORBextractor.minThFAST"];

	ORBextractor* pIniORBextractor = new ORBextractor(nFeatures,
		fScaleFactor, nLevels, fIniThFAST, fMinThFAST);

	//operator() to detect feature points
	vector<KeyPoint> mKeys;
	Mat mDescriptors;
	(*pIniORBextractor)(inputImage, cv::Mat(), mKeys, mDescriptors);
	delete pIniORBextractor;

	printf("feature point number: %d \n", mKeys.size());
	
	//convert to ImgFeature
	imgFeat.ht = sht;
	imgFeat.wd = swd;

	for (int i = 0; i<mKeys.size(); i++)
	{
		stPtFeature feat;

		feat.col = mKeys[i].pt.x / scale;
		feat.row = mKeys[i].pt.y / scale;

		feat.x = mKeys[i].pt.x / scale;
		feat.y = mKeys[i].pt.y / scale;

		feat.cx = feat.x - swd*0.5; //normalized coordinate
		feat.cy = sht*0.5 - feat.y; //normalized coordinate        

		feat.trackIdx = -1;
		feat.extra = -1;
		imgFeat.featPts.push_back(feat);
	}

	//added by xdh, 2018.1.2, for ORB descriptor specially
	imgFeat.mDescriptors = mDescriptors;

	return 0;
}

int CORBFeat::Detect(char* filePath, ImgFeature& imgFeat, int maxHt)
{
	Mat mImage = imread(filePath, CV_LOAD_IMAGE_UNCHANGED);
	if (mImage.channels() == 3)
	{
		cvtColor(mImage, mImage, CV_RGB2GRAY);
	}	
	int sht = mImage.rows;
	int swd = mImage.cols;
	if (mImage.empty())
		return -1;
	//printf("%d %d \n", sht, swd);

	//resize the image
	int dHt, dWd;
	double scale = 1.0;
	if (mImage.rows > maxHt)
	{
		dHt = maxHt;
		scale = double(maxHt) / double(mImage.rows);
		//printf("scale: %lf \n", scale);

		dWd = scale * double(mImage.cols);
		resize(mImage, mImage, Size(dWd, dHt));
	}

	int   nFeatures = 5000;   // fSettings["ORBextractor.nFeatures"];
	float fScaleFactor = 1.2; // fSettings["ORBextractor.scaleFactor"];
	int   nLevels = 3;        // fSettings["ORBextractor.nLevels"];
	int   fIniThFAST = 20;    // fSettings["ORBextractor.iniThFAST"];
	int   fMinThFAST = 7;     // fSettings["ORBextractor.minThFAST"];

	ORBextractor* pIniORBextractor = new ORBextractor(nFeatures,
		fScaleFactor, nLevels, fIniThFAST, fMinThFAST);
	
	//operator() to detect feature points
	vector<KeyPoint> mKeys;
	Mat mDescriptors;
	(*pIniORBextractor)(mImage, cv::Mat(), mKeys, mDescriptors);
	delete pIniORBextractor;
	
	printf("feature point number: %d \n", mKeys.size());
	

	//convert to ImgFeature
	imgFeat.ht = sht;
	imgFeat.wd = swd;

	for (int i = 0; i<mKeys.size(); i++)
	{
		stPtFeature feat;
		//feat.id = featPts[i].index;
		//feat.scl = featPts[i].scl;
		//feat.scl_octave = featPts[i].scl_octave;
		//feat.sub_intvl = featPts[i].sub_intvl;
		//feat.ori = featPts[i].ori;
		//feat.key_intvl = featPts[i].key_intvl;
		//feat.key_octave = featPts[i].key_octave;

		feat.col = mKeys[i].pt.x / scale;
		feat.row = mKeys[i].pt.y / scale;

		feat.x = mKeys[i].pt.x / scale;
		feat.y = mKeys[i].pt.y / scale;

		feat.cx = feat.x - swd*0.5; //normalized coordinate
		feat.cy = sht*0.5 - feat.y; //normalized coordinate        

		
		//for (j = 0; j<128; j++)
		//	feat.feat.push_back(featPts[i].descriptor[j]);
		
		feat.trackIdx = -1;
		feat.extra = -1;
		imgFeat.featPts.push_back(feat);
	}

	//added by xdh, 2018.1.2, for ORB descriptor specially
	imgFeat.mDescriptors = mDescriptors;
	
	////resize the feature points
	//for (int i = 0; i < mKeys.size(); i++)
	//{
	//	mKeys[i].pt.x /= scale;
	//	mKeys[i].pt.y /= scale;
	//}
	//imgFeat.mKeys = mKeys;
	
	/*
	//for debug, draw feature points
	Mat colorImg = imread(filePath, CV_LOAD_IMAGE_COLOR); // mImage.clone();
	for (int i = 0; i < mKeys.size(); i++)
	{
		//circle(colorImg, Point(Keys[i].pt.x, Keys[i].pt.y), 1, CV_RGB(255, 0, 0), 1);
		drawMarker(colorImg, Point(mKeys[i].pt.x, mKeys[i].pt.y), CV_RGB(255, 0, 0), MARKER_CROSS, 2);
	}
	imwrite("c:\\temp\\orb-feat.jpg", colorImg);
	*/
	
	return 0;
}


CORBFrame::CORBFrame()
{

}

CORBFrame::~CORBFrame()
{

}

int CORBFrame::Detect(string filePath, int maxHt)
{
	Mat mImage = imread(filePath, CV_LOAD_IMAGE_UNCHANGED);
	if (mImage.channels() == 3)
	{
		cvtColor(mImage, mImage, CV_RGB2GRAY);
	}
	int sht = mImage.rows;
	int swd = mImage.cols;
	if (mImage.empty())
		return -1;
	printf("%d %d \n", sht, swd);

	//resize the image
	int dHt, dWd;
	double scale = 1.0;
	if (mImage.rows > maxHt)
	{
		dHt = maxHt;
		scale = double(maxHt) / double(mImage.rows);
		printf("scale: %lf \n", scale);

		dWd = scale * double(mImage.cols);
		resize(mImage, mImage, Size(dWd, dHt));
	}

	int   nFeatures = 1600;   // fSettings["ORBextractor.nFeatures"];
	float fScaleFactor = 1.2; // fSettings["ORBextractor.scaleFactor"];
	int   nLevels = 8;        // fSettings["ORBextractor.nLevels"];
	int   fIniThFAST = 20;    // fSettings["ORBextractor.iniThFAST"];
	int   fMinThFAST = 7;     // fSettings["ORBextractor.minThFAST"];

	ORBextractor* pIniORBextractor = new ORBextractor(nFeatures,
		fScaleFactor, nLevels, fIniThFAST, fMinThFAST);

	(*pIniORBextractor)(mImage, cv::Mat(), mKeys, mDescriptors);
	delete pIniORBextractor;

	printf("feature point number: %d \n", mKeys.size());

	return 0;
}