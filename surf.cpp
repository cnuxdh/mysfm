
#include "surf.hpp"

#include <iostream>
#include <fstream>
#include <string>

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/features2d/features2d.hpp"
#include "opencv2/nonfree/features2d.hpp"
#include "opencv2/stitching/detail/matchers.hpp"

using namespace cv;
using namespace cv::detail;
using namespace std;


CSURF::CSURF()
{

}

CSURF::~CSURF()
{
}

int CSURF::Detect(char *filePath, char *featurePath)
{
	//load image into Mat
	Mat img = imread(filePath);

	//detection using Stitch module
	Ptr<FeaturesFinder> finder = new SurfFeaturesFinder();
    ImageFeatures feature;
	(*finder)(img, feature);

	cout<<"point number:"<<feature.keypoints.size()<<endl;

	//draw feature points on the image

	return 1;
}

int CSURF::Detect(char* filePath, vector<double>& px, vector<double>& py)
{
	//load image into Mat
	Mat img = imread(filePath);

	//detection using Stitch module
	Ptr<FeaturesFinder> finder = new SurfFeaturesFinder();
	ImageFeatures feature;
	
	const int64 start = getTickCount();
	(*finder)(img, feature);
	const double timeSec = (getTickCount() - start) / getTickFrequency();
	cout << "detect : " << timeSec << " sec" << endl;

	cout<<"point number:"<<feature.keypoints.size()<<endl;
	cout<<"feature discriptor"<<"rows: "<<feature.descriptors.rows<<" cols: "<<feature.descriptors.cols<<endl;

	//output feature points
	int numOfPt = feature.keypoints.size();
	px.resize(numOfPt);
	py.resize(numOfPt);

	for(int i=0; i<numOfPt; i++)
	{
		double x = feature.keypoints[i].pt.x;
		double y = feature.keypoints[i].pt.y;
		px.push_back(x);
		py.push_back(y);
		
		//cout<<feature.descriptors.row(i)<<endl;
	}

	return 1;
}

int CSURF::Detect(char* filePath, vector<PtFeature>& featPts)
{
	//load image into Mat
	Mat img = imread(filePath);

	//detection using Stitch module
	Ptr<FeaturesFinder> finder = new SurfFeaturesFinder();
	ImageFeatures feature;

	const int64 start = getTickCount();
	(*finder)(img, feature);
	const double timeSec = (getTickCount() - start) / getTickFrequency();
	cout << "detect : " << timeSec << " sec" << endl;

	//cv::drawKeypoints()

	int numOfPt = feature.keypoints.size();
	featPts.resize(numOfPt);

	int featDim = feature.descriptors.cols;
	uchar* pBuffer = NULL;
	for(int i=0; i<numOfPt; i++)
	{
		double x = feature.keypoints[i].pt.x;
		double y = feature.keypoints[i].pt.y;

		featPts[i].id = i;
		featPts[i].x = x;
		featPts[i].y = y;

		pBuffer = feature.descriptors.ptr<uchar>(i);
		for(int j=0; j<featDim; j++)
		{
			//double v = feature.descriptors.at<double>(i,j);
			double v = pBuffer[j];
			featPts[i].feat.push_back(v);
		}
	}

	return 1;
}

int CSURF::Detect(char* filePath, ImgFeature& imgFeat)
{

	//load image into Mat
	Mat img = imread(filePath);

	imgFeat.ht = img.rows;
	imgFeat.wd = img.cols;

	//detection using Stitch module
	Ptr<FeaturesFinder> finder = new SurfFeaturesFinder();
	ImageFeatures feature;

	const int64 start = getTickCount();
	(*finder)(img, feature);
	const double timeSec = (getTickCount() - start) / getTickFrequency();
	cout << "detect : " << timeSec << " sec" << endl;


	int numOfPt = feature.keypoints.size();
	imgFeat.featPts.resize(numOfPt);

	int featDim = feature.descriptors.cols;
	uchar* pBuffer = NULL;
	for(int i=0; i<numOfPt; i++)
	{
		double x = feature.keypoints[i].pt.x;
		double y = feature.keypoints[i].pt.y;

		imgFeat.featPts[i].id = i;
		imgFeat.featPts[i].x = x;
		imgFeat.featPts[i].y = y;

		pBuffer = feature.descriptors.ptr<uchar>(i);
		for(int j=0; j<featDim; j++)
		{
			//double v = feature.descriptors.at<double>(i,j);
			double v = pBuffer[j];
			imgFeat.featPts[i].feat.push_back(v);
		}
	}

	return 1;
}


CSURFKeyPoint::CSURFKeyPoint()
{

}

CSURFKeyPoint::~CSURFKeyPoint()
{

}

int CSURFKeyPoint::Detect(char *filePath, vector<PtFeature>& featPts)
{
	//load image into Mat
	Mat image = imread(filePath);

	// vector of keypoints
	std::vector<cv::KeyPoint> keypoints;
	
	// Construct the SURF feature detector object
	SurfFeatureDetector surf(2500.); // threshold 

	// Detect the SURF features
	surf.detect(image,keypoints);

	// Draw the keypoints with scale and orientation information
	Mat featImage = image;
	cv::drawKeypoints(image,       // original image
		keypoints,                 // vector of keypoints
		featImage,                     // the resulting image
		cv::Scalar(255,0,255),   // color of the points
		cv::DrawMatchesFlags::DRAW_RICH_KEYPOINTS); //flag

	imwrite("d:\\surf.jpg", featImage);

	int numOfPt = keypoints.size();
	featPts.resize(numOfPt);
	for(int i=0; i<numOfPt; i++)
	{
		double x = keypoints[i].pt.x;
		double y = keypoints[i].pt.y;
		featPts[i].id = i;
		featPts[i].x = x;
		featPts[i].y = y;
	}

	return 1;
}

