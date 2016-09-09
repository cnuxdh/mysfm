

#include "sift.hpp"

//siftdll
#include "siftfeature.h"

//coredll
#include "main.h"


//
#include "windows.h"
#include "mmsystem.h"



#ifdef OPENCV_1X 
	#include "cv.h"
	#include "highgui.h"
	#include "cxcore.h"
#else
	#include "opencv2/core/core.hpp"
	#include "opencv2/imgproc/imgproc_c.h"
	#include "opencv2/highgui/highgui.hpp"
	using namespace cv;
#endif

using namespace std;


CSIFT::CSIFT()
{

}


CSIFT::~CSIFT()
{

}


int CSIFT::Detect(char *filePath, char *featurePath)
{

	return 1;
}

int CSIFT::Detect(char* filePath, ImgFeature& imgFeat)
{    
	IplImage* pImage = cvLoadImage(filePath);
	int ht = pImage->height;
	int wd = pImage->width;

	//resize the image
	int dstHt,dstWd;
	GetResizeDimension(ht, wd, dstHt, dstWd);
	
	//sift feature detect
	struct feature* pFeat = NULL;
	int nFeat;
	nFeat = siftFeatures(filePath, dstHt, dstWd, &pFeat);	
	printf("feature Number: %d \n", nFeat);

	imgFeat.ht = dstHt;
	imgFeat.wd = dstWd;

	for(int i=0; i<nFeat; i++)
	{
		stPtFeature feat;
		feat.id         = i;   //pFeat[i].index;
		feat.scl        = pFeat[i].scl;
		feat.ori        = pFeat[i].ori;
		feat.x = pFeat[i].x;
		feat.y = pFeat[i].y;
		feat.cx = feat.x - dstWd*0.5; //normalized coordinate
		feat.cy = dstHt*0.5 - feat.y; //normalized coordinate        
		for(int j=0; j<128; j++)
			feat.feat.push_back( pFeat[i].descr[j]);
		feat.trackIdx = -1;
		feat.extra = -1;
		imgFeat.featPts.push_back(feat);
	}

	free(pFeat);
	cvReleaseImage(&pImage);

	return 1;
}

int CSIFT::Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat)
{
	//sift feature detect
	struct feature* pFeat = NULL;
	int nFeat;
	nFeat = siftFeatures(filePath, dstHt, dstWd, &pFeat);	
	printf("feature Number: %d \n", nFeat);

	imgFeat.ht = dstHt;
	imgFeat.wd = dstWd;

	for(int i=0; i<nFeat; i++)
	{
		stPtFeature feat;
		feat.id         = pFeat[i].index;
		feat.scl        = pFeat[i].scl;
		feat.ori        = pFeat[i].ori;
		feat.x = pFeat[i].x;
		feat.y = pFeat[i].y;
		feat.cx = feat.x - dstWd*0.5; //normalized coordinate
		feat.cy = dstHt*0.5 - feat.y; //normalized coordinate        
		for(int j=0; j<128; j++)
			feat.feat.push_back( pFeat[i].descr[j]);
		feat.trackIdx = -1;
		feat.extra = -1;
		imgFeat.featPts.push_back(feat);
	}

	free(pFeat);

	return 1;
}

//////////////////////////////////////////////////////////////////////////
CSIFTFloat::CSIFTFloat()
{

}


CSIFTFloat::~CSIFTFloat()
{
}

int CSIFTFloat::Detect(char *filePath, char *featurePath)
{
	return 1;
}

int CSIFTFloat::Detect(char* filePath, ImgFeature& imgFeat)
{
	//IplImage* pImage = cvLoadImage(filePath, 0);
	
	IplImage* pSrc = cvLoadImage(filePath, 0);
	int sht = pSrc->height;
	int swd = pSrc->width;

	//image resize, to improve the speed
	
	int dstHt, dstWd;
	GetResizeDimension(pSrc->height, pSrc->width, dstHt, dstWd);
	double sx = (double)(swd) / (double)(dstWd);
	double sy = (double)(sht) / (double)(dstHt);

	IplImage* pImage = cvCreateImage(cvSize(dstWd, dstHt), 8, 1);
	cvResize(pSrc, pImage);    
	cvReleaseImage(&pSrc);


	/*
	float* fImage;
	int ht,wd;
	IplImageToFloatImage(pImage, &fImage, &ht, &wd);
	*/
	
	int i,j;
	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;
	float* fImage = (float*)malloc(ht*wd*sizeof(float));
	float  scale = 1.0;//1.0/256.0;
	for(j=0; j<ht; j++)
		for(i=0; i<wd; i++)
		{
			fImage[j*wd+i] = float(  (unsigned char)(pImage->imageData[j*scanwd+i]) )*scale;
		}
	
	
	printf("ht: %d wd: %d \n", ht, wd);
	printf("Sift Feature detection ... \n");
	

	int keynumber = 0;
	unsigned long t1 = timeGetTime();
	//int64 t1 = getTickCount();
	Key_Point* featPts = SiftFeaturesFloat(fImage, wd, ht, keynumber);
	printf("Feature Number: %d \n", keynumber);
	unsigned long t2 = timeGetTime();	
	//int64 t2 = getTickCount();
	//printf("time for detection: %lf s \n", (double)(t2-t1)/getTickFrequency() );
	printf("time for detection: %lf s \n", (double)(t2-t1)/1000 );
   
	

	imgFeat.ht = sht;
	imgFeat.wd = swd;
	for(i=0; i<keynumber; i++)
	{
		stPtFeature feat;
		feat.id         = featPts[i].index;
		feat.scl        = featPts[i].scl;
		feat.scl_octave = featPts[i].scl_octave;
		feat.sub_intvl  = featPts[i].sub_intvl;
		feat.ori        = featPts[i].ori;
		feat.key_intvl  = featPts[i].key_intvl;
		feat.key_octave = featPts[i].key_octave;

		feat.col = featPts[i].initl_column*sx;
		feat.row = featPts[i].initl_row*sy;
		
		feat.x = featPts[i].initl_column*sx;
		feat.y = featPts[i].initl_row*sy;

		feat.cx = feat.x - swd*0.5; //normalized coordinate
		feat.cy = sht*0.5 - feat.y; //normalized coordinate        
		
		for(j=0; j<128; j++)
			feat.feat.push_back( featPts[i].descriptor[j]);
		feat.trackIdx = -1;
		feat.extra = -1;
		imgFeat.featPts.push_back(feat);
	}

	delete[] featPts;
	free(fImage);

	cvReleaseImage(&pImage);
	
	return 1;
}

int CSIFTFloat::Detect(char* filePath, int dstHt, int dstWd, ImgFeature& imgFeat)
{
	IplImage* pImage = cvLoadImage(filePath, 0);

	//image resize, to improve the speed
    IplImage* pSmallImg = cvCreateImage(cvSize(dstWd, dstHt), 8, 1);
	cvResize(pImage, pSmallImg);

	int i,j;
	int ht = pSmallImg->height;
	int wd = pSmallImg->width;
	int scanwd = pSmallImg->widthStep;
	float* fImage = (float*)malloc(ht*wd*sizeof(float));
	float  scale = 1.0;//1.0/256.0;
	for(j=0; j<ht; j++)
		for(i=0; i<wd; i++)
		{
			fImage[j*wd+i] = float(  (unsigned char)(pSmallImg->imageData[j*scanwd+i]) )*scale;
		}

	printf("Sift Feature detection ... \n");
	int keynumber = 0;
	Key_Point* featPts = SiftFeaturesFloat(fImage, wd, ht, keynumber);
	printf("Feature Number: %d \n", keynumber);

	imgFeat.ht = ht;
	imgFeat.wd = wd;

	for(i=0; i<keynumber; i++)
	{
		stPtFeature feat;
		feat.id         = featPts[i].index;
		feat.scl        = featPts[i].scl;
		feat.scl_octave = featPts[i].scl_octave;
		feat.sub_intvl  = featPts[i].sub_intvl;
		feat.ori        = featPts[i].ori;
		feat.key_intvl  = featPts[i].key_intvl;
		feat.key_octave = featPts[i].key_octave;
		feat.x = featPts[i].initl_column;
		feat.y = featPts[i].initl_row;
		feat.cx = feat.x - wd*0.5; //normalized coordinate
		feat.cy = ht*0.5 - feat.y; //normalized coordinate        
		for(j=0; j<128; j++)
			feat.feat.push_back( featPts[i].descriptor[j]);
		feat.trackIdx = -1;
		feat.extra = -1;
		imgFeat.featPts.push_back(feat);
	}

	delete[] featPts;
	free(fImage);

	cvReleaseImage(&pImage);
	cvReleaseImage(&pSmallImg);

	return 1;
}

