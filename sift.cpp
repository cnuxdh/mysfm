


//using namespace std;

#include "sift.hpp"
#include "panorama.hpp"
#include "funcs.hpp"


//corelib
#include "Matrix.h"


//siftdll
#include "siftfeature.h"

//coredll
#include "main.h"


#ifdef WIN32
#include "windows.h"
#include "mmsystem.h"
#endif


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

	double sx = (double)(wd) / (double)(dstWd);
	double sy = (double)(ht) / (double)(dstHt);
	
	//sift feature detect
	struct feature* pFeat = NULL;
	int nFeat;
	nFeat = siftFeatures(filePath, dstHt, dstWd, &pFeat);	
	printf("feature Number: %d \n", nFeat);

	imgFeat.ht = ht;
	imgFeat.wd = wd;

	for(int i=0; i<nFeat; i++)
	{
		stPtFeature feat;
		feat.id         = i;   //pFeat[i].index;
		feat.scl        = pFeat[i].scl;
		feat.ori        = pFeat[i].ori;
		
		feat.x = pFeat[i].x * sx;
		feat.y = pFeat[i].y * sy;

		feat.cx = feat.x - wd*0.5; //normalized coordinate
		feat.cy = ht*0.5 - feat.y; //normalized coordinate       

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

int CSIFT::Detect(IplImage* pImage, ImgFeature& imgFeat)
{
	return 0;
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
	
	#ifdef _WIN32
	//unsigned long t1 = timeGetTime();
	#endif
	
	//int64 t1 = getTickCount();
	Key_Point* featPts = SiftFeaturesFloat(fImage, wd, ht, keynumber);
	printf("Feature Number: %d \n", keynumber);
	
	#ifdef _WIN32
	//unsigned long t2 = timeGetTime();	
	//int64 t2 = getTickCount();
	//printf("time for detection: %lf s \n", (double)(t2-t1)/getTickFrequency() );
	//printf("time for detection: %lf s \n", (double)(t2-t1)/1000 );
  #endif
   
	

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


int CSIFTFloat::Detect(IplImage* pSrc, ImgFeature& imgFeat)
{
	//IplImage* pImage = cvLoadImage(filePath, 0);	
	//IplImage* pSrc = cvLoadImage(filePath, 0);

	int sht = pSrc->height;
	int swd = pSrc->width;

	//image resize, to improve the speed	
	int dstHt, dstWd;
	GetResizeDimension(pSrc->height, pSrc->width, dstHt, dstWd);
	double sx = (double)(swd) / (double)(dstWd);
	double sy = (double)(sht) / (double)(dstHt);
	IplImage* pImage = cvCreateImage(cvSize(dstWd, dstHt), 8, 1);
	cvResize(pSrc, pImage);    
	//cvReleaseImage(&pSrc);

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
	
	#ifdef _WIN32
	//unsigned long t1 = timeGetTime();
	#endif
	
	//int64 t1 = getTickCount();
	Key_Point* featPts = SiftFeaturesFloat(fImage, wd, ht, keynumber);
	printf("Feature Number: %d \n", keynumber);
	
	#ifdef _WIN32
	//unsigned long t2 = timeGetTime();	
	//int64 t2 = getTickCount();
	//printf("time for detection: %lf s \n", (double)(t2-t1)/getTickFrequency() );
	//printf("time for detection: %lf s \n", (double)(t2-t1)/1000 );
  #endif
   
	

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


int CSIFTFloat::Detect(char* filePath, ImgFeature& imgFeat, int maxHt)
{
	//IplImage* pImage = cvLoadImage(filePath, 0);
	
	IplImage* pSrc = cvLoadImage(filePath, 0);
	int sht = pSrc->height;
	int swd = pSrc->width;

	//image resize, to improve the speed	
	int dstHt, dstWd;
	dstHt = maxHt;
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
	
	#ifdef _WIN32
	//unsigned long t1 = timeGetTime();
	#endif
	
	//int64 t1 = getTickCount();
	Key_Point* featPts = SiftFeaturesFloat(fImage, wd, ht, keynumber);
	printf("Feature Number: %d \n", keynumber);
	
	#ifdef _WIN32
	//unsigned long t2 = timeGetTime();	
	//int64 t2 = getTickCount();
	//printf("time for detection: %lf s \n", (double)(t2-t1)/getTickFrequency() );
	//printf("time for detection: %lf s \n", (double)(t2-t1)/1000 );
  #endif
   
	

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

CSIFTPano::CSIFTPano()
{

}
CSIFTPano::~CSIFTPano()
{

}

int CSIFTPano::Detect(char* filePath, ImgFeature& imgFeat, int maxHt)
{
	IplImage* pImage = cvLoadImage(filePath, 0);

	Detect(pImage, imgFeat, maxHt);

	cvReleaseImage(&pImage);

	return 0;
}

int CSIFTPano::Detect(IplImage* pImage, ImgFeature& imgFeat, int maxHt)
{
	//if color, convert to gray image
	IplImage* pGray = cvCloneImage(pImage);
	if(pGray->nChannels==3)
	{
		IplImage* pTemp = cvCreateImage(cvGetSize(pGray), IPL_DEPTH_8U, 1);  
		cvCvtColor(pGray, pTemp, CV_BGR2GRAY);
		cvReleaseImage(&pGray);
		pGray = cvCloneImage(pTemp);
		cvReleaseImage(&pTemp);
	}	

	//the image parameter of original panorama
	double radius = (double)(pGray->width) / (2*PI);
	int pht = pGray->height;
	int pwd = pGray->width;

	//added by xiedonghai, 2017.1.29
	imgFeat.ht = pht;
	imgFeat.wd = pwd;

	//resize the image
	IplImage* pResizeImage = NULL;
	pResizeImage = ResizeImage(pGray, maxHt);
	cvReleaseImage(&pGray);
	
	//split the image
	double R[9] = {1,0,0,0,1,0,0,0,1};
	double T[3] = {0,0,0};
	vector<IplImage*>  projImages;
	vector<CameraPara> camParas;
	PanoToPlanes(pResizeImage, 90, 90, 90, 1, R, T, projImages, camParas);

	//detect the feature
	imgFeat.Clear();
	CFeatureBase* pSift = new CSIFTFloat();
	for(int i=0; i<projImages.size(); i++)
	{
		ImgFeature ifeat;
		pSift->Detect(projImages[i], ifeat);
		if(1)
		{
			char file[256];
			sprintf(file, "c:\\temp\\plane_%d.jpg",i);
			DrawFeatPt(ifeat, projImages[i]);
			cvSaveImage(file, projImages[i]);
		}
				
		double R[9];
		memcpy(R, camParas[i].R, sizeof(double)*9);
		invers_matrix(R, 3);

		double focalLen = camParas[i].focus;
		
		for(int k=0; k<ifeat.GetFeatPtSum(); k++)
		{
			Point2DDouble pt = ifeat.GetCenteredPt(k);

			//convert from perspective to panorama image
			double ipt[3];
			ipt[0] = pt.p[0];
			ipt[1] = pt.p[1];
			ipt[2] = -focalLen;
			double res[3];
			mult(R, ipt, res, 3, 3, 1);

			//from 3D to panorama image point
			double ix,iy;
			GrdToSphere_center(res[0], res[1], res[2], radius, ix, iy);
			stPtFeature fp = ifeat.GetFeatPt(k);
			fp.cx = ix;
			fp.cy = iy;
			fp.x  = ix + pwd*0.5;
			fp.y  = pht*0.5 - iy;
			fp.id = imgFeat.GetFeatPtSum(); //must be added, by xie donghai, 2017.1.16

			imgFeat.AddFeatPt( fp );
		}
	}
	delete pSift;


	//remove the same points



	for(int k=0; k<projImages.size(); k++)
		cvReleaseImage(&projImages[k]);	
	cvReleaseImage(&pResizeImage);

	return 0;
}