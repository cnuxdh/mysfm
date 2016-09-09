

#ifndef CV_MOSAIC_HPP
#define CV_MOSAIC_HPP

#include "export.hpp"
#include "defines.hpp"

/*
//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
using namespace cv;
*/

class DLL_EXPORT CSeamFinderBase
{
public:
	CSeamFinderBase(){}
	virtual ~CSeamFinderBase(){}

	//c-style interface
	virtual int FindPair(unsigned char* mask1, int ht1, int wd1, 
		unsigned char* mask2, int ht2, int wd2, 
		iPoint corner1, iPoint corner2){return 1;}

	virtual int FindPair(MatrixByte mask1, MatrixByte mask2, iPoint corner1, iPoint corner2){return 1;}

};


class DLL_EXPORT CVoroniSeamFinder: public CSeamFinderBase
{
public:
	CVoroniSeamFinder();
	~CVoroniSeamFinder();
private:
	int FindPair(unsigned char* mask1, int ht1, int wd1, 
				 unsigned char* mask2, int ht2, int wd2, 
				 iPoint corner1, iPoint corner2);
	int FindPair(MatrixByte mask1, MatrixByte mask2, iPoint corner1, iPoint corner2);
};


//////////////////////////////////////////////////////////////////////////
class DLL_EXPORT CBlenderBase
{
public:
	CBlenderBase(){}
	virtual ~CBlenderBase(){}

	//calculate the size of final mosaic
	virtual int Prepair(vector<MatrixByte>& images, vector<iPoint> corners){return 1;}
	//feed images one by one
	virtual int Feed( MatrixByte image, MatrixByte mask, iPoint corner){return 1;}
	//output the final blended image
	virtual int Blend( vector<MatrixByte> images, vector<MatrixByte> masks, vector<iPoint> corners, 
					   MatrixByte& dst, MatrixByte& dstMask){return 1;}

	virtual int Blend(char** filenames, int nFile, double blendWeight, char* outFilePath){return 1;}
	
	virtual int GeoBlend( vector<string> filenames, string outFilePath ){return 1;}
	virtual int GeoBlend( char** filenames, int nFile, char* outFilePath ){return 1;}
	virtual int GeoBlendSingle( char** filenames, int nFile, char* outFilePath ){return 1;}
	//virtual int GeoBlend( vector<Mat> images, char* outFilePath ){return 1;}
};

class DLL_EXPORT CMultiBandBlender: public CBlenderBase
{
public:
	CMultiBandBlender();
	~CMultiBandBlender();
	
	int Prepair(vector<MatrixByte>& images, vector<iPoint> corners);
	int Feed( MatrixByte image, MatrixByte mask, iPoint corner);
	int Blend(vector<MatrixByte> images, vector<MatrixByte> masks, vector<iPoint> corners, 
		MatrixByte& dst, MatrixByte& dstMask);
     
	int Blend(char** filenames, int nFile, double blendWeight, char* outFilePath);
    


	//for GeoTiff mosaic and blend
	int GeoBlend( vector<string> filenames, string outFilePath );
	int GeoBlend( char** filenames, int nFile, char* outFilePath );
	int GeoBlendSingle( char** filenames, int nFile, char* outFilePath );
	//int GeoBlend( vector<Mat> images, char* outFilePath );
	
private:


};



//int GenerateMask(char* filename, Mat& mask, double ratio);

//////////////////////////////////////////////////////////////////////////
// mosaic based on optical flow estimation, written by xdh, 2013.12.18
DLL_EXPORT int DirectMosaic(char** filenames, int nfile, char* mosaicFileName);



//DLL_EXPORT int MosaicWithDEM(char* nvmFile, char* outFile, char* mainPath);
//int GetUInt16GrayImage(char* filename, Mat& image, double ratio);

void AutoLevelRange(int* hist, int nHist, double ratio, int& nMin, int& nMax);
void CalculateVarianceAndMean(int* hist, int nHist, double& mean, double& std);





#endif
