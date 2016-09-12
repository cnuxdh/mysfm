
#include "blend.hpp"
#include "mosaic.hpp"

//corelib
#include "corelib/commondata.h"
#include "corelib/ImageFunc.h"
#include "corelib/image.h"
#include "Corelib/commonfile.h"
#include "Corelib/CommonFuncs.h"
#include "Corelib/Matrix.h"

//gdal
#include "gdal_priv.h"
#include "ogr_spatialref.h"

//coredll
#include"geotiff.h"

//matrix lib
//#include "matrix/matrix.h"



//opencv
#include "opencv/cv.h"
#include "opencv/cxcore.h"
#include "opencv/highgui.h"

//
#include <vector>
#include <algorithm>
using namespace std;


#define MAX_MASK_SIZE 900


/* generate mask for mosaic, ratio is default set to 1
*/
int GenerateImageMask(char* filename, double rawResolution, double outResolution,
					  unsigned char** pMask, int& mht, int& mwd, int fillValue)
{


	IplImage* pImage = cvLoadImage(filename, 0);
	
	
	double ratio = rawResolution / outResolution;

	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;
	
	mht = ht*ratio;
	mwd = wd*ratio;
	*pMask = (unsigned char*)malloc(mht*mwd);
	memset(*pMask, 0, mht*mwd);
	
	for(int j=0; j<mht; j++)
		for(int i=0; i<mwd; i++)
		{
			int y = j/ratio;
			int x = i/ratio;
			
			int imagevalue = (unsigned char)(pImage->imageData[y*scanwd+x]);
			if(imagevalue>0)
				(*pMask)[j*mwd+i] = fillValue;
		}
	
	cvReleaseImage(&pImage);

	return 0;
}


int ReadRawImage(unsigned char** pBuffer, int& ht, int& wd, char* filename)
{
	FILE* fp = fopen(filename, "rb");
	fread(&ht, sizeof(int), 1, fp);
	fread(&wd, sizeof(int), 1, fp);
	*pBuffer = (unsigned char*)malloc(ht*wd*sizeof(unsigned char));
	fread(*pBuffer, sizeof(unsigned char), ht*wd, fp);
	fclose(fp);

	return 0;
}


int SaveRawImage(unsigned char* pBuffer, int ht, int wd, char* filename)
{
	FILE* fp = fopen(filename, "wb");
	fwrite(&ht, sizeof(int), 1, fp);
	fwrite(&wd, sizeof(int), 1, fp);
	fwrite(pBuffer, sizeof(unsigned char), ht*wd, fp);
	fclose(fp);
	
	return 0;
}


int ResizeImage(unsigned char* pSrc, int sht, int swd,
				unsigned char* pDst, int dht, int dwd)
{
	memset(pDst, 0, dht*dwd);

	double dy = double(dht)/double(sht);
	double dx = double(dwd)/double(swd);

	for(int j=0; j<dht; j++)
		for(int i=0; i<dwd; i++)
		{
			int sj = min(sht-1, int(double(j)/dy) );
			int si = min(swd-1, int(double(i)/dx) );
			pDst[j*dwd+i] = pSrc[sj*swd+si];
		}

	return 0;
}



int SaveToJpg(unsigned char* pBuffer, int ht, int wd, char* filename)
{
	IplImage* pImage = cvCreateImage(cvSize(wd, ht), 8, 1);
	//int ht = pImage->height;
	//int wd = pImage->widthStep;
	int scanwd = pImage->widthStep;
    
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			pImage->imageData[j*scanwd+i] = pBuffer[j*wd+i];
		}
	
	cvSaveImage(filename, pImage);
	cvReleaseImage(&pImage);

	return 0;
}

int GetImageBuffer(char* filename, unsigned char** pBuffer, int& oht, int& owd)
{
	IplImage* pImage = cvLoadImage(filename, 0);

	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;

	oht = ht;
	owd = wd;
	*pBuffer = (unsigned char*)malloc(oht*owd);
	memset(*pBuffer, 0, oht*owd);

	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			//int imageValue =  (unsigned char)(pImage->imageData[j*scanwd+i]);
			//if(imageValue>0) (*pBuffer)[j*wd+i] = 255;
			(*pBuffer)[j*wd+i] = (unsigned char)(pImage->imageData[j*scanwd+i]);
		}

	cvReleaseImage(&pImage);

	return 0;
}

int CalculatePyramidLevel(int minWindowSize)
{
	int nL = 0;
	while( minWindowSize > 16 )
	{
		nL++;
		minWindowSize = minWindowSize >> 1;
	}
	return nL;
}


//dilate and then erode
int MorphCloseImage(unsigned char* pMask, int ht, int wd)
{
	//morphology processing
	IplImage* pMaskImage = cvCreateImage(cvSize(wd,ht), 8, 1);
	int scanWd = pMaskImage->widthStep;
	for(int kj=0; kj<ht; kj++)
		for(int ki=0; ki<wd; ki++)
		{
			pMaskImage->imageData[kj*scanWd+ki] = pMask[kj*wd+ki];
		}
	cvDilate(pMaskImage, pMaskImage);
	cvErode(pMaskImage, pMaskImage);
	//IplConvKernel* element = NULL;
	//element = cvCreateStructuringElementEx( 11,11, 5, 5, CV_SHAPE_RECT, 0 );		
	//cvErode(pMaskImage, pMaskImage, element);
	//cvReleaseStructuringElement(&element);   
	memset(pMask, 0, ht*wd);
	for(int kj=0; kj<ht; kj++)
		for(int ki=0; ki<wd; ki++)
		{
			pMask[kj*wd+ki] = pMaskImage->imageData[kj*scanWd+ki];
		}
	cvReleaseImage(&pMaskImage);
	//////////////////////////////////////////////////////////////////////////	

	return 0;
}


//dilate 
int MorphDilateImage(unsigned char* pMask, int ht, int wd, int nsize)
{
	//morphology processing
	IplImage* pMaskImage = cvCreateImage(cvSize(wd,ht), 8, 1);
	int scanWd = pMaskImage->widthStep;
	for(int kj=0; kj<ht; kj++)
		for(int ki=0; ki<wd; ki++)
		{
			pMaskImage->imageData[kj*scanWd+ki] = pMask[kj*wd+ki];
		}

	//cvDilate(pMaskImage, pMaskImage);
	
	IplConvKernel* element = NULL;
	element = cvCreateStructuringElementEx( nsize,nsize, nsize/2, nsize/2, CV_SHAPE_RECT, 0 );		
	cvDilate(pMaskImage, pMaskImage, element);
	cvReleaseStructuringElement(&element);   

	memset(pMask, 0, ht*wd);
	for(int kj=0; kj<ht; kj++)
		for(int ki=0; ki<wd; ki++)
		{
			pMask[kj*wd+ki] = pMaskImage->imageData[kj*scanWd+ki];
		}
	cvReleaseImage(&pMaskImage);
	//////////////////////////////////////////////////////////////////////////	

	return 0;
}


//dilate 
int MorphErodeImage(unsigned char* pMask, int ht, int wd, int nsize)
{
	//morphology processing
	IplImage* pMaskImage = cvCreateImage(cvSize(wd,ht), 8, 1);
	int scanWd = pMaskImage->widthStep;
	for(int kj=0; kj<ht; kj++)
		for(int ki=0; ki<wd; ki++)
		{
			pMaskImage->imageData[kj*scanWd+ki] = pMask[kj*wd+ki];
		}
		
	IplConvKernel* element = NULL;
	element = cvCreateStructuringElementEx( nsize,nsize, nsize/2, nsize/2, CV_SHAPE_RECT, 0 );		
	cvErode(pMaskImage, pMaskImage, element);
	cvReleaseStructuringElement(&element);   
	
	memset(pMask, 0, ht*wd);
	for(int kj=0; kj<ht; kj++)
		for(int ki=0; ki<wd; ki++)
		{
			pMask[kj*wd+ki] = pMaskImage->imageData[kj*scanWd+ki];
		}
	cvReleaseImage(&pMaskImage);
	//////////////////////////////////////////////////////////////////////////	

	return 0;
}


int    NormalizeImage(unsigned char* pImage, int ht, int wd, double dstStd, double dstMean)
{
	double srcMean = 0;
	double srcStd = 0;
	int hist[256];
	memset(hist, 0, sizeof(int)*256);
	for(int i=0; i<ht*wd; i++)
		hist[pImage[i]]++;

	CalculateVarianceAndMean(hist, 256, srcMean, srcStd);
	printf("%lf %lf \n", srcMean, srcStd);
	
	int nmin = 255;
	int nmax = 0;
	for(int i=0; i<ht*wd; i++)
	{
		//nmin = min(nmin, pImage[i]);
		//nmax = max(nmax, pImage[i]);
		if(nmin>pImage[i])
			nmin = pImage[i];
		if(nmax<pImage[i])
			nmax = pImage[i];
	}	

	//change to the set std and mean
	double sk = dstStd / srcStd;
	double sb = dstMean - sk*srcMean;
	printf("%lf %lf \n", sk, sb);

	for(int ii=0; ii<ht*wd; ii++)
	{
		if(pImage[ii]==0)
			continue;

		int enhanceValue = pImage[ii]; //(double)(pImage[ii]-nmin)/(double)(nmax-nmin)*255;
		//enhanceValue = min(255, max(0,enhanceValue));
		
		//fixed mean and std
		enhanceValue = sk*enhanceValue + sb;
		enhanceValue = min(255, max(0,enhanceValue));

		pImage[ii] = enhanceValue;
	}	
	return 0;
}

/* mask space is different with the mosaic space, 
   the seam line points found in the mask space must be converted to 
   the mosaic space, added by xiedonghai, 2015.11.10
*/
int EnlargeSeamLine(vector<iPoint>& srcSeamPts, int sht, int swd, 
					int dstHt, int dstWd, vector<iPoint>& dstSeamPts)
{
	dstSeamPts.clear();

	unsigned char* pSeamLineEdgeImage = (unsigned char*)malloc(sht*swd);
	memset(pSeamLineEdgeImage, 0, sht*swd);

	for(int i=0; i<srcSeamPts.size(); i++)
	{
		pSeamLineEdgeImage[ srcSeamPts[i].y*swd + srcSeamPts[i].x ] = 1;
	}

	double ry = double(dstHt)/double(sht);
	double rx = double(dstWd)/double(swd);
	double ratio = (rx+ry)*0.5;
	ratio = 1.0/ratio;
	
	for(int j=0; j<dstHt; j++)
		for(int i=0; i<dstWd; i++)
		{
			int sx = i*ratio;
			int sy = j*ratio;
			sx = min(swd-1, sx);
			sy = min(sht-1, sy);
			
			if( pSeamLineEdgeImage[sy*swd+sx]>0 )
			{
				iPoint p;
				p.x = i;
				p.y = j;
				dstSeamPts.push_back(p);
			}
		}

	free(pSeamLineEdgeImage);
	return 0;
}


int FindSeamLine(unsigned short* pMask, int ht, int wd, vector<iPoint>& seamPts )
{
	seamPts.clear();

	for(int j=1; j<ht-1; j++)
		for(int i=1; i<wd-1; i++)
		{
			if( pMask[j*wd+i]==0 )
				continue;

			int nZeros = 0;			
			for(int m=-1; m<=1; m++)
				for(int n=-1; n<=1; n++)
				{
					if(pMask[(j+m)*wd + (i+n)]==0)
					{
						nZeros = 1;
						break;
					}
				}
			if(nZeros>0)
				continue;
			
			int centerValue = pMask[j*wd+i];
			//int off[8] = {-1,1, -wd, -wd-1, -wd+1, wd, wd-1, wd+1};
			int off[4] = {-1,1, -wd, wd};
			int nDiff = 0;
			for(int k=0; k<4; k++)
			{
				int neighborValue = pMask[ j*wd+i +off[k] ];
				if(neighborValue!=centerValue)
				{
					nDiff = 1;
					break;
				}
			}

			if(nDiff>0)
			{
				iPoint pt;
				pt.x = i;
				pt.y = j;
				seamPts.push_back(pt);
			}
		}
	return 0;
}


int ZoomoutSeamLine(vector<iPoint>& inSeamPts, int ht, int wd, int ratio,
	vector<iPoint>& outSeamPts)
{
	outSeamPts.clear();

	int zht = ht / ratio;
	int zwd = wd / ratio;

	unsigned char* pMask = (unsigned char*)malloc(zht*zwd);
	memset(pMask, 0, zht*zwd);

	for(int i=0; i<inSeamPts.size(); i++)
	{
		int zx = inSeamPts[i].x / ratio;
		int zy = inSeamPts[i].y / ratio;

		zx = min(zwd-1, zx);
		zy = min(zht-1, zy);

		if( pMask[zy*zwd+zx] != 0 )
			continue;

		iPoint np;
		np.x = zx;
		np.y = zy;
		outSeamPts.push_back(np);
		pMask[zy*zwd+zx] = 1;
	}

	free(pMask);

	return 0;
}



int DirectBlend(char** filenames, int nFile,
				vector<stGeoInfo> geoArray, 
				double minx, double maxx, double miny, double maxy, 
				double outImageResolution,
				double maskResolution,
				unsigned short* pAllMask, int mht, int mwd,
				double* gainParas, //global intensity correction parameters
				vector<double> cp, //color correction parameters
				char* outFile
				)
{

	if(nFile<=1)
		return -1;

	int oht = (maxy-miny) / outImageResolution;
	int owd = (maxx-minx) / outImageResolution;

	//GDALAllRegister();
	//CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	stGeoInfo geoinfo;
	GetGeoInformation(filenames[0], geoinfo);
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	int zoneNumber = geoinfo.zoneNumber;
	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );      

	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}	
	printf("gdal creating file... \n");
	printf("%d %d \n", owd, oht);

	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );
	if(poDataset==NULL)
	{
		printf("Failed to create file using gdal ! \n");
		return -1;
	}

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = outImageResolution;
	geoTransform[5] = -outImageResolution;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(pszWKT);
	
	double dTime = 55;
	double dStep = dTime / (3*nFile);
	int    pi = 0;

	unsigned char* pMosaicBuffer = (unsigned char*)malloc(oht*owd*sizeof(unsigned char));
	
	
	bool isColorCorrection = false;
    if(cp.size()>0)
		isColorCorrection = true;

	double a1,b1,c1; //exposure 
	double v1,v2;    //vignetting

	double scale = 1; //outImageResolution/maskResolution;

	for(int bandId=0; bandId<3; bandId++)
	{
		memset(pMosaicBuffer, 0, oht*owd*sizeof(unsigned char));
        
		for(int i=0; i<nFile; i++)
		{
			scale = fabs(geoArray[i].dx) / maskResolution;

			if(isColorCorrection)
			{
				a1 = cp[i*9 + 3*bandId ];
				b1 = cp[i*9 + 3*bandId + 1];	
				c1 = cp[i*9 + 3*bandId + 2];
				v1 = cp[cp.size()-2];
				v2 = cp[cp.size()-1];
			}

			printf("band: %d  %s \n", bandId, filenames[i]);
			printf("gain: %lf \n \n ", gainParas[i]);

			//reading image 
			IplImage* pImage = cvLoadImage(filenames[i]);
			int ht = pImage->height;
			int wd = pImage->width;
			int scanWd = pImage->widthStep;

			//locate the rectangle of each file
			stGeoInfo sgeo;
			GetGeoInformation(filenames[i], sgeo);
			int sl = (sgeo.left - minx) / outImageResolution;
			int sr = sl + (sgeo.wd*fabs(sgeo.dx)) / outImageResolution;
			int st = (maxy - sgeo.top) / outImageResolution;
			int sb = st + (sgeo.ht*fabs(sgeo.dy)) / outImageResolution;
			
			//collect pixel index 
			vector<iPoint> ptIndex; 
			ptIndex.clear();
			for(int m=st; m<sb; m++)
				for(int n=sl; n<sr; n++)
				{
					double x = n*outImageResolution;
					double y = m*outImageResolution;	

					//mask space
					int mx = x / maskResolution;
					int my = y / maskResolution;
					mx = min(mx, mwd-1);
					my = min(my, mht-1);

					if( pAllMask[my*mwd+mx] == (i+1) )
					{
						iPoint pi;
						pi.x = n;
						pi.y = m;
						ptIndex.push_back(pi);
					}
				}
				printf("mask point: %d \n", ptIndex.size());

				//fill the image pixel
				for(int m=0; m<ptIndex.size(); m++)
				{
					//mosaic image space
					double x = minx + ptIndex[m].x*outImageResolution;
					double y = maxy - ptIndex[m].y*outImageResolution;

					//single image space
					int ix = (x - geoArray[i].left) / fabs(geoArray[i].dx) ;
					int iy = (geoArray[i].top - y)  / fabs(geoArray[i].dx) ;
					ix = max(0, min(ix, wd-1));
					iy = max(0, min(iy, ht-1));

					double cx = ix-wd*0.5;
					double cy = iy-ht*0.5;
					double radius = sqrt(cx*cx+cy*cy)*scale;				

					//color correction
					if(isColorCorrection)
					{
						int rawValue = (unsigned char)(pImage->imageData[ iy*scanWd + 3*ix + bandId]);
                        
						int correctValue = a1*rawValue + b1*rawValue*rawValue + c1 
							+ v1*rawValue*pow(radius,2) + v2*rawValue*pow(radius,4);

						//int correctValue = a1*rawValue + c1
						//	+ v1*rawValue*pow(radius,2) + v2*rawValue*pow(radius,4);

						correctValue = min(255, correctValue);

						pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = correctValue;
							
					}
					else if (gainParas!=NULL)
					{
						//global instensity correction
						int value = (unsigned char)(pImage->imageData[ iy*scanWd + 3*ix + bandId]);
                        value *= gainParas[i];
						pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = min(255, value);
					}
					else
					{
						pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = 
							(unsigned char)(pImage->imageData[ iy*scanWd + 3*ix + bandId]);
					}					
				}				

			poBand = poDataset->GetRasterBand( 3-bandId );
			poBand->RasterIO(GF_Write, 0, 0, owd, oht, pMosaicBuffer, owd, oht, GDT_Byte, 0, 0);

			cvReleaseImage(&pImage);
			WriteProgressValueToFile( dStep );		
		}
	}

	free(pMosaicBuffer);

	GDALClose( (GDALDatasetH) poDataset );
	
	return 0;
}

int LaplacianBlend(char** filenames, int nFile,
				   vector<stGeoInfo> geoArray, 
				   double minx, double maxx, double miny, double maxy, 
				   double outImageResolution,
				   double maskResolution,
				   int nPyramidLevel,
				   unsigned short* pAllMask, int mht, int mwd,
				   double* gainParas,
				   char* outFile)
{
	if(nFile<=1)
		return -1;

	int oht = (maxy-miny) / outImageResolution;
	int owd = (maxx-minx) / outImageResolution;
		
	//GDALAllRegister();
	//CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	stGeoInfo geoinfo;
	GetGeoInformation(filenames[0], geoinfo);
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	int zoneNumber = geoinfo.zoneNumber;
	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );      
	
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}	
	printf("gdal creating file... \n");
	printf("%d %d \n", owd, oht);
	
	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );
	if(poDataset==NULL)
	{
		printf("Failed to create file using gdal ! \n");
		return -1;
	}

    double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = outImageResolution;
	geoTransform[5] = -outImageResolution;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(pszWKT);
	//poDataset->SetProjection(geoinfo.projectRef);

		
	//prepare for the mosaic image
	PYRAMID_IMAGE_GENERAL<short> mosaicLaplacianPI;
	MallocPyramidGeneral(oht, owd, mosaicLaplacianPI, nPyramidLevel);
	
	
	double dTime = 55;
	double dStep = dTime / (3*nFile);
	int    pi = 0;
	//generate Laplacian pyramid image one by one
	for(int bandId=0; bandId<3; bandId++)
	{
		double dstMean,dstStd;
		for(int i=0; i<nFile; i++)
		{
			printf("band: %d  %s \n", bandId, filenames[i]);

			//reading image 
			IplImage* pImage = cvLoadImage(filenames[i]);

			//resize the image, added by xiedonghai
			double ratio = fabs(geoArray[i].dx) / outImageResolution;
			IplImage* pResizeImage = cvCreateImage(
				cvSize(pImage->width*ratio, pImage->height*ratio), 
				pImage->depth, pImage->nChannels);
			cvResize(pImage, pResizeImage);
			//cvSaveImage("d:\\resize.jpg", pResizeImage);
			cvReleaseImage(&pImage);
			
			int ht = pResizeImage->height;
			int wd = pResizeImage->width;
			int scanWd = pResizeImage->widthStep;
			//unsigned short* pImageBuffer = (unsigned short*)malloc(ht*wd*sizeof(unsigned short));
			unsigned char* pImageBuffer = (unsigned char*)malloc(ht*wd*sizeof(unsigned char));
			for(int kj=0; kj<ht; kj++)
				for(int ki=0; ki<wd; ki++)
				{
					pImageBuffer[kj*wd+ki] = (unsigned char)(pResizeImage->imageData[kj*scanWd + 3*ki + bandId]) ;
					pImageBuffer[kj*wd+ki] = min( 255, int(pImageBuffer[kj*wd+ki]*gainParas[i]) );
				}
			cvReleaseImage(&pResizeImage);
						
			//generate gaussian pyramid
			PYRAMID_IMAGE_GENERAL<unsigned char> pyramidImage; 
			ConstructPyramid(pImageBuffer, ht, wd, pyramidImage, nPyramidLevel);
			free(pImageBuffer);
			//WritePyramidImage(pyramidImage);

			//generate laplacian pyramid
			PYRAMID_IMAGE_GENERAL<short> laplacianPI;// = new PYRAMID_IMAGE_GENERAL();
			ConstructLaplacianPyramid(pyramidImage, laplacianPI);
			FreePyramidImageGeneral(pyramidImage);
			//WritePyramidImage(laplacianPI);

			double imageRatio = 1;
			for(int li=0; li<laplacianPI.nL; li++)
			{
				int lht = mosaicLaplacianPI.pHt[li];
				int lwd = mosaicLaplacianPI.pWd[li];

				int iht = laplacianPI.pHt[li];
				int iwd = laplacianPI.pWd[li];

				//points for file i
				vector<iPoint> ptIndex; 
				ptIndex.clear();
				for(int m=0; m<lht; m++)
					for(int n=0; n<lwd; n++)
					{
						//mosaic image space
						double x = n*outImageResolution*imageRatio;
						double y = m*outImageResolution*imageRatio;	

						//mask space
						int mx = x / maskResolution;
						int my = y / maskResolution;
						mx = min(mx, mwd-1);
						my = min(my, mht-1);

						if( pAllMask[my*mwd+mx] == (i+1) )
						{
							iPoint pi;
							pi.x = n;
							pi.y = m;
							ptIndex.push_back(pi);
						}
					}

					for(int m=0; m<ptIndex.size(); m++)
					{
						//mosaic image space
						double x = minx + ptIndex[m].x*outImageResolution*imageRatio;
						double y = maxy - ptIndex[m].y*outImageResolution*imageRatio;

						//single image space
						int ix = (x - geoArray[i].left) / outImageResolution / imageRatio;//fabs(geoArray[i].dx) / imageRatio ;
						int iy = (geoArray[i].top - y)  / outImageResolution / imageRatio;//fabs(geoArray[i].dx) / imageRatio ;
						ix = max(0, min(ix, iwd-1));
						iy = max(0, min(iy, iht-1));

						//short mvalue =  mosaicLaplacianPI.pLevelImage[li][ ptIndex[m].y*lwd + ptIndex[m].x];

						//if( mvalue != 0 )
						//	mosaicLaplacianPI.pLevelImage[li][ ptIndex[m].y*lwd + ptIndex[m].x] = 
						//		(mvalue + laplacianPI.pLevelImage[li][iy*iwd+ix]) * 0.5;
						//else
						mosaicLaplacianPI.pLevelImage[li][ ptIndex[m].y*lwd + ptIndex[m].x] = 
							laplacianPI.pLevelImage[li][iy*iwd+ix];
					}
					imageRatio *= 2;
			}
			FreePyramidImageGeneral(laplacianPI);
			
			pi++;

			WriteProgressValueToFile( dStep );
		}
		

		printf("seam line point average... \n");
		//deal with seam line
		vector<iPoint> seamPts;
		FindSeamLine(pAllMask, mht, mwd, seamPts);
		//enlarge the seam line
		vector<iPoint> enlargeSeamPts;
		EnlargeSeamLine(seamPts, mht, mwd, oht, owd, enlargeSeamPts);
		seamPts = enlargeSeamPts;

		double imageRatio = 2;
		for(int li=1; li<mosaicLaplacianPI.nL; li++)
		{
			int lht = mosaicLaplacianPI.pHt[li];
			int lwd = mosaicLaplacianPI.pWd[li];

			short* pAverageBuffer = (short*)malloc(lht*lwd*sizeof(short));
			memset(pAverageBuffer, 0, lht*lwd*sizeof(short));
			for(int k=0; k<lht*lwd; k++)
				pAverageBuffer[k] = mosaicLaplacianPI.pLevelImage[li][k];
			
			vector<iPoint> zoomSeamPts;
			ZoomoutSeamLine(seamPts, oht, owd, imageRatio, zoomSeamPts);
			
			for(int k=0; k<zoomSeamPts.size(); k++)
			{
				//for(int km=-1; km<=1; km++)
				//	for(int kn=-1; kn<=1; kn++)
					{
						//int zx = min(lwd-1, max(0, int(zoomSeamPts[k].x)+kn) );
						//int zy = min(lht-1, max(0, int(zoomSeamPts[k].y)+km) );
						
						//int zx = zoomSeamPts[k].x;//+kn;
						//int zy = zoomSeamPts[k].y;//+km;
						
						int zx = min(lwd-2, max(1, int(zoomSeamPts[k].x)) );
						int zy = min(lht-2, max(1, int(zoomSeamPts[k].y)) );
						
						//vector<int> neiValues;
						//neiValues.clear();
						int sv = 0;
						for(int m=-1; m<=1; m++)
							for(int n=-1; n<=1; n++)
							{
								//neiValues.push_back( mosaicLaplacianPI.pLevelImage[li][ (zy+m)*lwd+zx+n ] );
								//int ty = min(lht-1, max(0, zy+m));
								//int tx = min(lwd-1, max(0, zx+n));

								sv += mosaicLaplacianPI.pLevelImage[li][ (zy+m)*lwd + zx+n ];
							}
						//sort(neiValues.begin(), neiValues.end());
						//mosaicLaplacianPI.pLevelImage[li][ zy*lwd+zx ] = neiValues[4];	
						//mosaicLaplacianPI.pLevelImage[li][ zy*lwd+zx ] = sv/9;	
						pAverageBuffer[zy*lwd+zx] = sv/9.0;
					}
			}

			for(int k=0; k<lht*lwd; k++)
				mosaicLaplacianPI.pLevelImage[li][k] = pAverageBuffer[k];
			free(pAverageBuffer);

			imageRatio *= 2;
		}
		

		//WritePyramidImage(mosaicLaplacianPI);

		//reconstruction
		printf("ReconstrucFromLaplacian...\n");
		short* mosaicImage = NULL;
		int rht,rwd;
		ReconstrucFromLaplacian(mosaicLaplacianPI, &mosaicImage, &rht, &rwd );
		
		double imageToMaskRatio = outImageResolution/maskResolution;
		unsigned char* pMosaicImage = (unsigned char*)malloc(rht*rwd);
		for(int i=0; i<rht*rwd; i++)
		{
			//mosaicImage[i] = max(0,  min(mosaicImage[i], 255) );

			if( mosaicImage[i]<0 )
				mosaicImage[i] = 0;
			if( mosaicImage[i]>255 )
				mosaicImage[i] = 0;

			//set the pixel outside the mask to black
			int y = i/rwd;
			int x = i%rwd;
			int my = y*imageToMaskRatio;
			int mx = x*imageToMaskRatio;
			my = min(mht-1, my);
			mx = min(mwd-1, mx);
			if(pAllMask[my*mwd+mx]==0)
			{
				pMosaicImage[i] = 0;
				continue;
			}

			pMosaicImage[i] = mosaicImage[i];
			
			/*
			pMosaicImage[i] = mosaicImage[i];
			if( pMosaicImage[i]<0 )
				pMosaicImage[i] = 0;
			if( pMosaicImage[i]>255 )
				pMosaicImage[i] = 0;
				*/
		}
		
		//char outfile[256];
		//sprintf(outfile, "d:\\blend_%d.jpg", bandId);
		//SaveToJpgGeneral(pMosaicImage, rht, rwd, outfile);
		printf("save band buffer... \n");
		poBand = poDataset->GetRasterBand( 3-bandId );
		poBand->RasterIO(GF_Write, 0, 0, owd, oht, pMosaicImage, owd, oht, GDT_Byte, 0, 0);

		free(pMosaicImage);
		free(mosaicImage);
	}

	FreePyramidImageGeneral(mosaicLaplacianPI);
	
	GDALClose( (GDALDatasetH) poDataset );

	WriteProgressValueToFile(5.0);

	return 0;
}




/* 
   for 64bit, the limit for memory should also be considered, 
   so the dimension of image must be controlled 
   to avoid memory leak, written by xiedonghai
   input:
	srcHt, srcWd: the input ht and wd of image
	nByte: the number of bytes for one pixel
   return value:
	the ratio of input image zoomed out
*/
#define MAX_MEMORY_64BIT 1024
double CalculateRatioFor64bitOS(int srcHt, int srcWd, int nByte)
{
	double ratio = 1.0;	

	double dht = srcHt;
	double dwd = srcWd;

	printf("%d %d %lf \n", srcHt, srcWd, dwd*dht);

	double memWanted = (double)(dht*dwd*nByte) / 1024.0 / 1024.0; // MB

	printf("mem wanted: %lf \n", memWanted);
	if(memWanted>MAX_MEMORY_64BIT)
	{
		ratio = sqrt( (double)(MAX_MEMORY_64BIT) / memWanted );
	}
	return ratio;
}


//function to calculate the gain parameters
int CalculateGainOneByOne(char** maskNames, char** filenames, int nFile,
				  vector<MyRect> vecRect, double* pg)
{
	int nFirstImage = 0;
	int nSumImage = 0;
	
	vector<int> label;
	label.resize(nFile);
	for(int i=0; i<nFile; i++)
		label[i]=0;

	label[nFirstImage] = 1;
	pg[nFirstImage] = 1;

	MyRect curRect = vecRect[nFirstImage];

	//loading mask
	unsigned char* pMask1 = NULL;
	int mht1,mwd1;
	ReadRawImage(&pMask1, mht1, mwd1, maskNames[nFirstImage]);

	//loading image and resize
	IplImage* pImage1 = cvLoadImage(filenames[nFirstImage], 0);
	IplImage* pResizeImage1 = cvCreateImage( cvSize(mwd1, mht1), 8, 1 );
	cvResize(pImage1, pResizeImage1);
	//cvSmooth(pResizeImage1, pResizeImage1, CV_BLUR, 11, 11);
	int scanWd1 = pResizeImage1->widthStep;
	cvReleaseImage(&pImage1);	
	//cvSaveImage("d:\\1.jpg", pResizeImage1);

	nSumImage = 1;
	while(nSumImage<nFile)
	{
		//find the image with maximum intersection area
		int maxArea = 0;
		int nSeleIndex = 0;
		MyRect maxRect;
		for(int i=0; i<nFile; i++)
		{
			if(label[i]>0)
				continue;		

			MyRect iRect;
			if( RectIntersection(curRect, vecRect[i], iRect) < 0 )
				continue;
			int area = (iRect.bottom-iRect.top)*(iRect.right-iRect.left);
			if(area>maxArea)
			{
				area = maxArea;
				nSeleIndex = i;
				maxRect = iRect;
			}
		}
		
		//loading mask
		unsigned char* pMask2 = NULL;
		int mht2, mwd2;
		ReadRawImage( &pMask2, mht2, mwd2, maskNames[nSeleIndex] );

		//loading image and resize
		IplImage* pImage2 = cvLoadImage(filenames[nSeleIndex], 0);
		IplImage* pResizeImage2 = cvCreateImage( cvSize(mwd2, mht2), 8, 1 );
		int scanWd2 = pResizeImage2->widthStep;
		cvResize(pImage2, pResizeImage2);
		//cvSmooth(pResizeImage2, pResizeImage2, CV_BLUR, 11, 11);
		cvReleaseImage(&pImage2);	
		//cvSaveImage("d:\\2.jpg", pResizeImage2);

		double nPixel = 0;
		double nCurAvarageI = 0;
		double nNewAvarageI = 0;

#ifdef _DEBUG
		FILE* fp = fopen("d:\\overlap.txt", "w");
#endif

		for(int m=maxRect.top+5; m<maxRect.bottom-5; m+=2)
			for(int n=maxRect.left+5; n<maxRect.right-5; n+=2)
			{
				int nForeGround = 0;
				
				int index1 = ((m-curRect.top))*scanWd1+(n-curRect.left);
				int index2 = ((m-vecRect[nSeleIndex].top))*scanWd2+(n-vecRect[nSeleIndex].left);                    

				int mindex1 = (m-curRect.top)*mwd1+(n-curRect.left);
				int mindex2 = (m-vecRect[nSeleIndex].top)*mwd2+(n-vecRect[nSeleIndex].left);

				nForeGround += pMask1[mindex1];
				nForeGround += pMask2[mindex2];

				if( nForeGround>300 )
				{
					nPixel ++;
					nCurAvarageI += (unsigned char)(pResizeImage1->imageData[index1]);
					nNewAvarageI += (unsigned char)(pResizeImage2->imageData[index2]);
#ifdef _DEBUG
					fprintf(fp, "%d %d \n", (unsigned char)(pResizeImage1->imageData[index1]),
						(unsigned char)(pResizeImage2->imageData[index2]));
#endif
				}
			}

#ifdef _DEBUG
			fclose(fp);
#endif

			//mosaic
			nCurAvarageI /= nPixel;
			nNewAvarageI /= nPixel;
			double ratio = nCurAvarageI / nNewAvarageI;
			
			MyRect mosaicRect;
			mosaicRect.left   = min(curRect.left,   vecRect[nSeleIndex].left);
			mosaicRect.right  = max(curRect.right,  vecRect[nSeleIndex].right);
			mosaicRect.top    = min(curRect.top,    vecRect[nSeleIndex].top);
			mosaicRect.bottom = max(curRect.bottom, vecRect[nSeleIndex].bottom);

			IplImage* pMosaicImage = cvCreateImage(
				cvSize(mosaicRect.right-mosaicRect.left, mosaicRect.bottom-mosaicRect.top),
				8, 1);
			int mht = pMosaicImage->height;
			int mwd = pMosaicImage->width;
			int mscanWd = pMosaicImage->widthStep;
			memset(pMosaicImage->imageData, 0, mht*mscanWd);
			
			unsigned char* pMosaicMask = (unsigned char*)malloc(mht*mwd);
			memset(pMosaicMask, 0, mht*mwd);

			for(int j=curRect.top; j<curRect.bottom; j++)
				for(int i=curRect.left; i<curRect.right; i++)
				{
					int cy = j-mosaicRect.top;
					int cx = i-mosaicRect.left;
					int y = j-curRect.top;
					int x = i-curRect.left;

					pMosaicImage->imageData[cy*mscanWd+cx] = pResizeImage1->imageData[y*scanWd1+x];	
					pMosaicMask[cy*mwd+cx] = pMask1[y*mwd1+x];
				}
			//cvSaveImage("d:\\mosaic.jpg", pMosaicImage);

			for(int j=vecRect[nSeleIndex].top; j<vecRect[nSeleIndex].bottom; j++)
				for(int i=vecRect[nSeleIndex].left; i<vecRect[nSeleIndex].right; i++)
				{
					int cy = j-mosaicRect.top;
					int cx = i-mosaicRect.left;
					int y  = j-vecRect[nSeleIndex].top;
					int x  = i-vecRect[nSeleIndex].left;

					
					if( (unsigned char)pResizeImage2->imageData[y*scanWd2+x] == 0 )
						continue;
				
					/*
					if( (unsigned char)(pMosaicImage->imageData[cy*mscanWd+cx]) > 0 )
					{
						int mean = ( (unsigned char)(pMosaicImage->imageData[cy*mscanWd+cx]) + 
						         ratio*(unsigned char)pResizeImage2->imageData[y*scanWd2+x] ) * 0.5;

						pMosaicImage->imageData[cy*mscanWd+cx] = mean;
					}
					else
					{
						pMosaicImage->imageData[cy*mscanWd+cx] = ratio*pResizeImage2->imageData[y*scanWd2+x];		
					}*/

					pMosaicImage->imageData[cy*mscanWd+cx] = ratio*pResizeImage2->imageData[y*scanWd2+x];		
					//pMosaicImage->imageData[cy*mscanWd+cx] = pResizeImage2->imageData[y*scanWd2+x];		
					
					if( pMask2[y*mwd2+x]==255 )
						pMosaicMask[cy*mwd+cx] = pMask2[y*mwd2+x];
				}
            cvReleaseImage(&pResizeImage2);
			free(pMask2);
			cvReleaseImage(&pResizeImage1);
			pResizeImage1 = cvCloneImage(pMosaicImage);
			//cvSaveImage("d:\\mosaic1.jpg", pMosaicImage);
			cvReleaseImage(&pMosaicImage);
			//cvSaveImage("d:\\mosaic.jpg", pResizeImage1);
			
			free(pMask1);
			pMask1 = pMosaicMask;
			mht1 = mht;
			mwd1 = mwd;
			scanWd1 = mscanWd;
            curRect = mosaicRect;

			pg[nSeleIndex] = ratio;
			label[nSeleIndex] = 1;
			nSumImage ++;

			//SaveBmp("d:\\mosaicMask.bmp", pMask1, mht1, mwd1);
	}

	return 0;
}

//function to calculate the gain parameters
int CalculateGain(char** maskNames, char** filenames, int nFile, vector<MyRect>& vecRect, double* pg)
{
	////////////////////  gain compensation, added by Xie Donghai,2015.11.26/////////////////
	int* pMN = (int*)malloc(nFile*nFile*sizeof(int)); //number of pixel of overlapped area
	memset(pMN, 0, sizeof(int)*nFile*nFile);
	int* pMI = (int*)malloc(nFile*nFile*sizeof(int)); //mean intensity of overlapped area
	memset(pMI, 0, sizeof(int)*nFile*nFile);
	
	printf("Gain Compensation..... \n");
	for(int j=0; j<nFile;  j++)
	{	
		printf("%d %s \n", j, maskNames[j]);

		//loading mask
		unsigned char* pMask1 = NULL;
		int mht1,mwd1;
		ReadRawImage(&pMask1, mht1, mwd1, maskNames[j]);
		printf("%d %d \n", mht1, mwd1);
		
		//loading image and resize
		//printf("loading image: %s ... \n", filenames[j]);
		IplImage* pImage1 = cvLoadImage(filenames[j], 0);
		//if(pImage1==NULL)
		//	printf("failed to load image ! \n");
		IplImage* pResizeImage1 = cvCreateImage( cvSize(mwd1, mht1), 8, 1 );
		cvResize(pImage1, pResizeImage1);
		//cvSmooth(pResizeImage1, pResizeImage1);//, CV_BLUR, 11, 11);
		int scanWd1 = pResizeImage1->widthStep;
		cvReleaseImage(&pImage1);	

		//cvSaveImage("d:\\1.jpg", pResizeImage1);

		printf("gain calculation.... \n");
		for(int i=j+1; i<nFile; i++)
		{
			MyRect iRect;
			if( RectIntersection(vecRect[i], vecRect[j], iRect) < 0 )
				continue;

			//loading mask
			unsigned char* pMask2 = NULL;
			int mht2, mwd2;
			ReadRawImage( &pMask2, mht2, mwd2, maskNames[i] );
			

			//loading image and resize
			IplImage* pImage2 = cvLoadImage(filenames[i], 0);
			IplImage* pResizeImage2 = cvCreateImage( cvSize(mwd2, mht2), 8, 1 );
			int scanWd2 = pResizeImage2->widthStep;
			cvResize(pImage2, pResizeImage2);
			//cvSmooth(pResizeImage2, pResizeImage2, CV_BLUR, 11, 11);
			cvReleaseImage(&pImage2);	
			//cvSaveImage("d:\\2.jpg", pResizeImage2);

#ifdef _DEBUG
			FILE* fp = fopen("d:\\gain_input.txt", "w");
#endif
			for(int m=iRect.top+3; m<iRect.bottom-3; m++)
				for(int n=iRect.left+3; n<iRect.right-3; n++)
				{
					int nForeGround = 0;
					
					//int index1 = (mht1-(m-vecRect[j].top)-1)*scanWd1+(n-vecRect[j].left);
					//int index2 = (mht2-(m-vecRect[i].top)-1)*scanWd2+(n-vecRect[i].left);
					int index1 = ((m-vecRect[j].top))*scanWd1+(n-vecRect[j].left);
					int index2 = ((m-vecRect[i].top))*scanWd2+(n-vecRect[i].left);   

					int mindex1 = ((m-vecRect[j].top)*mwd1+(n-vecRect[j].left));
					int mindex2 = ((m-vecRect[i].top)*mwd2+(n-vecRect[i].left));

					nForeGround += pMask1[mindex1];
					nForeGround += pMask2[mindex2];
					
					if( nForeGround>10 )
					{
						pMN[j*nFile+i] ++;
						pMI[j*nFile+i] += (unsigned char)(pResizeImage1->imageData[index1]);
						pMI[i*nFile+j] += (unsigned char)(pResizeImage2->imageData[index2]);
#ifdef _DEBUG
						fprintf(fp, "%d %d \n", (unsigned char)(pResizeImage1->imageData[index1]),
							(unsigned char)(pResizeImage2->imageData[index2]));
#endif
					}
				}
#ifdef _DEBUG
			fclose(fp);
#endif

			cvReleaseImage(&pResizeImage2);
			free(pMask2);

			if( pMN[j*nFile+i] != 0 )
			{
				pMI[j*nFile+i] /= pMN[j*nFile+i];
				pMI[i*nFile+j] /= pMN[j*nFile+i];
			}
		}

		cvReleaseImage(&pResizeImage1);
		free(pMask1);
	}
    printf(" \n");

	for(int j=0; j<nFile; j++)
		for(int i=j+1; i<nFile; i++)
		{
			pMN[i*nFile + j] = pMN[j*nFile + i];
		}

#ifdef _DEBUG
	printf("Nij..... \n");
	PrintMatrix(pMN, nFile, nFile);
	printf("Iij..... \n");
	PrintMatrix(pMI, nFile, nFile);
#endif

	double* mA = (double*)malloc(nFile*nFile*sizeof(double));
	memset(mA, 0, nFile*nFile*sizeof(double));
	double* mB = (double*)malloc(nFile*sizeof(double));
	memset(mB, 0, nFile*sizeof(double));
	double sigmaN2 = 100;    //sigma*sigma of normalized intensity error
	double sigmaG2 = 0.01;   //sigma*sigma of gain
	for(int j=0; j<nFile; j++)
	{		
		for(int i=0; i<nFile; i++)
		{
			if(j==i) continue;

			mA[j*nFile+j]  += ( pMI[j*nFile+i]*pMI[j*nFile+i]/sigmaN2 + 1.0/sigmaG2 )*pMN[j*nFile+i];
			mA[j*nFile+i]  = -( pMI[j*nFile+i]*pMI[i*nFile+j]/sigmaN2 )*pMN[j*nFile+i];
			mB[j]          += ( 1.0/sigmaG2 )*pMN[j*nFile+i];			
		}
	}
	//printf("A.... \n");
	//PrintMatrix(mA, nFile, nFile);
	//printf("B.... \n");
	//PrintMatrix(mB, nFile, 1);

	//double* pg = (double*)malloc(nFile*sizeof(double));
	invers_matrix(mA, nFile);
	mult(mA, mB, pg, nFile, nFile, 1);
	printf("gain.... \n");
	for(int i=0; i<nFile; i++)
		printf("%s %lf \n", filenames[i], pg[i]);
	printf("\n");
		 
	free(pMN);
	free(pMI);
	free(mA);
	free(mB);

	return 0;
}

//function to calculate the gain and vignetting parameters
int CalculateGainAndVignettParas(char** maskNames, char** filenames, int nFile,
	vector<MyRect> vecRect, vector<double>& cp)
{
	int nGainPara = 9;
	int nParas = nGainPara*(nFile-1) + 2;  //2 paras for global vignetting, 9 for r,g,b of each file 
    
	double*  A    = (double*)malloc(nParas*sizeof(double));
    double*  AtA  = (double*)malloc(nParas*nParas*sizeof(double));

	double*  sAtL = (double*)malloc(nParas*sizeof(double));
	memset(sAtL, 0, nParas*sizeof(double));

	double*  sAtA = (double*)malloc(nParas*nParas*sizeof(double));
	memset(sAtA, 0, nParas*nParas*sizeof(double));
		
	for(int j=0; j<nFile;  j++)
	{	
		printf("%d ", j);
		//loading mask
		unsigned char* pMask1 = NULL;
		int mht1,mwd1;
		ReadRawImage(&pMask1, mht1, mwd1, maskNames[j]);

		//loading image and resize
		IplImage* pImage1 = cvLoadImage(filenames[j], 1);
		IplImage* pResizeImage1 = cvCreateImage( cvSize(mwd1, mht1), 8, 3 );
		cvResize(pImage1, pResizeImage1);
		int scanWd1 = pResizeImage1->widthStep;
		cvReleaseImage(&pImage1);	
		cvSmooth(pResizeImage1, pResizeImage1, CV_GAUSSIAN);

		//cvSaveImage("d:\\1.jpg", pResizeImage1);

		for(int i=j+1; i<nFile; i++)
		{
			MyRect iRect;
			if( RectIntersection(vecRect[i], vecRect[j], iRect) < 0 )
				continue;

			//loading mask
			unsigned char* pMask2 = NULL;
			int mht2, mwd2;
			ReadRawImage( &pMask2, mht2, mwd2, maskNames[i] );

			//loading image and resize
			IplImage* pImage2 = cvLoadImage(filenames[i], 1);
			IplImage* pResizeImage2 = cvCreateImage( cvSize(mwd2, mht2), 8, 3 );
			int scanWd2 = pResizeImage2->widthStep;
			cvResize(pImage2, pResizeImage2);
			cvReleaseImage(&pImage2);	
			cvSmooth(pResizeImage2, pResizeImage2, CV_GAUSSIAN);
			//cvSaveImage("d:\\2.jpg", pResizeImage2);

#ifdef _DEBUG
			int overlapLen = (iRect.bottom-iRect.top)*(iRect.right-iRect.left);
			unsigned char* pOverLap1 = (unsigned char*)malloc(overlapLen);
			unsigned char* pOverLap2 = (unsigned char*)malloc(overlapLen);
			FILE* fp = fopen("d:\\overlapColor.txt", "w");
#endif

			for(int m=iRect.top; m<iRect.bottom; m+=5)
				for(int n=iRect.left; n<iRect.right; n+=5)
				{
					int nForeGround = 0;
					int x1 = n-vecRect[j].left;
					int y1 = m-vecRect[j].top;
					int x2 = n-vecRect[i].left;
					int y2 = m-vecRect[i].top;
					
					double cx1 = x1 - mwd1*0.5;
					double cy1 = y1 - mht1*0.5;
					double cx2 = x2 - mwd2*0.5;
					double cy2 = y2 - mht2*0.5;

					double d1 = sqrt(cx1*cx1+cy1*cy1);
					double d2 = sqrt(cx2*cx2+cy2*cy2);
					
					//int index1 = y1*scanWd2+x1*3;
					//int index2 = y2*scanWd2+x2*3;

					int mindex1 = (m-vecRect[j].top)*mwd1+(n-vecRect[j].left);
					int mindex2 = (m-vecRect[i].top)*mwd2+(n-vecRect[i].left);

					nForeGround += pMask1[mindex1];
					nForeGround += pMask2[mindex2];
					
					if( nForeGround>300 )
					{
						double sr=0;
						double sg=0;
						double sb=0;
						int sumPixel = 0;

						int offset = 3;

						//calculate the average
						int left1   = max(0, x1-offset);
						int right1  = min(mwd1-1, x1+offset);
						int top1    = max(0, y1-offset);
						int bottom1 = min(mht1-1, y1+offset);
                        for(int jj=top1; jj<=bottom1; jj++)
							for(int ii=left1; ii<=right1; ii++)
							{
								if( pResizeImage1->imageData[jj*scanWd1+ii*3]==0 ||
									pResizeImage1->imageData[jj*scanWd1+ii*3+1]==0 ||
									pResizeImage1->imageData[jj*scanWd1+ii*3+2]==0)
									continue;
									
								sr += (unsigned char)pResizeImage1->imageData[jj*scanWd1+ii*3];
								sg += (unsigned char)pResizeImage1->imageData[jj*scanWd1+ii*3+1];
								sb += (unsigned char)pResizeImage1->imageData[jj*scanWd1+ii*3+2];
								
								sumPixel++;
							}

						if(sumPixel<16)
							continue;

						double r1 = sr / sumPixel;
						double g1 = sg / sumPixel;
						double b1 = sb / sumPixel;

						int left2   = max(0, x2-offset);
						int right2  = min(mwd2-1, x2+offset);
						int top2    = max(0, y2-offset);
						int bottom2 = min(mht2-1, y2+offset);
						sr = 0;
						sg = 0;
						sb = 0;
						sumPixel = 0;
						for(int jj=top2; jj<=bottom2; jj++)
							for(int ii=left2; ii<=right2; ii++)
							{
								if( pResizeImage2->imageData[jj*scanWd2+ii*3]==0 ||
									pResizeImage2->imageData[jj*scanWd2+ii*3+1]==0 ||
									pResizeImage2->imageData[jj*scanWd2+ii*3+2]==0)
									continue;

								sr += (unsigned char)pResizeImage2->imageData[jj*scanWd2+ii*3];
								sg += (unsigned char)pResizeImage2->imageData[jj*scanWd2+ii*3+1];
								sb += (unsigned char)pResizeImage2->imageData[jj*scanWd2+ii*3+2];

								sumPixel++;
							}
						
						if(sumPixel<16)
							continue;

						double r2 = sr / sumPixel;
						double g2 = sg / sumPixel;
						double b2 = sb / sumPixel;


						/*
						double r1 = (unsigned char)(pResizeImage1->imageData[index1]);
						double g1 = (unsigned char)(pResizeImage1->imageData[index1+1]);
						double b1 = (unsigned char)(pResizeImage1->imageData[index1+2]);
						
						double r2 = (unsigned char)(pResizeImage2->imageData[index2]);
						double g2 = (unsigned char)(pResizeImage2->imageData[index2+1]);
						double b2 = (unsigned char)(pResizeImage2->imageData[index2+2]);
						*/
#ifdef _DEBUG
						pOverLap1[ (m-iRect.top)*(iRect.right-iRect.left)+n-iRect.left] = r1;
						pOverLap2[ (m-iRect.top)*(iRect.right-iRect.left)+n-iRect.left] = r2;
						fprintf(fp, "%lf %lf %lf   %lf %lf %lf \n", r1, g1, b1, r2, g2, b2);
						//printf("%lf %lf %lf -- %lf %lf %lf \n", r1, g1, b1, r2, g2, b2);
#endif
												

						//red band
						memset(A, 0, nParas*sizeof(double));
						memset(AtA, 0, nParas*nParas*sizeof(double));
						//A[j*9]   = r1;		A[j*9+1] = r1*r1;		A[j*9+2]=1;
						A[(i-1)*9]   = r2;		A[(i-1)*9+1] = r2*r2;		A[(i-1)*9+2]=1;
						//A[(i-1)*nGainPara] = r2;		A[(i-1)*nGainPara+1]=1;
						A[nParas-2] = r2*pow(d2,2) - r1*pow(d1,2);
						A[nParas-1] = r2*pow(d2,4) - r1*pow(d1,4);	
						mult(A,A,AtA,nParas,1,nParas);
						//PrintMatrix(AtA, nParas, nParas);
						for(int ki=0; ki<nParas*nParas; ki++)
							sAtA[ki] += AtA[ki];
						for(int ki=0; ki<nParas; ki++)
							sAtL[ki] += A[ki]*r1;

						//PrintMatrix(sAtA, nParas, nParas);

						//green 
						memset(A, 0, nParas*sizeof(double));
						memset(AtA, 0, nParas*nParas*sizeof(double));
						//A[j*9+3]   = g1;		A[j*9+4] = g1*g1;		A[j*9+5]=1;
						A[(i-1)*9+3]   = g2;		A[(i-1)*9+4] = g2*g2;		A[(i-1)*9+5]=1;
						//A[(i-1)*nGainPara+2] = g2;		A[(i-1)*nGainPara+3] = 1;
						A[nParas-2] = g2*pow(d2,2) - g1*pow(d1,2);
						A[nParas-1] = g2*pow(d2,4) - g1*pow(d1,4);						
						mult(A,A,AtA,nParas,1,nParas);
						for(int ki=0; ki<nParas*nParas; ki++)
							sAtA[ki] += AtA[ki];
						for(int ki=0; ki<nParas; ki++)
							sAtL[ki] += A[ki]*g1;
						//PrintMatrix(sAtA, nParas, nParas);

						//blue
						memset(A, 0, nParas*sizeof(double));
						memset(AtA, 0, nParas*nParas*sizeof(double));
						//A[j*9+6]   = b1;		A[j*9+7] = b1*b1;		A[j*9+8]=1;
						A[(i-1)*9+6]   = b2;		A[(i-1)*9+7] = b2*b2;		A[(i-1)*9+8]=1;
						//A[(i-1)*nGainPara+4] = b2;			A[(i-1)*nGainPara+5] = 1;
						A[nParas-2] = b2*pow(d2,2) - b1*pow(d1,2);
						A[nParas-1] = b2*pow(d2,4) - b1*pow(d1,4);						
						mult(A,A,AtA,nParas,1,nParas);
						for(int ki=0; ki<nParas*nParas; ki++)
							sAtA[ki] += AtA[ki];	
						for(int ki=0; ki<nParas; ki++)
							sAtL[ki] += A[ki]*b1;

						//PrintMatrix(sAtA, nParas, nParas);
					}
				}				

#ifdef _DEBUG
				//SaveBmp("d:\\overlap1.bmp", pOverLap1, iRect.bottom-iRect.top, iRect.right-iRect.left);
				//SaveBmp("d:\\overlap2.bmp", pOverLap2, iRect.bottom-iRect.top, iRect.right-iRect.left);
				free(pOverLap2);
				free(pOverLap1);
				fclose(fp);
#endif

				cvReleaseImage(&pResizeImage2);
				free(pMask2);				
		}

		cvReleaseImage(&pResizeImage1);
		free(pMask1);
	}

	PrintMatrix(sAtA, nParas, nParas);
	PrintMatrix(sAtL, nParas, 1);
    

	/*
	//calculate the smallest eigenvector of AtA
	double* U  = (double*)malloc(nParas*nParas*sizeof(double));
	double* Vt = (double*)malloc(nParas*nParas*sizeof(double));
	//double* S  = (double*)malloc(nParas*sizeof(double));;
	//dgesvd_driver(nParas, nParas, sAtA, U, S, VT);
	double eps = 0.0000000001;
	int ka = nParas+1;
	bmuav(sAtA, nParas, nParas, U, Vt, eps, ka);
	
	double* V = (double*)malloc(nParas*nParas*sizeof(double));
	transpose(Vt, V, nParas, nParas);
	PrintMatrix(V, nParas, nParas);
	PrintMatrix(sAtA, nParas, nParas);
	//output the last column
	printf("solution.... \n");
	for(int i=0; i<nParas; i++)
	{
		printf("%lf ", V[i*nParas + nParas-1]);
		cp[i] = V[i*nParas + nParas-1];
	}*/

	invers_matrix(sAtA, nParas);
	double* solution = (double*)malloc(nParas*sizeof(double));
	mult(sAtA, sAtL, solution, nParas, nParas, 1);
	
	
	cp.resize(nParas+9);
	cp[0]=1; cp[1]=0; cp[2]=0;
	cp[3]=1; cp[4]=0; cp[5]=0;
	cp[6]=1; cp[7]=0; cp[8]=0;
	

	/*
	cp.resize(nParas+6);
	cp[0]=1; cp[1]=0; 
	cp[2]=1; cp[3]=0; 
	cp[4]=1; cp[5]=0; 
	*/

	for(int i=0; i<nParas; i++)
	{
		//printf("%lf ", solution[i]);
		//cp[i+6] = solution[i];
		cp[i+9] = solution[i];
	}
	
	printf("solution.... \n");
	for(int i=0; i<cp.size(); i++)
	{
		printf("%lf ", cp[i]);
	}

	free(A);
	free(AtA);
	free(sAtA);
	free(sAtL);
	free(solution);

	//free(U);
	//free(Vt);
	//free(V);

	printf(" Finished! \n");

	return 0;
}


//function to calculate the gain and vignetting parameters
int CalculateGainAndVignettParasUsingIntensity(char** maskNames, char** filenames, int nFile,
	vector<MyRect> vecRect, vector<double>& cp)
{

	int nGainPara = 2;
	int nParas = nGainPara*(nFile-1);// + 2; //2 paras for global vignetting, 2 for intensity of each file 
    
	double*  A    = (double*)malloc(nParas*sizeof(double));
    double*  AtA  = (double*)malloc(nParas*nParas*sizeof(double));

	double*  sAtL = (double*)malloc(nParas*sizeof(double));
	memset(sAtL, 0, nParas*sizeof(double));

	double*  sAtA = (double*)malloc(nParas*nParas*sizeof(double));
	memset(sAtA, 0, nParas*nParas*sizeof(double));
		
	for(int j=0; j<nFile;  j++)
	{	
		printf("%d ", j);
		//loading mask
		unsigned char* pMask1 = NULL;
		int mht1,mwd1;
		ReadRawImage(&pMask1, mht1, mwd1, maskNames[j]);

		//loading image and resize
		IplImage* pImage1 = cvLoadImage(filenames[j], 0);
		IplImage* pResizeImage1 = cvCreateImage( cvSize(mwd1, mht1), 8, 1 );
		cvResize(pImage1, pResizeImage1);
		int scanWd1 = pResizeImage1->widthStep;
		cvReleaseImage(&pImage1);	
		cvSmooth(pResizeImage1, pResizeImage1, CV_BLUR, 11, 11);
		cvSaveImage("d:\\1.jpg", pResizeImage1);

		for(int i=j+1; i<nFile; i++)
		{
			MyRect iRect;
			if( RectIntersection(vecRect[i], vecRect[j], iRect) < 0 )
				continue;

			//loading mask
			unsigned char* pMask2 = NULL;
			int mht2, mwd2;
			ReadRawImage( &pMask2, mht2, mwd2, maskNames[i] );

			//loading image and resize
			IplImage* pImage2 = cvLoadImage(filenames[i], 0);
			IplImage* pResizeImage2 = cvCreateImage( cvSize(mwd2, mht2), 8, 1 );
			int scanWd2 = pResizeImage2->widthStep;
			cvResize(pImage2, pResizeImage2);
			cvReleaseImage(&pImage2);	
			cvSmooth(pResizeImage2, pResizeImage2, CV_BLUR, 11, 11);
			cvSaveImage("d:\\2.jpg", pResizeImage2);

#ifdef _DEBUG
			int overlapLen = (iRect.bottom-iRect.top)*(iRect.right-iRect.left);
			unsigned char* pOverLap1 = (unsigned char*)malloc(overlapLen);
			unsigned char* pOverLap2 = (unsigned char*)malloc(overlapLen);
			FILE* fp = fopen("d:\\overlapColor.txt", "w");
#endif

			for(int m=iRect.top; m<iRect.bottom; m+=3)
				for(int n=iRect.left; n<iRect.right; n+=3)
				{
					int nForeGround = 0;
					int x1 = n-vecRect[j].left;
					int y1 = m-vecRect[j].top;
					int x2 = n-vecRect[i].left;
					int y2 = m-vecRect[i].top;
					
					double cx1 = x1 - mwd1*0.5;
					double cy1 = y1 - mht1*0.5;
					double cx2 = x2 - mwd2*0.5;
					double cy2 = y2 - mht2*0.5;

					double d1 = sqrt(cx1*cx1+cy1*cy1);
					double d2 = sqrt(cx2*cx2+cy2*cy2);
					
					int index1 = y1*scanWd1+x1;
					int index2 = y2*scanWd2+x2;

					int mindex1 = (m-vecRect[j].top)*mwd1+(n-vecRect[j].left);
					int mindex2 = (m-vecRect[i].top)*mwd2+(n-vecRect[i].left);

					nForeGround += pMask1[mindex1];
					nForeGround += pMask2[mindex2];
					
					if( nForeGround>300 )
					{						
						int i1,i2;

						i1 = (unsigned char)(pResizeImage1->imageData[index1]);
						i2 = (unsigned char)(pResizeImage2->imageData[index2]);					

#ifdef _DEBUG
						pOverLap1[ (m-iRect.top)*(iRect.right-iRect.left)+n-iRect.left] = i1;
						pOverLap2[ (m-iRect.top)*(iRect.right-iRect.left)+n-iRect.left] = i2;
						fprintf(fp, "%d %d \n", i1, i2);
#endif
												
						memset(A, 0, nParas*sizeof(double));
						memset(AtA, 0, nParas*nParas*sizeof(double));
						A[(i-1)*nGainPara] = i2;		A[(i-1)*nGainPara+1]=1;
						//A[nParas-2] = i2*pow(d2,2) - i1*pow(d1,2);
						//A[nParas-1] = i2*pow(d2,4) - i1*pow(d1,4);						
						
						mult(A,A,AtA,nParas,1,nParas);
						//PrintMatrix(AtA, nParas, nParas);
						for(int ki=0; ki<nParas*nParas; ki++)
							sAtA[ki] += AtA[ki];
						for(int ki=0; ki<nParas; ki++)
							sAtL[ki] += A[ki]*i1;

						//PrintMatrix(sAtA, nParas, nParas);

					}
				}				

#ifdef _DEBUG
				SaveBmp("d:\\overlap1.bmp", pOverLap1, iRect.bottom-iRect.top, iRect.right-iRect.left);
				SaveBmp("d:\\overlap2.bmp", pOverLap2, iRect.bottom-iRect.top, iRect.right-iRect.left);
				free(pOverLap2);
				free(pOverLap1);
				fclose(fp);
#endif

				cvReleaseImage(&pResizeImage2);
				free(pMask2);				
		}

		cvReleaseImage(&pResizeImage1);
		free(pMask1);
	}

	PrintMatrix(sAtA, nParas, nParas);
	PrintMatrix(sAtL, nParas, 1);
    

	invers_matrix(sAtA, nParas);
	double* solution = (double*)malloc(nParas*sizeof(double));
	mult(sAtA, sAtL, solution, nParas, nParas, 1);
	
	cp.resize(nParas+2);
	cp[0]=1; cp[1]=0; 
	
	for(int i=0; i<nParas; i++)
	{
		//printf("%lf ", solution[i]);
		cp[i+2] = solution[i];
	}
	
	printf("solution.... \n");
	for(int i=0; i<cp.size(); i++)
	{
		printf("%lf ", cp[i]);
	}

	free(A);
	free(AtA);
	free(sAtA);
	free(sAtL);
	free(solution);

	//free(U);
	//free(Vt);
	//free(V);

	printf(" Finished! \n");

	return 0;
}

/*  calculate the color correction parameters 
*/
int CalculateColorCorrectionParas(char** maskNames, char** filenames, int nFile,
								  vector<MyRect> vecRect)
{
	int referenceImageIndex = 0;

	//loading mask
	unsigned char* pMaskRef = NULL;
	int refHt,refWd;
	ReadRawImage(&pMaskRef, refHt, refWd, maskNames[referenceImageIndex]);

	//loading image and resize
	IplImage* pImage = cvLoadImage(filenames[referenceImageIndex], 3);
	IplImage* pImageRef = cvCreateImage( cvSize(refWd, refHt), 8, 3 );
	cvResize(pImage, pImageRef);
	cvReleaseImage(&pImage);

	MyRect refRect = vecRect[referenceImageIndex];

	//find the intersection image
	int nLeftImage = nFile-1;

	vector<int> visitLabel;
	visitLabel.resize(nFile);
	visitLabel[referenceImageIndex]=1;

	while(nLeftImage>0)
	{
		double maxArea = -1;
		int    selectIndex = -1;
		MyRect overlapRect;
		for(int i=0; i<nFile; i++)
		{
			if(visitLabel[i]==1)
				continue;
			
			MyRect iRect;
			if( RectIntersection(refRect, vecRect[i], iRect) < 0 )
				continue;

			double area = (iRect.bottom-iRect.top)*(iRect.right-iRect.left);
			if(area>maxArea)
			{
				maxArea=area;
				selectIndex = i;
				overlapRect = iRect;
			}
		}

		if(selectIndex<0)
			break;

		//reading selected mask and image
		//loading mask
		unsigned char* pSelectMask = NULL;
		int selectHt,selectWd;
		ReadRawImage(&pSelectMask, selectHt, selectWd, maskNames[selectIndex]);
		//loading image and resize
		IplImage* pImage = cvLoadImage(filenames[referenceImageIndex], 3);
		IplImage* pImageSelect = cvCreateImage( cvSize(selectWd, selectHt), 8, 3 );
		cvResize(pImage, pImageSelect);
		cvReleaseImage(&pImage);

		//calculate the parameters based on overlap 
		for(int j=overlapRect.top; j<overlapRect.bottom; j++)
			for(int i=overlapRect.left; i<overlapRect.right; i++)
			{
				int refx = i-refRect.left;
				int refy = j-refRect.top;
				int cx = i-vecRect[selectIndex].left;
				int cy = j-vecRect[selectIndex].top;
				
				if( pMaskRef[refy*refWd]>100 && pSelectMask[cy*selectWd+cx]>100 )
				{
					double rr,rg,rb,cr,cg,cb;	
					rr = pImageRef->imageData[refy*pImageRef->widthStep + refx*3];
					rg = pImageRef->imageData[refy*pImageRef->widthStep + refx*3 + 1];
					rb = pImageRef->imageData[refy*pImageRef->widthStep + refx*3 + 2];

					cr = pImageSelect->imageData[cy*pImageSelect->widthStep + cx*3];
					cg = pImageSelect->imageData[cy*pImageSelect->widthStep + cx*3 + 1];
					cb = pImageSelect->imageData[cy*pImageSelect->widthStep + cx*3 + 2];


				}
			}
		
		//update the reference image and mask



		nLeftImage--;
	}

	return 0;
}


/*
int GetGeoFileInfo( char** filenames, int nFile, 
					vector<stGeoInfo>& geoArray,
					double& minx, double& maxx, double& miny, double& maxy)
{
	//reading the geo-infomation from the geotiff files
	printf("reading the geo-infomation from the images \n");
	//vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);
	for(int i=0; i<nFile; i++)
	{	
		//printf("%s \n", filenames[i]);
		//GetGeoInformation( filenames[i], geoArray[i]);
		
		char geofile[256];
		strcpy(geofile, filenames[i]);
		int nLen = strlen(geofile);
		strcpy(geofile+nLen-4, "geo");
		
	}

	//double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		//resolution
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);

		//position
		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top  - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

	return 0;
}
*/

int GetImageGeoInfo(char** filenames, int nFile, 
					vector<stGeoInfo>& geoArray,
					double& minx, double& maxx, double& miny, double& maxy)
{
	//reading the geo-infomation from the geotiff files
	printf("reading the geo-infomation from the images \n");
	//vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);
	for(int i=0; i<nFile; i++)
	{	
		//printf("%s \n", filenames[i]);
		//GetGeoInformation( filenames[i], geoArray[i]);
		
		char geofile[256];
		strcpy(geofile, filenames[i]);
		int nLen = strlen(geofile);
		strcpy(geofile+nLen-4, "geo");
		ReadGeoFile(geofile, geoArray[i]);
	}

	//double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		//resolution
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);

		//position
		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top  - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

	return 0;
}


MyRect UnionRect(MyRect rect1, MyRect rect2)
{
	MyRect urect;
	urect.left   = min( rect1.left,   rect2.left );
    urect.right  = max( rect1.right,  rect2.right );
	urect.top    = min( rect1.top,    rect2.top );
	urect.bottom = max( rect1.bottom, rect2.bottom );
	return urect;
}


int SeamlinesFindingVoronoiOneByOne(char** maskNames, int nFile, 
							vector<MyRect> vecRect, int oht, int owd)
{

	//the whole mask
	MyRect allRect;
	allRect.top = 0;
	allRect.bottom = oht;
	allRect.left = 0;
	allRect.right = owd;
	unsigned short* pAllMask = (unsigned short*)malloc(oht*owd*sizeof(unsigned short)); //new unsigned short[oht*owd];
	memset(pAllMask, 0, oht*owd*sizeof(unsigned short));
	
	vector<int> label;
	label.resize(nFile);
	for(int i=0; i<nFile; i++)
		label[i] = 0;

	int nFirst = 0; //the index of first image
	label[nFirst] = 1;

	MyRect curRect = vecRect[nFirst];
	unsigned char* pMask1 = NULL;
	int mht1,mwd1;
	ReadRawImage(&pMask1, mht1, mwd1, maskNames[nFirst]);

	iPoint p1;		
	p1.x = curRect.left;
	p1.y = curRect.top;

	FillMaskValue(pAllMask, oht, owd, allRect, curRect, pMask1, mht1, mwd1, nFirst+1);

	CSeamFinderBase* pSeamFinder = new CVoroniSeamFinder();
	for(int i=1; i<nFile; i++ )
	{
		//find the image with maximum overlap 
		double maxArea = 0;
		int    nSele = -1;
		MyRect interSectRect;
		for(int k=1; k<nFile; k++)
		{	
			if( label[k]>0 )
				continue;

			MyRect iRect;
			if( RectIntersection(curRect, vecRect[k], iRect)<0 )
				continue;

			double area = (iRect.right-iRect.left)*(iRect.bottom-iRect.top);
			if(area>maxArea)
			{
				maxArea = area;
				nSele = k;
				interSectRect = iRect;
			}
		}
		if(nSele<0)
			break;

		label[nSele] = 1;

		//finding seam line between two images		
		unsigned char* pMask2 = NULL;
		int mht2,mwd2;
		ReadRawImage(&pMask2, mht2, mwd2, maskNames[nSele]);
		iPoint p2;		
		p2.x = vecRect[nSele].left;
		p2.y = vecRect[nSele].top;     
	

#ifdef _DEBUG
		SaveToJpg(pMask1, mht1, mwd1, "d:\\mask1-before.jpg");
		SaveToJpg(pMask2, mht2, mwd2, "d:\\mask2-before.jpg");
#endif
		
		pSeamFinder->FindPair(pMask1, mht1, mwd1, pMask2, mht2, mwd2, p1, p2);
		

		/*
		//clip the overlap area
		MyRect clipRect;
		clipRect.left   = interSectRect.left - 36;
		clipRect.top    = interSectRect.top - 36;
		clipRect.right  = interSectRect.right + 36;
		clipRect.bottom = interSectRect.bottom + 36;
		clipRect = RectIntersection(curRect, clipRect);
		int clipHt = clipRect.bottom - clipRect.top+1;
		int clipWd = clipRect.right  - clipRect.left+1;
		unsigned char* pClipMask = (unsigned char*)malloc(clipHt&clipWd);
		memset(pClipMask, 0, clipHt*clipWd);
		ReadRectValue(pClipMask, clipHt, clipWd, pMask1, mht1, mwd1, clipRect, curRect);
		//pSeamFinder->FindPair(pClipMask, clipHt, clipWd, pMask2, mht2, mwd2, p1, p2);
		*/

#ifdef _DEBUG
		SaveToJpg(pMask1, mht1, mwd1, "d:\\mask1.jpg");
		SaveToJpg(pMask2, mht2, mwd2, "d:\\mask2.jpg");
#endif
		
		FillMaskValue(pAllMask, oht, owd, allRect, vecRect[nSele], pMask2, mht2, mwd2, nSele+1);
	
		//union of two mask
		MyRect uRect = UnionRect(curRect, vecRect[nSele]);		
		int uHt = (uRect.bottom-uRect.top);
		int uWd = (uRect.right-uRect.left);
		unsigned char* pCurMask = (unsigned char*)malloc( uHt*uWd );
        memset(pCurMask, 0, uHt*uWd);
		FillMaskValue(pCurMask, uHt, uWd, uRect, curRect, pMask1, mht1, mwd1, 255);
		FillMaskValue(pCurMask, uHt, uWd, uRect, vecRect[nSele], pMask2, mht2, mwd2, 255);

		//update 
		free(pMask1);
		pMask1 = pCurMask;
		mht1 = uHt;
		mwd1 = uWd;
		curRect = uRect;
		p1.x = curRect.left;
		p1.y = curRect.top;
	
#ifdef _DEBUG
		SaveBmpGeneral(pAllMask, oht, owd, "d:\\allmask.bmp");
#endif
	
	}

	

	return 0;
}



int SeamlinesFindingVoronoi(char** maskNames, int nFile, 
							vector<MyRect> vecRect)
{
	printf("Finding seam lines by voronoi... \n");

	CSeamFinderBase* pSeamFinder = new CVoroniSeamFinder();
	
#ifdef _DEBUG
	//int maskHt = ;
	//int maskwd = ;
	//unsigned short* pAllMask = (unsigned short*)malloc(maskHt*maskWd*sizeof(unsigned short));
	//memset(pAllMask, 0, maskHt*maskWd*sizeof(unsigned short));
#endif

	for(int i=0; i<nFile; i++)
	{
		printf(" %d ", i);

		unsigned char* pMask1 = NULL;
		int mht1,mwd1;
		ReadRawImage(&pMask1, mht1,mwd1, maskNames[i]);

		iPoint p1;
		p1.x = vecRect[i].left;
		p1.y = vecRect[i].top;

		for(int j=i+1; j<nFile; j++)
		{						
			unsigned char* pMask2 = NULL;
			int mht2,mwd2;
			ReadRawImage(&pMask2, mht2, mwd2, maskNames[j]);

			iPoint p2;
			p2.x = vecRect[j].left;
			p2.y = vecRect[j].top;

			pSeamFinder->FindPair(pMask1, mht1, mwd1, pMask2, mht2, mwd2, p1, p2);

			SaveRawImage(pMask2, mht2, mwd2, maskNames[j]);

			free(pMask2);
		}

#ifdef _DEBUG
		char maskfile[256];
		strcpy(maskfile, maskNames[i]);
		strcat(maskfile, ".jpg");
		SaveToJpg(pMask1, mht1, mwd1, maskfile);
#endif

		SaveRawImage(pMask1, mht1, mwd1, maskNames[i]);
        
		free(pMask1);
	}

	delete pSeamFinder;	

	return 0;
}



int MultiBandBlend(char** filenames, char** maskNames, int nFile,
				   vector<stGeoInfo> geoArray,
				   vector<MyRect>    vecRect,
				   double outImageResolution,
				   double minx, double maxx, double miny, double maxy,
				   char* outFile)
{
	if(nFile<=1)
		return -1;

	int oht = (maxy-miny) / outImageResolution;
	int owd = (maxx-minx) / outImageResolution;

	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	stGeoInfo geoinfo;
	GetGeoInformation(filenames[0], geoinfo);
	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	int zoneNumber = geoinfo.zoneNumber;
	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );      

	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}	
	printf("gdal creating file... \n");
	printf("%d %d \n", owd, oht);

	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );
	if(poDataset==NULL)
	{
		printf("Failed to create file using gdal ! \n");
		return -1;
	}

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = outImageResolution;
	geoTransform[5] = -outImageResolution;
	poDataset->SetGeoTransform(geoTransform);
	//poDataset->SetProjection(pszWKT);


	short* pMosaicImage = (short*)malloc(oht*owd*sizeof(short));
	memset(pMosaicImage, 0, oht*owd*sizeof(short));
	//unsigned char* pMosaicImage = (unsigned char*)malloc(oht*owd*sizeof(unsigned char));
	//memset(pMosaicImage, 0, oht*owd);

	int    nLevel = 5;
	double sigma  = 11; 
    
	int blockSize = 2000;
	for(int bandId=0; bandId<3; bandId++)
	{
		memset(pMosaicImage, 0, oht*owd*sizeof(short));
		
		int by = oht/blockSize + 1;
		int bx = owd/blockSize + 1;

		for(int j=0; j<by; j++)
			for(int i=0; i<bx; i++)
			{
				MyRect rect;
				rect.left   = i*blockSize;
				rect.right  = min(owd-1, (i+1)*blockSize);
				rect.top    = j*blockSize;
				rect.bottom = min(oht-1, (j+1)*blockSize);

				//finding the images overlap with the current block
				vector<int>    imageIndex;
				vector<MyRect> vecIRects;
				for(int k=0; k<nFile; k++)
				{
					MyRect iRect;
					if( RectIntersection(rect, vecRect[k], iRect) >= 0 )
					{
						vecIRects.push_back( vecRect[k] );
						imageIndex.push_back(k);
					}
				}
				
				//loading the original images and resize
				vector<IplImage*> vecImagePointer;
				vector<IplImage*> vecMaskPointer;
				vecImagePointer.resize( imageIndex.size() );
				vecMaskPointer.resize( imageIndex.size() );
				for(int k=0; k<imageIndex.size(); k++)
				{
					int imageId = imageIndex[k];

					//loading mask
					unsigned char* pMask;
					int mht,mwd;
					ReadRawImage(&pMask, mht,mwd, maskNames[imageId]);
					IplImage* pMaskImage = cvCreateImage(cvSize(mwd,mht),8,1);
					int scanWd = pMaskImage->widthStep;
					for(int kj=0; kj<mht; kj++)
						for(int ki=0; ki<mwd; ki++)
							pMaskImage->imageData[kj*scanWd+ki] = pMask[kj*mwd+ki];
					free(pMask);
					vecMaskPointer[k] = pMaskImage;

					//loading image 
					IplImage* pImage = cvLoadImage( filenames[imageId] );
					//resize the image
					//double ratio = fabs(geoArray[i].dx) / outImageResolution;
					IplImage* pResizeImage = cvCreateImage(
						cvSize(mwd, mht), pImage->depth, pImage->nChannels);
					cvResize(pImage, pResizeImage);
					//cvSaveImage("d:\\resize.jpg", pResizeImage);
					cvReleaseImage(&pImage);
					vecImagePointer[k] = pResizeImage;					
				}

				//blending
				for(int k=nLevel-1; k>=0; k--)
				{
					vector<IMAGE_GENERAL*> vecSIP;
					vector<IplImage*> vecSMP;

					vecSIP.resize(imageIndex.size());
					vecSMP.resize(imageIndex.size());

					int currentSigma = sigma*k;
					if(currentSigma%2==0)
						currentSigma ++;

					//smooth mask images
					for(int ki=0; ki<vecMaskPointer.size(); ki++)
					{
						vecSMP[ki] = cvCloneImage( vecMaskPointer[ki]);
						cvSmooth(vecMaskPointer[ki], vecSMP[ki], CV_GAUSSIAN, currentSigma);
						//cvSaveImage("d:\\smoothMask.jpg", vecSMP[ki]);
					}

					if( k==(nLevel-1) )
					{						
						for(int ki=0; ki<vecImagePointer.size(); ki++)
						{
							
							int ht = vecImagePointer[ki]->height;
							int wd = vecImagePointer[ki]->width;
							int swd = vecImagePointer[ki]->widthStep;
							IplImage* pSmooth = cvCloneImage( vecImagePointer[ki]);
							cvSmooth(vecImagePointer[ki], pSmooth, CV_GAUSSIAN, currentSigma);
							cvSaveImage("d:\\smoothProvious.jpg", pSmooth);
							
							vecSIP[ki] = (IMAGE_GENERAL*)malloc(sizeof(IMAGE_GENERAL));
							vecSIP[ki]->pImage = (short*)malloc(ht*wd*sizeof(short));
							vecSIP[ki]->nHt = ht;
							vecSIP[ki]->nWd = wd;
							memset(vecSIP[ki]->pImage, 0, ht*wd*sizeof(short));
						
							for(int mj=0; mj<ht; mj++)
								for(int mi=0; mi<wd; mi++)
								{
									vecSIP[ki]->pImage[mj*wd + mi] = 
										(unsigned char)(pSmooth->imageData[mj*swd + mi*3 + bandId]);
								}
							cvReleaseImage(&pSmooth);			
							SaveBmp("d:\\smooth.bmp", vecSIP[ki]->pImage, ht, wd );							
						}						
					}
					else
					{
						for(int ki=0; ki<vecImagePointer.size(); ki++)
						{							
							IplImage* pSmoothCurrent  = cvCloneImage( vecImagePointer[ki]);
							IplImage* pSmoothPrevious = cvCloneImage( vecImagePointer[ki]);
							
							cvSmooth(vecImagePointer[ki], pSmoothCurrent, CV_GAUSSIAN, currentSigma);
							int preSigma = currentSigma+sigma;
							if(preSigma%2==0)
								preSigma++;
							cvSmooth(vecImagePointer[ki], pSmoothPrevious, CV_GAUSSIAN, preSigma);
							
							int ht = pSmoothCurrent->height;
							int wd = pSmoothCurrent->width;
							int swd = pSmoothCurrent->widthStep;

							vecSIP[ki] = (IMAGE_GENERAL*)malloc(sizeof(IMAGE_GENERAL));
							vecSIP[ki]->pImage = (short*)malloc(ht*wd*sizeof(short));
							vecSIP[ki]->nHt = ht;
							vecSIP[ki]->nWd = wd;
							memset(vecSIP[ki]->pImage, 0, ht*wd*sizeof(short));

							for(int mj=0; mj<ht; mj++)
								for(int mi=0; mi<wd; mi++)
								{
									vecSIP[ki]->pImage[mj*wd + mi] = 
										(unsigned char)(pSmoothCurrent->imageData[mj*swd + mi*3 + bandId])-
										(unsigned char)(pSmoothPrevious->imageData[mj*swd + mi*3 + bandId]);
								}

							cvReleaseImage(&pSmoothCurrent);
							cvReleaseImage(&pSmoothPrevious);
						}						
					}

					for(int kj=rect.top; kj<rect.bottom; kj++)
						for(int ki=rect.left; ki<rect.right; ki++)
						{
							double sumValue = 0;
							double sumWeight = 0;
							for(int id=0; id<imageIndex.size(); id++)
							{
								int imageId = imageIndex[id];
								
								int ht = vecSIP[id]->nHt;
								int wd = vecSIP[id]->nWd;
								
								int maskSwd = vecSMP[id]->widthStep;

								int cy = kj - vecRect[imageId].top;
								int cx = ki - vecRect[imageId].left;
								
								if( cy<0 || cx<0 || cx>=wd || cy>= ht )
									continue;

								sumValue += (vecSIP[id]->pImage[ cy*wd + cx]) *
									        (unsigned char)(vecSMP[id]->imageData[cy*maskSwd+cx]);
								sumWeight += (unsigned char)(vecSMP[id]->imageData[cy*maskSwd+cx]);								
							}

							if(sumWeight==0)
								continue;
							
							pMosaicImage[kj*owd+ki] += sumValue/sumWeight;
						}
						//SaveToJpg( "d:\\mosaicImage.jpg")
						for(int k=0; k<imageIndex.size(); k++)
						{
							free( vecSIP[k]->pImage );
							cvReleaseImage( &vecSMP[k] );
						}
					}			
				
					for(int k=0; k<imageIndex.size(); k++)
					{
						cvReleaseImage( &vecImagePointer[k] );
						cvReleaseImage( &vecMaskPointer[k] );
					}
				}

				unsigned char* pBuffer = (unsigned char*)malloc(oht*owd);
				for(int ki=0; ki<oht*owd; ki++)
				{
					if(pMosaicImage[ki]>255)
						pMosaicImage[ki]=255;
					if(pMosaicImage[ki]<0)
						pMosaicImage[ki]=0;
					pBuffer[ki] = pMosaicImage[ki];
				}
				poBand = poDataset->GetRasterBand( 3-bandId );
				poBand->RasterIO(GF_Write, 0, 0, owd, oht, pBuffer, owd, oht, GDT_Byte, 0, 0);
				
				free(pBuffer);
				
	}

	GDALClose( (GDALDatasetH) poDataset );	

	return 0;
}

int WeightedLoGBlend( char** filenames, char** masknames, int nFile, 
					 vector<stGeoInfo> geoArray,
					 vector<MyRect> vecRect,
					 double outImageResolution, 
					 double maskResolution,
					 double minx,double maxx,double miny, double maxy,
					 int nPyramidLevel,	
					 double* gainParas,
					 char* outFile )
{
	if(nFile<=1)
		return -1;

	int oht = (maxy-miny) / outImageResolution;
	int owd = (maxx-minx) / outImageResolution;
	printf("oht: %d  owd: %d \n", oht, owd);

	//generate the whole mask
	printf("generate the whole mask....\n");
	int mht = (maxy-miny) / maskResolution;
	int mwd = (maxx-minx) / maskResolution;
	unsigned char* pWholeMask = (unsigned char*)malloc(mht*mwd);
	if(pWholeMask==NULL)
	{
		printf("failed to malloc memory for whole mask .... \n");
		return -1;
	}
	memset(pWholeMask, 0, mht*mwd);
	double maskRatio = outImageResolution / maskResolution;
	for(int i=0; i<nFile; i++)
	{		
		//reading mask
		unsigned char* pMask = NULL;
		int ht,wd;
		ReadRawImage(&pMask, ht, wd, masknames[i]);
		printf("%d %s  %d %d \n", i, masknames[i], ht, wd);

		MyDilate(pMask, ht, wd, 3);

		int l = max(0, min(owd-1, int(vecRect[i].left)))*maskRatio;
		int r = max(0, min(owd-1, int(vecRect[i].right)))*maskRatio;
		int t = max(0, min(oht-1, int(vecRect[i].top)))*maskRatio;
		int b = max(0, min(oht-1, int(vecRect[i].bottom)))*maskRatio;

		l = min(mwd-1, l);
		r = min(mwd-1, r);
		t = min(mht-1, t);
		b = min(mht-1, b);
		
		printf("%d %d %d %d \n", l, r, t, b);

		for(int kj=t; kj<=b; kj++)
			for(int ki=l; ki<=r; ki++)
			{
				int my = max(0, min(ht-1, (kj-t)));
				int mx = max(0, min(wd-1, (ki-l)));

				if(pMask[my*wd + mx]==0)
					continue;

				pWholeMask[kj*mwd + ki] = pMask[my*wd + mx];
			}

		free(pMask);
	}
	MyErode(pWholeMask, mht, mwd, 3);

#ifdef _DEBUG
	SaveToJpg(pWholeMask, mht, mwd, "d:\\wholemask.jpg");
#endif

	//GDALAllRegister();
	//CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	//stGeoInfo geoinfo;
	//GetGeoInformation(filenames[0], geoinfo);
	
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");


	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	int zoneNumber = geoArray[0].zoneNumber; //geoinfo.zoneNumber;
	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );      
	
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}	
	printf("gdal creating file... \n");
	printf("%d %d \n", owd, oht);
	
	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );
	if(poDataset==NULL)
	{
		printf("Failed to create file using gdal ! \n");
		return -1;
	}

    double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = outImageResolution;
	geoTransform[5] = -outImageResolution;
	poDataset->SetGeoTransform(geoTransform);
#ifndef _DEBUG
	poDataset->SetProjection(pszWKT);
#endif	
	
	double dTime = 55;
	double dStep = dTime / (3*nFile);
	int    pi = 0;
	//generate Laplacian pyramid image one by one
	for(int bandId=0; bandId<3; bandId++)
	{
		//prepare for the mosaic image
		PYRAMID_IMAGE_GENERAL<short> mosaicLaplacianPI;
		MallocPyramidGeneral(oht, owd, mosaicLaplacianPI, nPyramidLevel);

		PYRAMID_IMAGE_GENERAL<short> mosaicLaplacianWeight;
		MallocPyramidGeneral(oht, owd, mosaicLaplacianWeight, nPyramidLevel);
	
		double dstMean,dstStd;
		for(int i=0; i<nFile; i++)
		{
			printf("band: %d  %s \n", bandId, filenames[i]);
						
			//reading image 
			IplImage* pImage = cvLoadImage(filenames[i]);
			//resize the image, added by xiedonghai
			double ratio = fabs(geoArray[i].dx) / outImageResolution;
			int rht = ratio * geoArray[i].ht;
			int rwd = ratio * geoArray[i].wd;
			IplImage* pResizeImage = cvCreateImage(
				cvSize(rwd, rht), pImage->depth, pImage->nChannels);
			cvResize(pImage, pResizeImage);
			//cvSaveImage("d:\\resize.jpg", pResizeImage);
			cvReleaseImage(&pImage);
			
			int ht = pResizeImage->height;
			int wd = pResizeImage->width;
			int scanWd = pResizeImage->widthStep;
			//unsigned short* pImageBuffer = (unsigned short*)malloc(ht*wd*sizeof(unsigned short));
			unsigned char* pImageBuffer = (unsigned char*)malloc(ht*wd*sizeof(unsigned char));
			for(int kj=0; kj<ht; kj++)
				for(int ki=0; ki<wd; ki++)
				{
					pImageBuffer[kj*wd+ki] = (unsigned char)(pResizeImage->imageData[kj*scanWd + 3*ki + bandId]) ;
					pImageBuffer[kj*wd+ki] = min( 255, int(pImageBuffer[kj*wd+ki]*gainParas[i]) );
				}
			cvReleaseImage(&pResizeImage);		


			//reading mask
			unsigned char* pMask = NULL;
			int mht1,mwd1;
			ReadRawImage(&pMask, mht1, mwd1, masknames[i]);
			//assert(mht==ht && mwd==wd);
			//resize the mask
			unsigned char* pResizeMask = (unsigned char*)malloc(ht*wd);
			ResizeImage(pMask, mht1, mwd1, pResizeMask, ht, wd);


			//generate gaussian pyramid and smoothed mask
			PYRAMID_IMAGE_GENERAL<unsigned char> pyramidImage; 
			ConstructPyramidWithMask(pImageBuffer, pResizeMask, ht, wd, pyramidImage, nPyramidLevel);
			free(pImageBuffer);
			free(pMask);
			free(pResizeMask);

#ifdef _DEBUG
			//WritePyramidImage(pyramidImage);
#endif

			//generate laplacian pyramid of each image
			PYRAMID_IMAGE_GENERAL<short> laplacianPI; 
			ConstructLaplacianPyramid(pyramidImage, laplacianPI);

#ifdef _DEBUG
			//WritePyramidImage(laplacianPI);
#endif
			double imageRatio = 1;
			for(int li=0; li<laplacianPI.nL; li++)
			{
				int lht = mosaicLaplacianPI.pHt[li];
				int lwd = mosaicLaplacianPI.pWd[li];

				int iht = laplacianPI.pHt[li];
				int iwd = laplacianPI.pWd[li];				
               
				for(int kj=0; kj<iht; kj++)
					for(int ki=0; ki<iwd; ki++)
					{
						int ay = min( lht-1, kj+ int(vecRect[i].top*imageRatio));
						int ax = min( lwd-1, ki+ int(vecRect[i].left*imageRatio));
						
						int index = ay*lwd + ax;

						mosaicLaplacianPI.pLevelImage[li][ index ] +=
							laplacianPI.pLevelImage[li][kj*iwd+ki] * pyramidImage.pMask[li][kj*iwd+ki];


						mosaicLaplacianWeight.pLevelImage[li][index ] += 
											pyramidImage.pMask[li][kj*iwd+ki];
					}
				imageRatio *= 0.5;

				//WritePyramidImage(mosaicLaplacianPI);
			}

			FreePyramidImageGeneral(pyramidImage);
			FreePyramidImageGeneral(laplacianPI);
			
			//WritePyramidImage(mosaicLaplacianPI);

			WriteProgressValueToFile( dStep );
		}

#ifdef _DEBUG
		//WritePyramidImage(mosaicLaplacianPI);
#endif
		
		//calculate the weighted 
		for(int li=0; li<mosaicLaplacianPI.nL; li++)
		{
			int lht = mosaicLaplacianPI.pHt[li];
			int lwd = mosaicLaplacianPI.pWd[li];
		
			for(int kj=0; kj<lht*lwd; kj++)
				{
					if(mosaicLaplacianWeight.pLevelImage[li][kj] == 0)
						continue;
					mosaicLaplacianPI.pLevelImage[li][kj] /= mosaicLaplacianWeight.pLevelImage[li][kj];
				}
		}

#ifdef _DEBUG
		WritePyramidImage(mosaicLaplacianPI);
#endif

		//erode the last level image
		int nL = mosaicLaplacianPI.nL;
		//MyErode(mosaicLaplacianPI.pLevelImage[nL-1], mosaicLaplacianPI.pHt[nL-1], mosaicLaplacianPI.pWd[nL-1], 2);

		/*
		//smooth the last level image, but may brighten the constructed image
		SmoothImage(mosaicLaplacianPI.pLevelImage[nL-1], 
			mosaicLaplacianPI.pHt[nL-1], mosaicLaplacianPI.pWd[nL-1], 
			2, 2);
		*/
		
#ifdef _DEBUG
		WritePyramidImage(mosaicLaplacianPI);
#endif        

		//reconstruction
		printf("ReconstrucFromLaplacian...\n");
		short* mosaicImage = NULL;
		int rht,rwd;
		ReconstrucFromLaplacian(mosaicLaplacianPI, &mosaicImage, &rht, &rwd );
		//double imageToMaskRatio = outImageResolution/maskResolution;
		unsigned char* pMosaicImage = (unsigned char*)malloc(rht*rwd);
		memset(pMosaicImage, 0, rht*rwd);
		
		for(int j=0; j<rht; j++)
		{
			int my = min(mht-1, int(j*maskRatio+0.5));
			for(int i=0; i<rwd; i++)
			{
				int index = j*rwd+i;

				int mx = min(mwd-1, int(i*maskRatio+0.5));

				if( mosaicImage[index]<0 )
					mosaicImage[index] = 0;
				if( mosaicImage[index]>255 )
					mosaicImage[index] = 255;	

				if(pWholeMask[my*mwd+mx]>0)
					pMosaicImage[index] = mosaicImage[index];		
			}
		}
		
		//char outfile[256];
		//sprintf(outfile, "d:\\blend_%d.jpg", bandId);
		//SaveToJpgGeneral(pMosaicImage, rht, rwd, outfile);
		printf("save band buffer... \n");
		poBand = poDataset->GetRasterBand( 3-bandId );
		poBand->RasterIO(GF_Write, 0, 0, owd, oht, pMosaicImage, owd, oht, GDT_Byte, 0, 0);

		free(pMosaicImage);
		free(mosaicImage);

		FreePyramidImageGeneral(mosaicLaplacianPI);
		FreePyramidImageGeneral(mosaicLaplacianWeight);
	}

	GDALClose( (GDALDatasetH) poDataset );

	free(pWholeMask);

	return 0;
}

int CalculateSeqsPyramidLevel(vector<stGeoInfo> geoArray, double outresolution)
{
	int nPyramidLevel = -1;

	double ratio = geoArray[0].dx / outresolution;

	vector<int> imageMedianSize;
	for(int i=0; i<geoArray.size(); i++)
	{
		imageMedianSize.push_back( (geoArray[i].wd+geoArray[i].ht)*0.5*ratio );
	} 

	//sort the median image size
	sort(imageMedianSize.begin(), imageMedianSize.end());
	nPyramidLevel = CalculatePyramidLevel( imageMedianSize[geoArray.size()*0.5] );

	return nPyramidLevel;
}




/* mosaic the "jpeg" files using mutiple band blend, created by xiedonghai, 2015.12.21
*/
int GeoTiffBlend(char** filenames, int nFile, char* outFile, double outResolution)
{

	if(nFile<2)
		return -1;

	WriteProgressValueToFile(0.0);

	//////////////////  1. get the geoinfo from the geotiff files //////////////
	vector<stGeoInfo> geoArray;
	double            minx,maxx,miny,maxy;
	GetImageGeoInfo(filenames, nFile, 
					geoArray, minx, maxx, miny, maxy);

	//get the proper resolution
	double resolution = outResolution; 
	int oht = (maxy-miny) / resolution;
	int owd = (maxx-minx) / resolution;
	double memRatio = CalculateRatioFor64bitOS(oht, owd, 1);
	oht *= memRatio;
	owd *= memRatio;
	outResolution = resolution / memRatio;
	printf("%lf %d %d \n", memRatio, oht, owd);
	printf(" restricted resolution: %lf \n", resolution);
	WriteProgressValueToFile(5.0);
	
	
		
	//////////////////  2. generate the original mask image and save them into files//////////////
	char** maskNames = f2c(nFile, 256);
	for(int i=0; i<nFile; i++)
	{
		//GenerateProductFile(filenames[i], "product", "msk", &(maskNames[i]) );
		strcpy(maskNames[i], filenames[i]);
		int len = strlen(maskNames[i]);
		strcpy(maskNames[i]+len-3, "msk");
		printf("%s \n", maskNames[i]);
	}
	//set the resolution for the mask
	int iht = geoArray[0].ht;
	int iwd = geoArray[0].wd;
	int nMaskSize = max(iht, iwd);
	double maskResolution = fabs(geoArray[0].dx);
	if( nMaskSize>MAX_MASK_SIZE )
	{
		double maskRatio = (double)(MAX_MASK_SIZE) / (double)(nMaskSize);
		maskResolution = maskResolution / maskRatio;
	}
	if(maskResolution<outResolution)
		maskResolution = outResolution;

	vector<MyRect> vecRectMosaic;
	vecRectMosaic.resize(nFile);
	vector<MyRect> vecRectMask;
	vecRectMask.resize(nFile);
	printf("generating original mask... \n");
	for(int i=0; i<nFile; i++)
	{
		printf("%d \n", i);
		unsigned char* pMask = NULL;
		int mht,mwd;
		GenerateImageMask(filenames[i], 
			fabs(geoArray[i].dx), 
			maskResolution,
			&pMask, mht, mwd, 10);	    
        
		//calculate the translation of each image in mask space
		double mosaicRatio =  fabs(geoArray[i].dx) / outResolution;		
		int x = (geoArray[i].left - minx) / outResolution;
		int y = (maxy - geoArray[i].top)  / outResolution;
		vecRectMosaic[i].top  = y;
		vecRectMosaic[i].left = x;
		vecRectMosaic[i].bottom = y + geoArray[i].ht*mosaicRatio + 0.5;
		vecRectMosaic[i].right  = x + geoArray[i].wd*mosaicRatio + 0.5;


		double maskRatio = fabs(geoArray[i].dx) / maskResolution;		
		x = (geoArray[i].left - minx) / maskResolution;
		y = (maxy - geoArray[i].top)  / maskResolution;
		vecRectMask[i].top  = y;
		vecRectMask[i].left = x;
		vecRectMask[i].bottom = y + geoArray[i].ht*maskRatio + 0.5;
		vecRectMask[i].right  = x + geoArray[i].wd*maskRatio + 0.5;
			
		SaveRawImage(pMask, mht, mwd, maskNames[i]);
		free(pMask);
	}
	WriteProgressValueToFile(10.0);
	



	//////////////////////   3. gain compensation ////////////////////////////////////
	double* pg = (double*)malloc(nFile*sizeof(double));
	for(int i=0; i<nFile; i++)
		pg[i] = 1;
	CalculateGain(maskNames, filenames, nFile, vecRectMask, pg);
    WriteProgressValueToFile(10.0);
	



	//////////////////////   4. finding seams ////////////////////////////////////
	SeamlinesFindingVoronoi(maskNames, nFile, vecRectMask);
	//SeamlinesFindingVoronoiOneByOne(maskNames, nFile, vecRect, oht, owd);
	WriteProgressValueToFile(20.0);
    
	
	///////////////////////  5. blend  ////////////////////////////////////
	//MultiBandBlend(filenames, maskNames, nFile, geoArray, vecRect, 
	//	outResolution, minx, maxx, miny, maxy, outFile);
	
	//if(type==0)
	{
		int nPyramidLevel = CalculateSeqsPyramidLevel(geoArray, outResolution);
		WeightedLoGBlend(filenames, maskNames, nFile, geoArray, vecRectMosaic,
			outResolution, maskResolution, minx, maxx, miny, maxy, nPyramidLevel, pg, outFile );
	}
	/*else
	{	
		//direct blending
		vector<double> cp;
		DirectBlend(filenames, nFile, 
			geoArray, 
			minx,maxx,miny,maxy, 
			resolution, 
			maskResolution, 
			pAllMask, maskHt, maskWd, pg, cp, outFile);	
	}*/

	printf("Finished ! \n");

	return 0;
}

int DirectMosaicGeoTiff(char** filenames, int nFile, char* outFile, double outResolution)
{
	if(nFile<1)
		return -1;

	WriteProgressValueToFile(0.0);

	//reading the geo-infomation from the geotiff files
	printf("reading the geo-infomation from the images \n");
	vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);
	for(int i=0; i<nFile; i++)
	{	
		//printf("%s \n", filenames[i]);
		char* postfix;
		GetPostfix(filenames[i], &postfix);
		printf("postfix: %s \n", postfix);

		if( strcmp(postfix, "tif")==0 || strcmp(postfix, "TIF")==0 )
		{
			GetGeoInformation( filenames[i], geoArray[i]);
		}
		else if(strcmp(postfix, "jpeg")==0 || strcmp(postfix, "JPEG")==0)
		{
			char geofile[256];
			strcpy(geofile, filenames[i]);
			strcpy(geofile+strlen(geofile)-4, "geo");
			ReadGeoFile(geofile, geoArray[i]);
		}
	}
		
	double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		//resolution
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);
		//position
		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top  - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

	//set the output mosaic resolution
	double resolution = outResolution; 
	int oht = (maxy-miny) / resolution;
	int owd = (maxx-minx) / resolution;
	double memRatio = CalculateRatioFor64bitOS(oht, owd, 1);
	oht *= memRatio;
	owd *= memRatio;
	resolution = resolution / memRatio;
	printf("%lf %d %d \n", memRatio, oht, owd);
	printf(" restricted resolution: %lf \n", resolution);
			
	//set the resolution for Mask 
	double maskResolution = resolution*4;
	int maskHt = (maxy-miny) / maskResolution;
	int maskWd = (maxx-minx) / maskResolution;
	int nMaskSize = max(maskHt, maskWd);
	if( nMaskSize>MAX_MASK_SIZE )
	{
		double ratio = (double)(nMaskSize) / (double)(MAX_MASK_SIZE) ;
		maskResolution = maskResolution * ratio;
		maskHt = (maxy-miny) / maskResolution;
		maskWd = (maxx-minx) / maskResolution;
	}
		
	//calculate the translation of each image in mask space
	double mosaicRatio =  fabs(geoArray[0].dx) / resolution;
	vector<POINT2> pts;
	pts.resize(nFile);
	
	vector<MyRect> vecRect;
	vecRect.resize(nFile);

	vector<int> imageMedianSize;
	for(int i=0; i<nFile; i++)
	{
		int x = (geoArray[i].left - minx) / maskResolution;
		int y = (maxy - geoArray[i].top)  / maskResolution;
		pts[i].x = x;
		pts[i].y = y;
		vecRect[i].top = y;
		vecRect[i].left = x;

		imageMedianSize.push_back( (geoArray[i].wd+geoArray[i].ht)*0.5*mosaicRatio );
	} 
	
	//sort the median image size
	sort(imageMedianSize.begin(), imageMedianSize.end());
	int nPyramidLevel = CalculatePyramidLevel( imageMedianSize[nFile*0.5] );
	//nPyramidLevel = min(3, nPyramidLevel);
	printf("pyramid level: %d \n", nPyramidLevel);
	
	char** maskNames = f2c(nFile, 256);
	for(int i=0; i<nFile; i++)
	{
		strcpy(maskNames[i], filenames[i]);
		int len = strlen(maskNames[i]);
		strcpy(maskNames[i]+len-3, "msk");
	}
	
	///////////////////  generate mask files and save into the disk //////////////////////
	printf("generate mask files .... \n");
	double dTime = 30;
	double dStep = dTime / nFile;
	vector<MyRect> maskSize;
	maskSize.resize(nFile);
	for(int i=0; i<nFile; i++)
	{
		unsigned char* pMask = NULL;
		int mht,mwd;

		double ratio = fabs(geoArray[i].dx) / maskResolution;
		stGeoInfo geoInfo;
		GetGeoInformation(filenames[i], geoInfo);
		int ht = geoInfo.ht;
		int wd = geoInfo.wd;
		mht = ht*ratio;
		mwd = wd*ratio;
		ReadGeoFileByte(filenames[i], 1, ratio, &pMask, mht, mwd);
		//convert to mask value
		for(int k=0; k<mht*mwd; k++)
		{
			if(pMask[k]>0)
				pMask[k]=255;
		}

		/*
		GenerateImageMask(filenames[i], 
			fabs(geoArray[i].dx), 
			maskResolution,
			&pMask, mht, mwd);*/

		vecRect[i].right  = vecRect[i].left + mwd;
		vecRect[i].bottom = vecRect[i].top  + mht;

		MorphErodeImage(pMask, mht, mwd, 9);

		MyRect rect;
		rect.height = mht;
		rect.width  = mwd;
		maskSize.push_back(rect);

		//SaveToJpg(pMask, mht, mwd, maskfile);	
		SaveRawImage(pMask, mht, mwd, maskNames[i]);
		free(pMask);

		WriteProgressValueToFile(dStep);
	}
	
	
	//gain compensation
	double* pg = (double*)malloc(nFile*sizeof(double));
	for(int i=0; i<nFile; i++)
		pg[i] = 1;
	//CalculateGain(maskNames, filenames, nFile, vecRect, pg);
	//CalculateGainOneByOne(maskNames, filenames, nFile, vecRect, pg);
	
	//color correction
	vector<double> cp;
	//CalculateGainAndVignettParas(maskNames, filenames, nFile, vecRect, cp);
	//CalculateGainAndVignettParasUsingIntensity(maskNames, filenames, nFile, vecRect, cp);
	
	//////////////////////////// finding seam line and generate the mask for mosaic ///////////////
	char maskfile[256];
	CSeamFinderBase* pSeamFinder = new CVoroniSeamFinder();
	unsigned short* pAllMask = (unsigned short*)malloc(maskHt*maskWd*sizeof(unsigned short));
	memset(pAllMask, 0, maskHt*maskWd*sizeof(unsigned short));
	
	printf("Finding Seam line.... \n");
	for(int i=0; i<nFile; i++)
	{
		printf(" %d ", i);

		unsigned char* pMask1 = NULL;
		int mht1,mwd1;
		ReadRawImage(&pMask1, mht1,mwd1, maskNames[i]);
		//GetImageBuffer(maskfile, &pMask1, mht1, mwd1 );
		//SaveBmp("d:\\mask-i.bmp", pMask1, mht1, mwd1);

		iPoint p1;
		p1.x = pts[i].x;
		p1.y = pts[i].y;

		for(int j=i+1; j<nFile; j++)
		{						
			unsigned char* pMask2 = NULL;
			int mht2,mwd2;
			ReadRawImage(&pMask2, mht2, mwd2, maskNames[j]);
			//GetImageBuffer(maskfile, &pMask2, mht2, mwd2 );
			//SaveBmp("d:\\mask-j.bmp", pMask2, mht2, mwd2);

			iPoint p2;
			p2.x = pts[j].x;
			p2.y = pts[j].y;

			//SaveBmp("d:\\mask1-before.bmp", pMask1, mht1, mwd1);
			//SaveBmp("d:\\mask2-before.bmp", pMask2, mht2, mwd2);
			pSeamFinder->FindPair(pMask1, mht1, mwd1, pMask2, mht2, mwd2, p1, p2);
			//SaveBmp("d:\\mask1-after.bmp", pMask1, mht1, mwd1);
			//SaveBmp("d:\\mask2-after.bmp", pMask2, mht2, mwd2);

			//char maskfile[256];
			strcpy(maskfile, filenames[j]);
			int len = strlen(maskfile);
			strcpy(maskfile+len-3, "msk");
			
			//strcat(maskfile, ".msk");
			//SaveToJpg(pMask2, mht2, mwd2, maskfile);
			SaveRawImage(pMask2, mht2, mwd2, maskNames[j]);

			free(pMask2);
		}

		//char maskfile[256];
		//MorphCloseImage(pMask1, mht1, mwd1);
		//MorphDilateImage(pMask1, mht1, mwd1, 3);
		MyDilate(pMask1, mht1, mwd1, 3);
	
		//SaveToJpg(pMask1, mht1, mwd1, maskfile);
		SaveRawImage(pMask1, mht1, mwd1, maskNames[i]);
		
		for(int kj=0; kj<mht1; kj++)
			for(int ki=0; ki<mwd1; ki++)
			{
				int x = min(maskWd-1, (int)(pts[i].x+ki) );
				int y = min(maskHt-1, (int)(pts[i].y+kj) );
				if( pMask1[kj*mwd1+ki] > 200 )
					pAllMask[ y*maskWd + x ] = i+1;
			}
		free(pMask1);
	}
	delete pSeamFinder;	
	printf("\n");
	WriteProgressValueToFile(10.0);

#ifdef _DEBUG
	GdalWriteImageUShort("d:\\allmask.tif", pAllMask, maskHt, maskWd);
#endif	
	
	////////////////////////   blend ////////////////////////////
	{	
		vector<double> cp;
		
		//get the type of data
		GDALDataType ntype = GetDataType(filenames[0]);
		int nByte = 1;
		unsigned char* pByteBuffer   = NULL;
		unsigned short* pUShortBuffer = NULL;
		short* pShortBuffer = NULL;
		unsigned int* pUIntBuffer = NULL;
		int*    pIntBuffer = NULL;
		float*  pFloatBuffer = NULL;
		double* pDoubleBuffer = NULL;

		switch (ntype)
		{
		case 1: //byte
			nByte = 1;
			printf("unsigned char ...\n");
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pByteBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);					
			break;
		case 2: //unsigned short
			printf("unsigned short ...\n");
			nByte = 2;
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pUShortBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);	
			break;
		case 3: //short
			printf("short ...\n");
			nByte = 2;
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pShortBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);	
			break;
		case 4: //unsigned int
			nByte = 4;
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pUIntBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);	
			break;
		case 5: //int 
			nByte = 4;
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pIntBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);	
			break;
		case 6: //float
			nByte = 4;
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pFloatBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);	
			break;
		case 7: //double
			nByte = 8;
			DirectBlendTemplate(filenames, nFile, geoArray, 
				pDoubleBuffer,ntype, nByte,
				minx,maxx,miny,maxy, 
				outResolution, 
				maskResolution, 
				pAllMask, maskHt, maskWd, pg, cp, outFile);	
			break;
		}


		/*
		DirectBlend(filenames, nFile, 
			geoArray, 
			minx,maxx,miny,maxy, 
			outResolution, 
			maskResolution, 
			pAllMask, maskHt, maskWd, pg, cp, outFile);
			*/
	}
	
	free(pg);
	free(pAllMask);

	printf("Finished ! \n");


	return 0;
}

int BlendMosaicGeoTiff(char** filenames, int nFile, char* outFile, double outResolution)
{
	WriteProgressValueToFile(0.0);

	//reading the geo-infomation from the geotiff files
	printf("reading the geo-infomation from the images \n");
	vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);
	for(int i=0; i<nFile; i++)
	{	
		//printf("%s \n", filenames[i]);
		char* postfix;
		GetPostfix(filenames[i], &postfix);
		printf("postfix: %s \n", postfix);

		if( strcmp(postfix, "tif")==0 || strcmp(postfix, "TIF")==0 )
		{
			GetGeoInformation( filenames[i], geoArray[i]);
		}
		else if(strcmp(postfix, "jpeg")==0 || strcmp(postfix, "JPEG")==0)
		{
			char geofile[256];
			strcpy(geofile, filenames[i]);
			strcpy(geofile+strlen(geofile)-4, "geo");
			ReadGeoFile(geofile, geoArray[i]);
		}
	}
		
	double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		//resolution
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);
		//position
		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top  - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

	
	//get the proper resolution
	double resolution = outResolution; 
	int oht = (maxy-miny) / resolution;
	int owd = (maxx-minx) / resolution;
	double memRatio = CalculateRatioFor64bitOS(oht, owd, 1);
	oht *= memRatio;
	owd *= memRatio;
	outResolution = resolution / memRatio;
	printf("%lf %d %d \n", memRatio, oht, owd);
	printf(" restricted resolution: %lf \n", resolution);
	WriteProgressValueToFile(5.0);
		
	//////////////////  2. generate the original mask image and save them into files//////////////
	char** maskNames = f2c(nFile, 256);
	for(int i=0; i<nFile; i++)
	{
		//GenerateProductFile(filenames[i], "product", "msk", &(maskNames[i]) );
		strcpy(maskNames[i], filenames[i]);
		int len = strlen(maskNames[i]);
		strcpy(maskNames[i]+len-3, "msk");
		printf("%s \n", maskNames[i]);
	}
	
	//set the resolution for the mask
	int iht = geoArray[0].ht;
	int iwd = geoArray[0].wd;
	int nMaskSize = max(iht, iwd);
	double maskResolution = fabs(geoArray[0].dx);
	if( nMaskSize>MAX_MASK_SIZE )
	{
		double maskRatio = (double)(MAX_MASK_SIZE) / (double)(nMaskSize);
		maskResolution = maskResolution / maskRatio;
	}
	if(maskResolution<outResolution)
		maskResolution = outResolution;
	
	vector<MyRect> vecRectMosaic;
	vecRectMosaic.resize(nFile);
	vector<MyRect> vecRectMask;
	vecRectMask.resize(nFile);
	printf("generating original mask... \n");
	for(int i=0; i<nFile; i++)
	{
		printf("%d \n", i);
		unsigned char* pMask = NULL;
		int mht,mwd;
		GenerateImageMask(filenames[i], 
			fabs(geoArray[i].dx), 
			maskResolution,
			&pMask, mht, mwd, 10);	    

		//calculate the translation of each image in mask space
		double mosaicRatio =  fabs(geoArray[i].dx) / outResolution;		
		int x = (geoArray[i].left - minx) / outResolution;
		int y = (maxy - geoArray[i].top)  / outResolution;
		vecRectMosaic[i].top  = y;
		vecRectMosaic[i].left = x;
		vecRectMosaic[i].bottom = y + geoArray[i].ht*mosaicRatio + 0.5;
		vecRectMosaic[i].right  = x + geoArray[i].wd*mosaicRatio + 0.5;


		double maskRatio = fabs(geoArray[i].dx) / maskResolution;		
		x = (geoArray[i].left - minx) / maskResolution;
		y = (maxy - geoArray[i].top)  / maskResolution;
		vecRectMask[i].top  = y;
		vecRectMask[i].left = x;
		vecRectMask[i].bottom = y + geoArray[i].ht*maskRatio + 0.5;
		vecRectMask[i].right  = x + geoArray[i].wd*maskRatio + 0.5;

		SaveRawImage(pMask, mht, mwd, maskNames[i]);
		free(pMask);
	}
	WriteProgressValueToFile(10.0);

	//getchar();

	//////////////////////   3. gain compensation ////////////////////////////////////
	double* pg = (double*)malloc(nFile*sizeof(double));
	for(int i=0; i<nFile; i++)
		pg[i] = 1;
	CalculateGain(maskNames, filenames, nFile, vecRectMask, pg);
	WriteProgressValueToFile(10.0);

	
	//////////////////////   4. finding seams ////////////////////////////////////
	SeamlinesFindingVoronoi(maskNames, nFile, vecRectMask);
	//SeamlinesFindingVoronoiOneByOne(maskNames, nFile, vecRect, oht, owd);
	WriteProgressValueToFile(20.0);
	
		
	{
		int nPyramidLevel = CalculateSeqsPyramidLevel(geoArray, outResolution);
		WeightedLoGBlend(filenames, maskNames, nFile, geoArray, vecRectMosaic,
			resolution, maskResolution, minx, maxx, miny, maxy, nPyramidLevel, pg, outFile);
	}

	free(pg);

	printf("Finished ! \n");

	return 0;
}


//mosaic geotiff files
int MosaicGeoTiff(char** filenames, int nFile, char* outFile, double outResolution, int type)
{
	WriteProgressValueToFile(0.0);
	
	////////////////////////   blend ////////////////////////////
	if(type==0) //direct blend
	{		
		DirectMosaicGeoTiff(filenames, nFile, outFile, outResolution);
	}
	else //multi-band 
	{
		BlendMosaicGeoTiff(filenames, nFile, outFile, outResolution);
	}
		
	printf("Finished ! \n");

	return 0;
}

int LoGBlend(char** filenames, int nFile, char* outFile, int nLevel )
{	
	//reading the geo-infomation from the images
	printf("reading the geo-infomation from the images \n");
	vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);
	for(int i=0; i<nFile; i++)
	{	
		printf("%s \n", filenames[i]);
		GetGeoInformation( filenames[i], geoArray[i]);
	}

	double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		//resolution
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);
		//position
		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top  - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);



	//set the output resolution
	double resolution = maxrx * nLevel; //max(maxrx, maxry);
	int oht = (maxy-miny) / resolution;
	int owd = (maxx-minx) / resolution;


	//set the resolution for Mask 
	double maskResolution = resolution*0.2;

	/*
	double memRatio = CalculateRatioFor64bitOS(oht, owd, 1);
	oht *= memRatio;
	owd *= memRatio;
	resolution = resolution / memRatio;
	*/

	//calculate the translation of each image
	vector<POINT2> pts;
	pts.resize(nFile);
	vector<int> imageMedianSize;
	for(int i=0; i<nFile; i++)
	{
		int x = (geoArray[i].left - minx) / resolution;
		int y = (maxy - geoArray[i].top)  / resolution;
		pts[i].x = x;
		pts[i].y = y;
		imageMedianSize.push_back( (geoArray[i].wd+geoArray[i].ht)*0.5 );
	} 
	
	//sort the median image size
	sort(imageMedianSize.begin(), imageMedianSize.end());
	int nPyramidLevel = CalculatePyramidLevel( imageMedianSize[nFile*0.5] );
	printf("pyramid level: %d \n", nPyramidLevel);
	nPyramidLevel = min(6, nPyramidLevel);

	WriteProgressValueToFile(5);

	//generate mask files
	double dTime = 20;
	double dStep = dTime / nFile;
	vector<MyRect> maskSize;
	maskSize.resize(nFile);
	for(int i=0; i<nFile; i++)
	{
		char maskfile[256];
		strcpy(maskfile, filenames[i]);
		strcat(maskfile, ".msk");
		
		unsigned char* pMask = NULL;
		int mht,mwd;
		GenerateImageMask(filenames[i], fabs(geoArray[i].dx), resolution, &pMask, mht, mwd);
		
		MorphErodeImage(pMask, mht, mwd, 9);
		
		MyRect rect;
		rect.height = mht;
		rect.width  = mwd;
		maskSize.push_back(rect);

		//SaveToJpg(pMask, mht, mwd, maskfile);	
		SaveRawImage(pMask, mht, mwd, maskfile);
		free(pMask);

		WriteProgressValueToFile( int(5+i*dStep));
	}
	WriteProgressValueToFile(25);

	char maskfile[256];
	//generate masks for mosaic
	CSeamFinderBase* pSeamFinder = new CVoroniSeamFinder();
	
	unsigned short* pAllMask = (unsigned short*)malloc(oht*owd*sizeof(unsigned short));
	memset(pAllMask, 0, oht*owd*sizeof(unsigned short));
	//unsigned short* pAllMask = (unsigned short*)malloc(oht*owd*sizeof(unsigned short));
	//memset(pAllMask, 0, oht*owd*sizeof(unsigned short));

	printf("Finding Seam line.... \n");
	dTime = 10;
	dStep = dTime / nFile;

	for(int i=0; i<nFile; i++)
	{
		printf(" %d ", i);
		
		//char maskfile[256];
		strcpy(maskfile, filenames[i]);
		strcat(maskfile, ".msk");
		unsigned char* pMask1 = NULL;
		int mht1,mwd1;
		ReadRawImage(&pMask1, mht1,mwd1, maskfile);
		//GetImageBuffer(maskfile, &pMask1, mht1, mwd1 );
		//SaveBmp("d:\\mask-i.bmp", pMask1, mht1, mwd1);
		
		iPoint p1;
		p1.x = pts[i].x;
		p1.y = pts[i].y;

		for(int j=i+1; j<nFile; j++)
		{						
			//char maskfile[256];
			strcpy(maskfile, filenames[j]);
			strcat(maskfile, ".msk");
			unsigned char* pMask2 = NULL;
			int mht2,mwd2;
			ReadRawImage(&pMask2, mht2, mwd2, maskfile);
			//GetImageBuffer(maskfile, &pMask2, mht2, mwd2 );
			//SaveBmp("d:\\mask-j.bmp", pMask2, mht2, mwd2);

			iPoint p2;
			p2.x = pts[j].x;
			p2.y = pts[j].y;

			//SaveBmp("d:\\mask1-before.bmp", pMask1, mht1, mwd1);
			//SaveBmp("d:\\mask2-before.bmp", pMask2, mht2, mwd2);
			pSeamFinder->FindPair(pMask1, mht1, mwd1, pMask2, mht2, mwd2, p1, p2);
			//SaveBmp("d:\\mask1-after.bmp", pMask1, mht1, mwd1);
			//SaveBmp("d:\\mask2-after.bmp", pMask2, mht2, mwd2);

			//char maskfile[256];
			strcpy(maskfile, filenames[j]);
			strcat(maskfile, ".msk");
			//SaveToJpg(pMask2, mht2, mwd2, maskfile);
			SaveRawImage(pMask2, mht2, mwd2, maskfile);

			free(pMask2);
		}
		
		//char maskfile[256];
		//MorphCloseImage(pMask1, mht1, mwd1);
		MorphDilateImage(pMask1, mht1, mwd1, 3);

		strcpy(maskfile, filenames[i]);
		strcat(maskfile, ".msk");
		//SaveToJpg(pMask1, mht1, mwd1, maskfile);
		SaveRawImage(pMask1, mht1, mwd1, maskfile);


		for(int kj=0; kj<mht1; kj++)
			for(int ki=0; ki<mwd1; ki++)
			{
				int x = min(owd-1, (int)(pts[i].x+ki) );
				int y = min(oht-1, (int)(pts[i].y+kj) );
				
				if( pMask1[kj*mwd1+ki] > 200 )
					pAllMask[ y*owd + x ] = i+1;
			}

		free(pMask1);

		WriteProgressValueToFile( int(25+i*dStep));
	}
	delete pSeamFinder;	

	WriteProgressValueToFile( 35 );

	printf("\n");

	//GdalWriteImageUShort("d:\\allmask.tif", pAllMask, oht, owd);

	LaplacianBlend(filenames, nFile, 
		geoArray, 
		//pts, 
		minx,maxx,miny,maxy, 
		resolution, maskResolution, nPyramidLevel, 
		pAllMask, oht, owd, NULL, outFile);

	free(pAllMask);
	printf("Finished ! \n");

	
	return 0;
}