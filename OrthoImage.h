

#ifndef ORTHO_IMAGE_HEADER
#define ORTHO_IMAGE_HEADER


#include <vector>
using namespace std;

//coredll
#include "defs.h"
#include "geotiff.h"

//cvlib
#include "bundlerio.hpp"

//gdal
#include "gdal_priv.h"



//same as ReadBundlerOutFile(...)
int ReadBundlerOutFileBak(char* filename, vector<stPOS>& camParas, vector<stTrack>& tracks );

void GrdToImg1(double* R, double* T, double f, double k1, double k2, 
			   double gx, double gy, double gz, double* ix, double* iy, int ht, int wd);

void ImgtoGrd1(double* R, double* T, double f, double ix, double iy, double gz, double* gx, double* gy);
double CalculateRatioFor32bitOS(int srcHt, int srcWd, int nByte);
void FindGroundPlane(stAbsPOS& absPosParams, vector<stPOS>& camParas, vector<stTrack>& tracks);
void CalculateGroudArea(vector<stTrack> tracks, 
						double& minx, double& maxx, double& miny, double& maxy, double& meanHeight);

double CalculateResolution(int ht, int wd, double meanHeight, vector<stPOS> camParas );
int    DLL_EXPORT GetZoneNumber(vector<stPOS> camParas);
int    GenerateMosaicImage(int ht, int wd, 
						 char**  filenames,
						 double  rawResolution,
						 double  outResolution, 
						 int  demHt,
						 int  demWd,
						 double  demScale,
						 float* iz,
						 vector<stPOS>& camParas, stAbsPOS absPosParams,
						 double  minx, double maxx, double miny, double maxy,
						 double  meanHeight,
						 char* outFile);
int ReadingPosFromJpegImages(char* listfile, vector<stPOS>& camParas);
int ReadingPosFile(char* posFile, char* imageListFile, vector<stPOS>& camParas);
int AbsPosEstimation(stAbsPOS& absPosParams, vector<stPOS>& camParas, vector<stTrack>& tracks);


void DLL_EXPORT  FitPlane1(double* px, double* py, double* pz, int np, double* R);

int  DLL_EXPORT  MosaicOnBundleWithDEM(char* imageListFile, char* bundleFile, 
									   char* smoothedFile, 
									   char* posFile, 
									   int nIsGeoCorrection, char* outFile, int nLevel);


int  DLL_EXPORT  MosaicOnBundleWithWYOut(char* imageListFile, char* camfile,  char* ptfile, 
										 int nLevel, char* outFile);


DLL_EXPORT int VignettingCorrectFromFile(char* filename, char* outfile);


template<typename T>
int OrthoResampleGeneral(T* pSrc, 
						 int nByte,    //the number of bytes for image pixel
						 GDALDataType nType,						 
						 char*   filename,
						 double  rawResolution,
						 double  outResolution, 
						 int	 demHt,
						 int     demWd,
						 double  demScale,
						 float*  iz,
						 double  minx,
						 double  miny,
						 double  maxx,
						 double  maxy,
						 stPOS&  camParas,
						 double  meanHeight,
						 int     zoneNumber,
						 char*   outfile)
{

	printf("%s \n", filename);

	GDALAllRegister();

	//handle the Chinese characters
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	stGeoInfo ginfo;
	GetGeoInformation(filename, ginfo);
	int ht = ginfo.ht;
	int wd = ginfo.wd;
	int nBand = ginfo.nband;
	printf("src: %d %d \n", ht, wd);

	double minx1,miny1,maxx1,maxy1;
	minx1 = camParas.l;
	maxx1 = camParas.r;
	miny1 = camParas.t;
	maxy1 = camParas.b;

	int oht = (maxy1 - miny1) / outResolution;
	int owd = (maxx1 - minx1) / outResolution;
	printf("ortho: %d %d \n", oht, owd);

	double tResolution    = 1.0 / outResolution;
	double tDemResolution = 1.0 / (rawResolution*demScale);

	//generate the mapping of resampling
	int* srcIndex = (int*)malloc(owd*oht*sizeof(int));
	memset(srcIndex, -1, owd*oht*sizeof(int));
     

	for(int j=0; j<oht; j++)
	{
		for(int i=0; i<owd; i++)
		{
			double gy = j*outResolution + miny1;
			double gx = i*outResolution + minx1;

			//DEM image coordinate
			int demY = (gy-miny) * tDemResolution;
			int demX = (gx-minx) * tDemResolution;				
			demY = max(0, min(demHt-1, demY) );
			demX = max(0, min(demWd-1, demX) );

			//3D point
			double gP[3];
			gP[0] = gx; 
			gP[1] = gy; 
			gP[2] = iz[demY*demWd+demX];
			//gP[2] = meanHeight;

			//3D to 2D projection
			double ix,iy;
			GrdToImg1(camParas.R,  camParas.T,  camParas.f, 
				camParas.k1, camParas.k2, gP[0], gP[1], gP[2], &ix, &iy, ht, wd);
			ix = int(ix);
			iy = int(iy);

			//resample
			if(ix>=0 && ix<wd && iy>=0 && iy<ht)
			{
				//int index1 = (oht-j-1)*owd + i;
				int index2 = (ht-iy-1)*wd + ix;

				srcIndex[j*owd+i] = index2;
			}
		}
	}


	//reading file
	GDALDataset *pSrcDataSet = (GDALDataset*)GDALOpen(filename, GA_ReadOnly); 
	if(pSrcDataSet==NULL)
	{
		printf("Open file using gdal failed ! \n");
		return 0;
	}

	//create the output file to write the buffer
	char **papszOptions = NULL;
	char* drvName="GTIFF";
	GDALDriver *pDstDriver = GetGDALDriverManager()->GetDriverByName(drvName);
	GDALDataset* pDstDataSet = pDstDriver->Create(outfile, owd, oht, nBand, nType, papszOptions);

	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT ); 
	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx1;
	geoTransform[3] = maxy1;
	geoTransform[1] = outResolution;
	geoTransform[5] = -outResolution;
	pDstDataSet->SetGeoTransform(geoTransform);
#ifndef _DEBUG
	pDstDataSet->SetProjection(pszWKT);
#endif

	T* pDst = (T*)malloc(oht*owd*nByte);	

	double gx,gy;	

	printf("ortho resampling.... \n");
	for(int k=0; k<nBand; k++)
	{
		//reading the source buffer of band
		GDALRasterBand *poBand  = pSrcDataSet->GetRasterBand(k+1);
		wd = poBand->GetXSize();
		ht = poBand->GetYSize();
		printf("band:%d %d %d \n", k+1, ht, wd);

		pSrc = (T*)malloc(nByte*wd*ht);
		if(pSrc==NULL)
		{
			printf("malloc failed! \n ");
			break;
		}
		poBand->RasterIO(GF_Read,0,0,wd,ht,pSrc,wd,ht,nType,0,0);
        
		//clear the destiny buffer
		memset(pDst, 0, oht*owd*nByte);

		printf("fill data ...\n");
		for(int j=0; j<oht; j++)
		{
			for(int i=0; i<owd; i++)
			{				
				//resample
				if( srcIndex[j*owd+i]>0 )
				{
					int index1 = (oht-j-1)*owd + i;
					int index2 = srcIndex[j*owd+i];
					//printf("%d ", index2);
					//assert(index2<(wd*ht));
					//if(index2>(wd*ht))
					//	printf("%d ", index2);

					pDst[index1] = pSrc[index2];
				}
			}
		}
		
		printf("saving band buffer ... \n");
		GDALRasterBand *outBand  = pDstDataSet->GetRasterBand(k+1);
		outBand->RasterIO(GF_Write,0,0,owd,oht,pDst,owd,oht,nType,0,0);	
		
		free(pSrc);
	}

	free(pDst);
	free(srcIndex);

	//close the input file
	GDALClose( (GDALDatasetH) pDstDataSet );
	GDALClose( (GDALDatasetH) pSrcDataSet );

	return 0;
}



#endif

