
#include "fusion.hpp"


//corelib
#include "corelib/commondata.h"
#include "corelib/ImageFunc.h"

//gdal
#include "gdal_priv.h"

//geotiff
#include "geotiff.h"

//opencv
#include "cv.h"
#include "cxcore.h"
#include "highgui.h"





//using RGB to calculate the intensity
double CalculateIntensity(int r, int g, int b)
{
	double ratio = 1.0 / 3.0;
	double intensity = (r  + g  + b ) * ratio;
	return intensity;
}

//
int  CalculateMeanAndSigma(unsigned short* pBuffer, int ht, int wd, double& mean, double& sigma)
{
	mean = 0;
	sigma = 0;

	for(int i=0; i<ht*wd; i++)
	{
		mean += pBuffer[i];
	}
	mean /= (double)(ht*wd);

	for(int i=0; i<ht*wd; i++)
	{
		sigma += ( pBuffer[i] - mean )*( pBuffer[i] - mean );
	}
	sigma = sqrt( sigma / (ht*wd-1) );


	return 0;
}


/*
inputs:
	spectralFile: low resolution multiple spectral file
	panfile:      high resolution panchromatic file
*/
int GF_IHSFusion(char* spectralFile, char* panFile)
{
	//reading MSS image and resized to the pan dimensions
	unsigned short* pBlue = NULL;
	unsigned short* pGreen = NULL;
	unsigned short* pRed = NULL;

	stGeoInfo geoInfo;
	GetGeoInformation(spectralFile, geoInfo);
	
	GDALAllRegister();

	GDALDataset *poDataset = (GDALDataset*)GDALOpen(spectralFile, GA_ReadOnly); 
	GDALRasterBand *poBand  = poDataset->GetRasterBand(1);
	int wd = poBand->GetXSize();
	int ht = poBand->GetYSize();
	int srcWd = wd*4;
	int srcHt = ht*4;
	pBlue = (unsigned short*)malloc(sizeof(unsigned short)*srcWd*srcWd);
	poBand->RasterIO(GF_Read,0,0,wd,ht,pBlue,srcWd,srcHt,GDT_UInt16,0,0);
	
	poBand  = poDataset->GetRasterBand(2);
	pGreen = (unsigned short*)malloc(sizeof(unsigned short)*srcWd*srcWd);
	poBand->RasterIO(GF_Read,0,0,wd,ht,pGreen,srcWd,srcHt,GDT_UInt16,0,0);
	
	poBand  = poDataset->GetRasterBand(3);
	pRed = (unsigned short*)malloc(sizeof(unsigned short)*srcWd*srcWd);
	poBand->RasterIO(GF_Read,0,0,wd,ht,pRed,srcWd,srcHt,GDT_UInt16,0,0);
	
	GDALClose( (GDALDatasetH) poDataset );

	geoInfo.dx = geoInfo.dx / 4;
	geoInfo.dy = geoInfo.dy / 4;
	geoInfo.ht *= 4;
	geoInfo.wd *= 4;	

	//calculate the Intensity
	unsigned short* pIntensity = (unsigned short*)malloc(sizeof(unsigned short)*srcWd*srcWd);
	for(int j=0; j<srcHt; j++)
		for(int i=0; i<srcWd; i++)
		{
			int index = j*srcWd+i;
			pIntensity[index] = CalculateIntensity(pRed[index], pGreen[index], pBlue[index]);
		}
	
	//save intensity
	GdalWriteUShort("d:\\intensity.tif", pIntensity, srcHt, srcWd, geoInfo);


	//motion between intensity and pan 


	//reading PAN image
	unsigned short* pPan;
	ReadGeoFileUShort(panFile, 1, &pPan, srcHt, srcWd);
	GdalWriteUShort("d:\\pan.tif", pPan, srcHt, srcWd, geoInfo);

		
	//histogram match
	double meanP, meanI;
	double sigmaP, sigmaI;
	CalculateMeanAndSigma(pIntensity, srcHt, srcWd, meanI, sigmaI);
	CalculateMeanAndSigma(pPan, srcHt, srcWd, meanP, sigmaP);
	for(int i=0; i<srcHt*srcWd; i++)
	{
		pPan[i] = (sigmaI/sigmaP) * (pPan[i]-meanP) + meanI;
	}

	//save previous image
	GdalWriteColorImageUShort("d:\\MSS-Resize.tif", pRed, pGreen, pBlue, 
		srcHt, srcWd, geoInfo);

	//fusion
	for(int i=0; i<srcHt*srcWd; i++)
	{
		double dif = ( pPan[i]-pIntensity[i] );
		pBlue[i]  += dif; 
		pGreen[i] += dif; 
		pRed[i]   += dif; 
	}

	//save fusion image
	GdalWriteColorImageUShort("d:\\fusion.tif", pRed, pGreen, pBlue, 
		srcHt, srcWd, geoInfo);
	
	free(pBlue);
	free(pGreen);
	free(pRed);
	free(pPan);

	return 0;
}




