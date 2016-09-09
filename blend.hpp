
#ifndef BLEND_HEAD 
#define BLEND_HEAD


#include "corelib/commondata.h"

#include "geotiff.h"



/* blend based on Laplacian images of Pyramid
*/
_declspec(dllexport) int LoGBlend(char** filenames, int nfile, char* outFile, int nLevel=1);

/*
   type: 0:direct blend, 1-multiband blend 
*/
//_declspec(dllexport) int BlendMosaic(char** filenames, int nfile, char* outFile, 
//	double outResolution, int type=1);
_declspec(dllexport) int MosaicGeoTiff(char** filenames, int nFile, char* outFile, double outResolution, int type);


//
_declspec(dllexport) int GeoTiffBlend(char** filenames, int nFile, char* outFile, 
						double outResolution);


_declspec(dllexport) int SaveToJpg(unsigned char* pBuffer, int ht, int wd, char* filename);


int GenerateImageMask(char* filename, double rawResolution, double outResolution,
					  unsigned char** pMask, int& mht, int& mwd, int fillValue=255);


template<typename T>
int ReadRectValue(T* pDstMask, int dht, int dwd,
				  T* pSrcMask, int sht, int swd,
				  MyRect dstRect,
				  MyRect srcRect
				  )
{
	for(int j=dstRect.top; j<dstRect.bottom; j++)
		for(int i=dstRect.left; i<dstRect.right; i++)
		{
			pDstMask[ (j-dstRect.top)*dwd+(i-dstRect.left)] 
			= pSrcMask[ (j-srcRect.top)*dwd+(i-srcRect.left) ]; 			
		}

		return 0;
}


template<typename T>
int FillMaskValue(T* pAllMask, int oht, int owd,
				  MyRect allRect,
				  MyRect rect,
				  unsigned char* pMask, int mht, int mwd, 
				  int fillvalue)
{
	for(int j=rect.top; j<rect.bottom; j++)
		for(int i=rect.left; i<rect.right; i++)
		{
			if(pMask[ (j-rect.top)*mwd + (i-rect.left)]>0)
				pAllMask[(j-allRect.top)*owd+i-allRect.left] = fillvalue;
		}

	return 0;
}

template<typename T>
int SaveToJpgGeneral(T* pImage, int ht, int wd, char* filename)
{
	int minv = 1000000;
	int maxv = -100000;
	for(int k=0; k<ht*wd; k++)
	{
		if(minv>pImage[k])
			minv = pImage[k];
		if(maxv<pImage[k])
			maxv = pImage[k];
	}

	printf("minv: %d maxv: %d \n", minv, maxv);

	unsigned char* pImageBuffer = (unsigned char*)malloc(ht*wd);
	for(int k=0; k<ht*wd; k++)
	{
		pImageBuffer[k] = (double)(pImage[k] - minv) / (double)(maxv-minv) * 255;
	}

	//SaveToJpg(pImageBuffer, ht, wd, file);
	SaveToJpg(pImageBuffer, ht, wd, filename);
	
	free(pImageBuffer);

	return 0;
}

template<typename T>
int WritePyramidImage( PYRAMID_IMAGE_GENERAL<T>& pyraImage)
{
	for(int i=0; i<pyraImage.nL; i++)
	{
		char file[256];
		sprintf(file, "d:\\pyramid_%d.jpg", i);
		
		int minv = 1000000;
		int maxv = -100000;
		int ht =  pyraImage.pHt[i];
		int wd =  pyraImage.pWd[i];
		for(int k=0; k<ht*wd; k++)
		{
			if(minv>pyraImage.pLevelImage[i][k])
				minv = pyraImage.pLevelImage[i][k];
			if(maxv<pyraImage.pLevelImage[i][k])
				maxv = pyraImage.pLevelImage[i][k];
		}

		unsigned char* pImageBuffer = (unsigned char*)malloc(ht*wd);
		for(int k=0; k<ht*wd; k++)
		{
			pImageBuffer[k] = (double)(pyraImage.pLevelImage[i][k] - minv) / (double)(maxv-minv+0.001) * 255;
		}
		
		SaveToJpg(pImageBuffer, ht, wd, file);
		free(pImageBuffer);

		//char file[256];
		if(pyraImage.pMask != NULL)
		{
			sprintf(file, "d:\\mask_%d.jpg", i);
			SaveToJpg(pyraImage.pMask[i], pyraImage.pHt[i], pyraImage.pWd[i], file);
		}
	}

	return 0;
}

//direct mosaic fro *.tif files
template<typename T>
int DirectBlendTemplate(char** filenames, int nFile,
	vector<stGeoInfo> geoArray, 
	T* pBuffer, GDALDataType ntype, int nByte,
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
	
	stGeoInfo geoinfo;
	GetGeoInformation(filenames[0], geoinfo);
	int nBand = geoinfo.nband;

	int oht = (maxy-miny) / outImageResolution;
	int owd = (maxx-minx) / outImageResolution;

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
	poDataset = poDriver->Create(outFile, owd, oht, nBand, ntype, papszOptions );
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
	double dStep = dTime / (nBand*nFile);
	int    pi = 0;

	T* pMosaicBuffer = (T*)malloc(oht*owd*nByte);
	
	bool isColorCorrection = false;
	if(cp.size()>0)
		isColorCorrection = true;

	double a1,b1,c1;  //exposure 
	double v1,v2;     //vignetting

	double scale = 1; //outImageResolution/maskResolution;

	for(int bandId=0; bandId<nBand; bandId++)
	{
		memset(pMosaicBuffer, 0, oht*owd*nByte);

		for(int i=0; i<nFile; i++)
		{
			scale = fabs(geoArray[i].dx) / maskResolution;

			/*
			if(isColorCorrection)
			{
				a1 = cp[i*9 + 3*bandId ];
				b1 = cp[i*9 + 3*bandId + 1];	
				c1 = cp[i*9 + 3*bandId + 2];
				v1 = cp[cp.size()-2];
				v2 = cp[cp.size()-1];
			}
			*/

			//printf("band: %d  %s \n", bandId, filenames[i]);
			//printf("gain: %lf \n \n ", gainParas[i]);

			//reading file
			GDALDataset *pSrcDataSet = (GDALDataset*)GDALOpen(filenames[i], GA_ReadOnly); 
			if(pSrcDataSet==NULL)
			{
				printf("Open file using gdal failed ! \n");
				return 0;
			}
			GDALRasterBand *pSrcBand  = pSrcDataSet->GetRasterBand(bandId+1);
			int wd = pSrcBand->GetXSize();
			int ht = pSrcBand->GetYSize();
			T* pSrc = (T*)malloc(nByte*wd*ht);
			pSrcBand->RasterIO(GF_Read,0,0,wd,ht,pSrc,wd,ht,ntype,0,0);

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
				//printf("mask point: %d \n", ptIndex.size());

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

					int rawValue = pSrc[iy*wd+ix];

					/*//color correction
					if(isColorCorrection)
					{
						int correctValue = a1*rawValue + b1*rawValue*rawValue + c1 
							+ v1*rawValue*pow(radius,2) + v2*rawValue*pow(radius,4);		
						correctValue = min(255, correctValue);
						pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = correctValue;
					}
					else if (gainParas!=NULL)
					{
						//global instensity correction
						int value = rawValue; 
						value *= gainParas[i];
						pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = min(255, value);
					}
					else*/
					{
						pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = rawValue;
					}					
				}				

				poBand = poDataset->GetRasterBand( bandId+1 );
				poBand->RasterIO(GF_Write, 0, 0, owd, oht, pMosaicBuffer, owd, oht, ntype, 0, 0);

				free(pSrc);
				GDALClose( (GDALDatasetH) pSrcDataSet );

				WriteProgressValueToFile( dStep );		
		}
	}

	free(pMosaicBuffer);

	GDALClose( (GDALDatasetH) poDataset );
	
	return 0;
}



#endif