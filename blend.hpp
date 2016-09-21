
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

	printf("mosaic size - ht:%d  wd:%d \n", oht, owd);

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
	//double scale = 1; //outImageResolution/maskResolution;

	for(int bandId=0; bandId<nBand; bandId++)
	{
		memset(pMosaicBuffer, 0, oht*owd*nByte);

		for(int i=0; i<nFile; i++)
		{
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
			int sl = (geoArray[i].left - minx) / outImageResolution;
			int sr = sl + (geoArray[i].wd*fabs(geoArray[i].dx)) / outImageResolution + 0.5;
			int st = (maxy - geoArray[i].top) / outImageResolution;
			int sb = st + (geoArray[i].ht*fabs(geoArray[i].dy)) / outImageResolution + 0.5;

			//collect pixel index 
			vector<iPoint> ptIndex; 
			ptIndex.clear();
			for(int m=st; m<sb; m++)
			//for(int m=0; m<oht; m++)
			{
				for(int n=sl; n<sr; n++)
				//for( int n=0; n<owd; n++ )
				{
					double x = n*outImageResolution;
					double y = m*outImageResolution;

					//mask space
					int mx = x / maskResolution;
					int my = y / maskResolution;
					mx = max(0, min(mx, mwd-1));
					my = max(0, min(my, mht-1));

					//if( pAllMask[my*mwd+mx] == (i+1) )
					if( pAllMask[my*mwd+mx] > 0 )
					{
						iPoint pi;
						pi.x = n;
						pi.y = m;
						ptIndex.push_back(pi);
					}
				}
			}

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

				//double cx = ix-wd*0.5;
				//double cy = iy-ht*0.5;
				//double radius = sqrt(cx*cx+cy*cy)*scale;				

				int rawValue = pSrc[iy*wd+ix];		

				if(rawValue==0)
					continue;
				
				pMosaicBuffer[ ptIndex[m].y*owd + ptIndex[m].x  ] = rawValue;									
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


//
template<typename T>
int LoGBlendGeneral(char** filenames, char** masknames, int nFile, 
					vector<stGeoInfo> geoArray,
					vector<MyRect> vecRect,
					double outImageResolution, 
					double maskResolution,					
					double minx,double maxx,double miny, double maxy,
					int nPyramidLevel,	
					double* gainParas,
					char* outFile,
					unsigned char* pWholeMask, double maskRatio, int mht, int mwd,
					T* pBuffer, GDALDataType ntype, int nByte,
					int nBand,
					int oht, int owd
				    )
{	
	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");
	GDALDriver *poDriver   = NULL;
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
#ifndef _DEBUG
	poDataset->SetProjection(pszWKT);
#endif	
	
	double dTime = 55;
	double dStep = dTime / (nBand*nFile);
	int    pi = 0;
	//generate Laplacian pyramid image one by one
	for(int bandId=0; bandId<nBand; bandId++)
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
			//IplImage* pImage = cvLoadImage(filenames[i]);
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
			double ratio = fabs(geoArray[i].dx) / outImageResolution;
			int rht = ratio * geoArray[i].ht;
			int rwd = ratio * geoArray[i].wd;			

			T* pSrc = (T*)malloc(nByte*rwd*rht);
			pSrcBand->RasterIO(GF_Read,0,0,wd,ht,pSrc,rwd,rht,ntype,0,0);
			GDALClose( (GDALDatasetH) pSrcDataSet );
						
			for(int kj=0; kj<rht*rwd; kj++)
			{
					pSrc[kj] = pSrc[kj]*gainParas[i] ;
			}
			
			//reading mask
			unsigned char* pMask = NULL;
			int mht1,mwd1;
			ReadRawImage(&pMask, mht1, mwd1, masknames[i]);
			//assert(mht==ht && mwd==wd);
			//resize the mask
			unsigned char* pResizeMask = (unsigned char*)malloc(rht*rwd);
			ResizeImage(pMask, mht1, mwd1, pResizeMask, rht, rwd);


			//generate gaussian pyramid and smoothed mask
			PYRAMID_IMAGE_GENERAL<T> pyramidImage; 
			ConstructPyramidWithMask(pSrc, pResizeMask, ht, wd, pyramidImage, nPyramidLevel);
			
			free(pSrc);
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
		T* pMosaicImage = (T*)malloc(rht*rwd*nByte);
		memset(pMosaicImage, 0, rht*rwd*nByte);
		
		for(int j=0; j<rht; j++)
		{
			int my = min(mht-1, int(j*maskRatio+0.5));
			for(int i=0; i<rwd; i++)
			{
				int index = j*rwd+i;

				int mx = min(mwd-1, int(i*maskRatio+0.5));

				if(pWholeMask[my*mwd+mx]>0)
					pMosaicImage[index] = mosaicImage[index];		
			}
		}
		
		//char outfile[256];
		//sprintf(outfile, "d:\\blend_%d.jpg", bandId);
		//SaveToJpgGeneral(pMosaicImage, rht, rwd, outfile);
		printf("save band buffer... \n");
		poBand = poDataset->GetRasterBand( bandId );
		poBand->RasterIO(GF_Write, 0, 0, owd, oht, pMosaicImage, owd, oht, ntype, 0, 0);

		free(pMosaicImage);
		free(mosaicImage);

		FreePyramidImageGeneral(mosaicLaplacianPI);
		FreePyramidImageGeneral(mosaicLaplacianWeight);
	}

	GDALClose( (GDALDatasetH) poDataset );
	
	return 0;
}


#endif