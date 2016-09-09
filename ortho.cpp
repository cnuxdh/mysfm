#include "ortho.hpp"


//corelib
#include "CalcAngle.h"
#include "image.h"

//coredll
#include "geotiff.h"

//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

//gdal
#include "gdal_priv.h"
#include "ogr_spatialref.h"


COrthoUAV::COrthoUAV()
{

}

COrthoUAV::~COrthoUAV()
{

}

int COrthoUAV::GenerateOrthoImage(CameraPara camInfo, char* inputFile, char* outFile)
{   

	if( fabs(camInfo.ax)>3 || fabs(camInfo.ay)>3 )
		return 0;

	//camera position vector
	camInfo.t[0] = camInfo.xs;
	camInfo.t[1] = camInfo.ys;
	camInfo.t[2] = camInfo.zs;
    
	//rotation matrix
	GenerateRMatrix_dpc2_Y(camInfo.ax, camInfo.ay, -camInfo.az, camInfo.R);

	//read image
	IplImage* pImage = cvLoadImage(inputFile, 1);
	int nband = pImage->nChannels;
	int ht = pImage->height;
	int wd = pImage->width;
	int scanwd = pImage->widthStep;
	COLOR* pColor = (COLOR*)malloc(ht*wd*sizeof(COLOR));
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			if(nband==3)
			{
				pColor[j*wd+i].r = (unsigned char)(pImage->imageData[j*scanwd+i*3]);
				pColor[j*wd+i].g = (unsigned char)(pImage->imageData[j*scanwd+i*3+1]);
				pColor[j*wd+i].b = (unsigned char)(pImage->imageData[j*scanwd+i*3+2]);
			}
			if(nband==1)
			{
				pColor[j*wd+i].r = (unsigned char)(pImage->imageData[j*scanwd+i]);
				pColor[j*wd+i].g = (unsigned char)(pImage->imageData[j*scanwd+i]);
				pColor[j*wd+i].b = (unsigned char)(pImage->imageData[j*scanwd+i]);
			}
		}
	cvReleaseImage(&pImage);
	//SaveBmp("d:\\input.bmp", pColor, ht, wd);


	//calculate the geo-transform matrix for GeoTiff file
	double l,r,t,b;
	GenerateScanArea(camInfo.R, camInfo.t, camInfo.focus, camInfo.u0, camInfo.v0, ht, wd, 0, &l, &r, &t, &b);			
    	
	//construct the GeoTiff file 
	double resolution = (r-l)/wd;
    
	//
	double LongCen = camInfo.lon;
	int nZoneNumber = int((LongCen + 180)/6) + 1;
	
	
	//char    *pszWKT ="PROJCS["WGS 84 / UTM zone 50N",GEOGCS["Ellipse_WGS_84",DATUM["unknown",
	//                  SPHEROID["WGS 84",6378137,298.257223563,AUTHORITY["EPSG","7030"]]],PRIMEM["Greenwich",0],UNIT[,0.0174532925199433]],PROJECTION["Transverse_Mercator"],PARAMETER["latitude_of_origin",0],PARAMETER["central_meridian",117],PARAMETER["scale_factor",0.9996],PARAMETER["false_easting",500000],PARAMETER["false_northing",0],UNIT["metre",1,AUTHORITY["EPSG","9001"]],AUTHORITY["EPSG","32650"]]";  
	OGRSpatialReference oSRS;
	//oSRS.SetProjCS("WGS 84 / UTM zone 51N");
	oSRS.SetUTM(nZoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );  
	printf( "%s\n", pszWKT );  
	
	//generate ortho-image
	int oht = (b-t)/resolution;
	int owd = (r-l)/resolution;
	unsigned char* ir = (unsigned char*)malloc(oht*owd);
	unsigned char* ig = (unsigned char*)malloc(oht*owd);
	unsigned char* ib = (unsigned char*)malloc(oht*owd);
	memset(ir, 0, oht*owd);
	memset(ig, 0, oht*owd);
	memset(ib, 0, oht*owd);
	GenerateOrtho( camInfo.R, camInfo.t, camInfo.focus, 
				   camInfo.u0, camInfo.v0, l, r, t, b, resolution, pColor, ht, wd, ir, ig, ib, oht, owd);
 	
	GDALAllRegister();
	GDALDriver* poDriver = NULL;
	GDALDataset *poDataset = NULL;   //GDALÊý¾Ý¼¯
	GDALRasterBand *poBand = NULL;
	char **papszOptions = NULL;
	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTIFF");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	
	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = l;
	geoTransform[3] = b;
	geoTransform[1] = resolution;
	geoTransform[5] = -resolution;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(pszWKT);
	
	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Write, 0, 0, owd, oht, ib, owd, oht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(2);
	poBand->RasterIO(GF_Write, 0, 0, owd, oht, ig, owd, oht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(3);
	poBand->RasterIO(GF_Write, 0, 0, owd, oht, ir, owd, oht, GDT_Byte, 0, 0);

	//close
	GDALClose( (GDALDatasetH) poDataset );

	//GdalWriteImageByteColor("d:\\ortho.jpg", ir, ig, ib, oht, owd);
	GdalWriteJpgCopy("d:\\ortho.jpg", ir, oht, owd);

	free(ir);
	free(ig);
	free(ib);
	free(pColor);

	return 1;
}



//////////////////////////////////////////////////////////////////////////
void  GenerateOrtho(double* R, double* T, 
						 double f, double x0, double y0,
						 double l, double r, double t, double b,
						 double resolution,
						 COLOR* image, int ht, int wd,						 
						 unsigned char* ir, 
						 unsigned char* ig, 
						 unsigned char* ib, 
						 int oht, int owd)
{

	int oj,oi;
	int ix,iy;
	double dx,dy;

	oj = 0;
	oi = 0;
	for(double dj=b-resolution; dj>(t+resolution); dj-=resolution)
	{
		oi = 0;
		for(double di=l+resolution; di<r-resolution; di+=resolution)
		{
			GrdToImg(di, dj, 0, &dx, &dy, R, T, f, x0, y0, ht, wd);

			ix = dx;
			iy = dy;

			//ty =  oht - (int)( (dj-at)*iResolution) - 1;
			//tx =  (int)( (di-al)*iResolution );

			if(ix>=0 && ix<wd && iy>=0 && iy<ht)
			{
				ir[oj*owd + oi] = image[iy*wd+ix].r;
				ig[oj*owd + oi] = image[iy*wd+ix].g;
				ib[oj*owd + oi] = image[iy*wd+ix].b;
			}			
			oi ++;
		}
		oj ++;
	}
}