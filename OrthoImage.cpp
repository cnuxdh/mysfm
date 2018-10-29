
#include"stdio.h"

#include "OrthoImage.h"

//opencv dll
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>

//gdal dll
#include "gdal_priv.h"
#include "ogr_spatialref.h"
  
//corelib 
#include "Matrix.h"
#include "FitObject.h"
#include "commondata.h"
#include "commonfile.h"
#include "CommonFuncs.h"
#include "LatLong-UTMconversion.h"
#include "ImageFunc.h"

//cvdll
#include "absOri.hpp"
#include "ColorCorrection.hpp"

//coredll
#include "geotiff.h"

//mba dll
#include "mbaExports.h"


int ReadSmoothedPointCloud(char* filename, vector<stTrack>& tracks)
{
	tracks.clear();

	int nPtNum = 0;
	FILE* fp = fopen(filename, "r");
	fscanf(fp, "%d", &nPtNum);
	tracks.resize(nPtNum);
	for(int i=0; i<nPtNum; i++)
	{
		double x,y,z;
		fscanf(fp, "%lf %lf %lf ", &x, &y, &z);
		tracks[i].x = x;
		tracks[i].y = y;
		tracks[i].z = z;
	}
	fclose(fp);

	return 0;
}




/* reading bundler *.out file, including camera pose information, ground point and its projections
*/
//int ReadBundlerOutFileBak(char* filename, vector<stPOS>& camParas, vector<stTrack>& tracks )
int ReadBundlerOutFileBak(char* filename, vector<stPOS>& camParas)
{
	char strLine[256];
	int numCamera, numPt;
	double f,k1,k2;
	double R[9];
	double T[3];
	int i,j;

	FILE* fp =fopen(filename, "r");
	if(fp==0)
	{
		printf(" failed to open %s !\n ", filename );
		return -1;
	}

	fgets(strLine,256, fp);
	printf("%s\n", strLine);

	fscanf(fp, "%d %d", &numCamera, &numPt);
	printf("%d %d\n", numCamera, numPt);

	printf("reading camera information ..... \n");
	stPOS cam;
	for(i=0; i<numCamera; i++)
	{
		

		fscanf(fp, "%lf %lf %lf", &cam.f, &cam.k1, &cam.k2);
		for(j=0; j<9; j++)	
			fscanf(fp, "%lf", &(cam.R[j]) );
		for(j=0; j<3; j++)	
			fscanf(fp, "%lf", &(cam.T[j]) );

		/*
		//calculate the rotation angle
		cam.pitch = atan( cam.R[5]/cam.R[8] )/PI*180; 
		cam.roll  = asin( -cam.R[2] )/PI*180;
		cam.yaw   = atan( cam.R[1]/cam.R[0])/PI*180;
		printf("rotation angle: %lf %lf %lf \n", cam.pitch, cam.roll, cam.yaw);
		*/

		camParas.push_back(cam);
	
		printf("camera: %d  \n", i);
	}

	/*
	//read 3d point
	double rgb[3];
	int nPt;
	int nImageIndex, nPtIndex;
	double ix,iy;

	for( i=0; i<numPt; i++)
	{
		stTrack singleTrack;

		fscanf(fp, "%lf %lf %lf", &singleTrack.x, &singleTrack.y, &singleTrack.z );
		fscanf(fp, "%d %d %d", &singleTrack.r, &singleTrack.g, &singleTrack.b );
		fscanf(fp, "%d", &nPt);

		for(int k=0; k<nPt; k++)
		{
			POINT2 pt;				
			fscanf(fp,"%d %d %lf %lf", &nImageIndex, &nPtIndex, &pt.x, &pt.y);
			singleTrack.imgpt.push_back(pt);
			singleTrack.imgid.push_back(nImageIndex);
			singleTrack.ptid.push_back(nPtIndex);
		}		
		tracks.push_back(singleTrack);
	}
	fclose(fp);
	*/

	return 0;
}


//bundler ground to image model
//model: RX+T
void GrdToImg1(double* R, double* T, double f, double k1, double k2, 
			   double gx, double gy, double gz, double* ix, double* iy, int ht, int wd)
{
	double res[3];
	double P[3];
	double px,py;
	double rp;

	P[0] = gx;
	P[1] = gy;
	P[2] = gz;

	//mult(R, P, res, 3, 3, 1);
	res[0]=R[0]*P[0]+R[1]*P[1]+R[2]*P[2];
	res[1]=R[3]*P[0]+R[4]*P[1]+R[5]*P[2];
	res[2]=R[6]*P[0]+R[7]*P[1]+R[8]*P[2];
	
	res[0] += T[0];
	res[1] += T[1];
	res[2] += T[2];

	double tr = 1.0/res[2];

	px = -res[0]*tr;
	py = -res[1]*tr;			

	double r2 = (px*px+py*py);

	rp = f*(1+k1*r2+k2*r2*r2);
    
	(*ix) = px*rp + (wd>>1);
	(*iy) = py*rp + (ht>>1);
}


//x=RX+T, when the Z is known
void ImgtoGrd1(double* R, double* T, double f, double ix, double iy, double gz, double* gx, double* gy)
{	
	double iR[9];
	double ip[3];
	double iT[3];
	double res[3];

	memcpy(iR, R, sizeof(double)*9);	
	invers_matrix(iR, 3);

	ip[0] = ix;
	ip[1] = iy;
	ip[2] = -f;   
	mult(iR, ip, res, 3, 3, 1);
	mult(iR, T, iT, 3, 3, 1);

	//*gx = -res[0]/res[2]*(iT[2]+gz) - iT[0];
	//*gy = -res[1]/res[2]*(iT[2]+gz) - iT[1];	

	//revised by Xie Donghai, 2015.11.10
	*gx = res[0]/res[2]*(iT[2]+gz) - iT[0];
	*gy = res[1]/res[2]*(iT[2]+gz) - iT[1];	
}


/* 
   for 32bit, the limit for memory is critical, so the dimension of image must be controlled 
   to avoid memory leak, written by xiedonghai
   input:
	srcHt, srcWd: the input ht and wd of image
	nByte: the number of bytes for one pixel
   return value:
	the ratio of input image zoomed out
*/
#define MAX_MEMORY_32BIT 640
double CalculateRatioFor32bitOS(int srcHt, int srcWd, int nByte)
{
	double ratio = 1.0;

	double dht = srcHt;
	double dwd = srcWd;
	double memWanted = (double)(dht*dwd*nByte) / 1024.0 / 1024.0; // MB

	if(memWanted>MAX_MEMORY_32BIT)
	{
		ratio = sqrt( (double)(MAX_MEMORY_32BIT) / memWanted );
	}
	return ratio;
}


void FindGroundPlane(stAbsPOS& absPosParams, vector<stPOS>& camParas, vector<stTrack>& tracks)
{
	int np = tracks.size();
	double* px = (double*)malloc( np*sizeof(double) );
	double* py = (double*)malloc( np*sizeof(double) );
	double* pz = (double*)malloc( np*sizeof(double) );
	for(int i=0; i<tracks.size(); i++)
	{
		px[i] = tracks[i].x;
		py[i] = tracks[i].y;
		pz[i] = tracks[i].z;
	}

	//double planeR[9];
	FitPlane1(px, py, pz, np, absPosParams.R);
	free(px);
	free(py);
	free(pz);
	
	//transform each camera parameters from current coordinate to new coordinate
	double Rg[9];
	memcpy(Rg, absPosParams.R, 9*sizeof(double) );
	invers_matrix(Rg, 3);
    double Tg[3];
	memcpy(Tg, absPosParams.T, 3*sizeof(double) );

	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f == 0)
			continue;

		//new R
		double newR[9];
		mult(camParas[i].R, Rg, newR, 3, 3, 3);

		//new T
		double newT[3];
		mult(newR, Tg, newT, 3, 3, 1);
		for(int k=0; k<3; k++)
			newT[k] = absPosParams.scale*camParas[i].T[k] - newT[k];

		memcpy( camParas[i].R, newR, sizeof(double)*9 );
		memcpy( camParas[i].T, newT, sizeof(double)*3 );

		//rotation angle
		camParas[i].pitch = atan(  camParas[i].R[5]/camParas[i].R[8] )/PI*180; 
		camParas[i].roll  = asin( -camParas[i].R[2] )/PI*180;
		camParas[i].yaw   = atan(  camParas[i].R[1]/camParas[i].R[0])/PI*180;
				
		//convert form RX+T to R( X - (-inv(R)*T) )
		double t1[3];
		double R1[9];
		memcpy(R1, camParas[i].R, 9*sizeof(double) );
		invers_matrix(R1, 3);
		mult(R1, camParas[i].T, t1, 3, 3, 1);
		camParas[i].xs = -t1[0]; 
		camParas[i].ys = -t1[1]; 
		camParas[i].zs = -t1[2];
	}

	//transform each point from free coordinate to absolute coordinate
	//Xg = lamda.R.X + T
	for(int i=0; i<tracks.size(); i++)
	{
		double fP[3];
		fP[0] = tracks[i].x;
		fP[1] = tracks[i].y;
		fP[2] = tracks[i].z;
        
		double tp[3];
		mult(absPosParams.R, fP, tp, 3, 3, 1);

		for(int k=0; k<3; k++)
			fP[k] = absPosParams.scale*tp[k] + absPosParams.T[k];

		tracks[i].x = fP[0];
		tracks[i].y = fP[1];
		tracks[i].z = fP[2];
	}
}

void FindGroundPlane(stAbsPOS& absPosParams, vector<CameraPara>& camParas, vector<TrackInfo>& tracks)
{
	int np = tracks.size();
	double* px = (double*)malloc(np*sizeof(double));
	double* py = (double*)malloc(np*sizeof(double));
	double* pz = (double*)malloc(np*sizeof(double));
	for (int i = 0; i<tracks.size(); i++)
	{
		px[i] = tracks[i].grd.p[0];
		py[i] = tracks[i].grd.p[1];
		pz[i] = tracks[i].grd.p[2];
	}

	//double planeR[9];
	FitPlane1(px, py, pz, np, absPosParams.R);
	free(px);
	free(py);
	free(pz);

	//transform each camera parameters from current coordinate to new coordinate
	double Rg[9];
	memcpy(Rg, absPosParams.R, 9 * sizeof(double));
	invers_matrix(Rg, 3);
	double Tg[3];
	memcpy(Tg, absPosParams.T, 3 * sizeof(double));

	for (int i = 0; i<camParas.size(); i++)
	{
		if (!camParas[i].bIsAddedIntoNet)
			continue;

		//new R
		double newR[9];
		mult(camParas[i].R, Rg, newR, 3, 3, 3);

		//new T
		double newT[3];
		mult(newR, Tg, newT, 3, 3, 1);
		for (int k = 0; k<3; k++)
			newT[k] = absPosParams.scale*camParas[i].T[k] - newT[k];

		memcpy(camParas[i].R, newR, sizeof(double) * 9);
		memcpy(camParas[i].T, newT, sizeof(double) * 3);

		//rotation angle
		camParas[i].pitch = atan(camParas[i].R[5] / camParas[i].R[8]) / PI * 180;
		camParas[i].roll = asin(-camParas[i].R[2]) / PI * 180;
		camParas[i].yaw = atan(camParas[i].R[1] / camParas[i].R[0]) / PI * 180;

		//convert form RX+T to R( X - (-inv(R)*T) )
		double t1[3];
		double R1[9];
		memcpy(R1, camParas[i].R, 9 * sizeof(double));
		invers_matrix(R1, 3);
		mult(R1, camParas[i].T, t1, 3, 3, 1);
		camParas[i].xs = -t1[0];
		camParas[i].ys = -t1[1];
		camParas[i].zs = -t1[2];
	}

	//transform each point from free coordinate to absolute coordinate
	//Xg = lamda.R.X + T
	for (int i = 0; i<tracks.size(); i++)
	{
		double fP[3];
		fP[0] = tracks[i].grd.p[0];
		fP[1] = tracks[i].grd.p[1];
		fP[2] = tracks[i].grd.p[2];

		double tp[3];
		mult(absPosParams.R, fP, tp, 3, 3, 1);

		for (int k = 0; k<3; k++)
			fP[k] = absPosParams.scale*tp[k] + absPosParams.T[k];

		tracks[i].grd.p[0] = fP[0];
		tracks[i].grd.p[1] = fP[1];
		tracks[i].grd.p[2] = fP[2];
	}
}

void CalculateGroudArea(vector<TrackInfo> tracks,
	double& minx, double& maxx, double& miny, double& maxy, double& meanHeight)
{

	double dHeight = 0;
	minx = 100000000;
	maxx = -10000000;
	miny = 100000000;
	maxy = -10000000;

	int np = tracks.size();

	for (int i = 0; i<np; i++)
	{
		double tp[3];
		tp[0] = tracks[i].grd.p[0];
		tp[1] = tracks[i].grd.p[1];
		tp[2] = tracks[i].grd.p[2];

		dHeight += tp[2];
		minx = min(minx, tp[0]);
		maxx = max(maxx, tp[0]);
		miny = min(miny, tp[1]);
		maxy = max(maxy, tp[1]);
	}
	dHeight /= np;
	printf("height: %lf \n", dHeight);
	printf("range: %lf %lf %lf %lf \n", minx, maxx, miny, maxy);

	meanHeight = dHeight;
}

void CalculateGroudArea(vector<stTrack> tracks, 
						double& minx, double& maxx, double& miny, double& maxy, double& meanHeight)
{
	double dHeight = 0;
	minx = 100000000;
	maxx = -10000000;
	miny = 100000000;
	maxy = -10000000;
	
	int np = tracks.size();

	for(int i=0; i<np; i++)
	{
		double tp[3];
		tp[0] = tracks[i].x;
		tp[1] = tracks[i].y;
		tp[2] = tracks[i].z;
        
		dHeight += tp[2];
		minx = min(minx, tp[0]);
		maxx = max(maxx, tp[0]);
		miny = min(miny, tp[1]);
		maxy = max(maxy, tp[1]);
	}
	dHeight /= np;
	printf("height: %lf \n", dHeight);
	printf("range: %lf %lf %lf %lf \n", minx, maxx, miny, maxy);

	meanHeight = dHeight;	
}

//meanHeight: elevation from the sea surface
double CalculateResolution(double meanHeight, vector<CameraPara> camParas)
{
	double sInterval = 0;
	int num = 0;

	for (int i = 0; i<camParas.size(); i++)
	{
		if (!camParas[i].bIsAddedIntoNet)
			continue;

		int ht = camParas[i].rows;
		int wd = camParas[i].cols;

		double height = fabs(meanHeight - camParas[i].zs);
		double grdWd  = (double)(wd) / camParas[i].focalLen * height;
		double interval = grdWd / (double)(wd);
		//printf(" %lf ", interval);
		sInterval += interval;
		num++;
	}
	sInterval /= (double)(num);

	return sInterval;
}

double CalculateResolution(int ht, int wd, double meanHeight, vector<stPOS> camParas )
{
	double sInterval = 0;	
	int num = 0;
	for(int i=0; i<camParas.size(); i++)
	{
		if( !camParas[i].bHavingPOS )
			continue;

		double height = fabs(meanHeight-camParas[i].zs);
		double grdWd = (double)(wd)/camParas[i].f * height;
		double interval = grdWd/(double)(wd);
		printf(" %lf ", interval);
		sInterval += interval;
		num ++;
	}
	sInterval /= (double)( num );
    
	return sInterval;
}


int GetZoneNumber(vector<CameraPara> camParas)
{
	double slon = 0;
	int n = 0;
	for (int i = 0; i<camParas.size(); i++)
	{
		if (camParas[i].bIsAddedIntoNet)
		{
			slon += camParas[i].lon;
			n++;
		}
	}

	slon /= (double)(n);
	int zoneNumber = int((slon + 180) / 6) + 1;

	return zoneNumber;
}

int GetZoneNumber(vector<stPOS> camParas)
{
	double slon = 0;
	int n = 0;
	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f != 0)
		{
			slon += camParas[i].lon;
			n ++;
		}
	}

	slon /= (double)(n);
	int zoneNumber = int((slon + 180)/6) + 1;

	return zoneNumber;
}

/*
int VignettingCorrect(IplImage* pImage)
{
	int ht = pImage->height;
	int wd = pImage->width;
	int scanWd = pImage->widthStep;
	int nChannels = pImage->nChannels;

	double ratio = 1;
	if(wd>64)
		ratio = 64.0 / double(wd);

	//resize the image
	int sht = ht*ratio  + 0.5;
	int swd = wd*ratio  + 0.5;
	IplImage* pSmallImage = cvCreateImage( cvSize(swd, sht), 8, nChannels );
	cvResize(pImage, pSmallImage);

	
	//convert from image to gray
	IplImage* pGrayImage = NULL; 
	if(nChannels==3)
	{
		pGrayImage = cvCreateImage( cvSize(swd, sht), 8, 1 );
		cvCvtColor(pSmallImage, pGrayImage, CV_BGR2GRAY);
	}
	else
	{
		pGrayImage = cvCloneImage(pSmallImage);	
	}
    
	unsigned char* pImageBuffer = (unsigned char*)malloc( sht*swd );
	for(int j=0; j<sht; j++)
		for(int i=0; i<swd; i++)
		{
			pImageBuffer[j*swd+i] 
			= (unsigned char)(pGrayImage->imageData[j*pGrayImage->widthStep+i]);
		}
	
	//vignetting correction
	vector<double> vp;
	VignettingCorrectionUsingRG(pImageBuffer, sht, swd, vp);

	
	int nV = vp.size();
	for(int i=0; i<nV; i++)
	{
		vp[i] = exp(vp[i]);
		//printf("%lf \n", vp[i]);
	}
	

	//int nV = vp.size();
	double maxVr = vp[0];
	for(int i=0; i<nV; i++)
	{
		vp[i] = vp[i]/maxVr;
		printf("%lf ", vp[i]);
	}
	
	//apply the vignetting correction
	int halfHt = ht*0.5;
	int halfWd = wd*0.5;
	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			double cx = i-halfWd;
			double cy = j-halfHt;

			double radius = sqrt(cx*cx+cy*cy) + 0.5;
			radius *= ratio;
            
			//linear interpolation
			double vValue = 1;
			int nR = int(radius);
			if( nR==0 )
			{
				vValue = vp[0];					
			}
			else if(nR<nV)
			{
				double dr = radius - nR;
				vValue = vp[nR-1]*(1-dr) + vp[nR]*dr;
			}
			else
			{
				vValue = vp[nV-1];
			}

			//radius = max(0, min(nV-1, radius) );			
			//double scale = 1.0 / vp[radius];
			double scale = 1.0/vValue;
						
			int r = (unsigned char)(pImage->imageData[j*scanWd+i*3]);
			int g = (unsigned char)(pImage->imageData[j*scanWd+i*3+1]);
			int b = (unsigned char)(pImage->imageData[j*scanWd+i*3+2]);
			r *= scale;
			g *= scale;
			b *= scale;
			
			//pImageColor->imageData[(ht-j-1)*scanWd+i*3]   = min(255,r);
			//pImageColor->imageData[(ht-j-1)*scanWd+i*3+1] = min(255,g);
			//pImageColor->imageData[(ht-j-1)*scanWd+i*3+2] = min(255,b);
			pImage->imageData[(j)*scanWd+i*3]   = min(255,r);
			pImage->imageData[(j)*scanWd+i*3+1] = min(255,g);
			pImage->imageData[(j)*scanWd+i*3+2] = min(255,b);
		}

	free(pImageBuffer);
	cvReleaseImage(&pGrayImage);
	cvReleaseImage(&pSmallImage);

	return 0;
}

int VignettingCorrectFromFile(char* filename, char* outfile)
{
	IplImage* pImage = cvLoadImage(filename);
	VignettingCorrect(pImage);
	cvSaveImage(outfile, pImage);
	cvReleaseImage(&pImage);

	return 0;
}*/


/* generate one orthoimage as jpeg and header file, written by xiedonghai, 2016.1.5
*/
int  GenerateSingleZipOrthoRGB(char*   filename,  // jpg 
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
						 char*   outImageFile, //.jpeg
						 char*   outGeoFile    //.geo
						 )
{
	if(camParas.f==0)
		return -1;
	
	IplImage* pImage = cvLoadImage(filename, 1);
	int ht = pImage->height;
	int wd = pImage->width;

	///////// Single Image Color Correction, added by Xie Donghai, 2015.12.15//////////////////
    //VignettingCorrect(pImage);
	///////////////////////////////////////////////////////////////////////////////////////////

	double bx[4],by[4];
	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	double minx1 = 100000000;
	double maxx1 = -10000000;
	double miny1 = 100000000;
	double maxy1 = -10000000;

	double gz = meanHeight;
	for(int j=0; j<4; j++)
	{
		double ix = bx[j];
		double iy = by[j];
		double gx,gy;
		ImgtoGrd1(camParas.R, camParas.T, camParas.f, ix, iy, gz, &gx, &gy);

		if(minx1>gx) minx1=gx;
		if(maxx1<gx) maxx1=gx;
		if(miny1>gy) miny1=gy;
		if(maxy1<gy) maxy1=gy;
	}
	camParas.l = minx1;
	camParas.r = maxx1;
	camParas.t = miny1;
	camParas.b = maxy1;

	int oht = (maxy1 - miny1) / outResolution;
	int owd = (maxx1 - minx1) / outResolution;

	//the ortho-images of some photos may be very large because of rotation angle, 
	//just discard them, added by Xie Donghai, 2015.11.28
	int maxImageSize = sqrt( double(ht*ht+wd*wd) );
	if( max(oht,owd)>maxImageSize )
	{
		cvReleaseImage(&pImage);
		return -1;
	}

	FILE* fp = fopen(outGeoFile, "w");
	fprintf(fp, "%d \n", zoneNumber);
	fprintf(fp, "%d %d \n", oht, owd);
	fprintf(fp, "%lf %lf %lf \n", minx1, maxy1, outResolution);
	fclose(fp);


	double tResolution = 1.0 / outResolution;
	double tDemResolution = 1.0 / (rawResolution*demScale);

    unsigned char* pOrtho1 = (unsigned char*)malloc(oht*owd);
	unsigned char* pOrtho2 = (unsigned char*)malloc(oht*owd);
	unsigned char* pOrtho3 = (unsigned char*)malloc(oht*owd);
	memset(pOrtho1, 0, oht*owd);
	memset(pOrtho2, 0, oht*owd);
	memset(pOrtho3, 0, oht*owd);

	double gx,gy;	
	for(int j=0; j<oht; j++)
	{
		for(int i=0; i<owd; i++)
		{
			gy = j*outResolution + miny1;
			gx = i*outResolution + minx1;

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
			if(ix>0 && ix<wd && iy>0 && iy<ht)
			{
				int index1 = (oht-j-1)*owd+i;
				int index2 = (ht-iy-1)*pImage->widthStep + ix*3;

				pOrtho1[ index1 ] = pImage->imageData[ index2 ];
				pOrtho2[ index1 ] = pImage->imageData[ index2+1 ];
				pOrtho3[ index1 ] = pImage->imageData[ index2+2 ];		

			}
		}
	}

	cvReleaseImage(&pImage);		    
	GdalWriteImageByteColor(outImageFile, pOrtho3, pOrtho2, pOrtho1, oht, owd);
	free(pOrtho1);
	free(pOrtho2);
	free(pOrtho3);
}

//generate the ground coordinates for each image, 2018.1.5
int GenerateDOMGeoData(string srcFile, double meanHeight, 
	CameraPara camPara, double outResolution, stGeoInfo& geoinfo)
{
	if (!camPara.bIsAddedIntoNet)
		return -1;

	char* filename = const_cast<char*>(srcFile.c_str());

	stGeoInfo ginfo;
	GetGeoInformation(filename, ginfo);
	int ht = ginfo.ht;
	int wd = ginfo.wd;
	int nBand = ginfo.nband;

	double bx[4], by[4];
	bx[0] = -wd / 2;	by[0] = -ht / 2;
	bx[1] = wd / 2;		by[1] = -ht / 2;
	bx[2] = -wd / 2;	by[2] = ht / 2;
	bx[3] = wd / 2;		by[3] = ht / 2;

	double minx1 = 100000000;
	double maxx1 = -10000000;
	double miny1 = 100000000;
	double maxy1 = -10000000;

	double gz = meanHeight;
	for (int j = 0; j<4; j++)
	{
		double ix = bx[j];
		double iy = by[j];
		double gx, gy;
		ImgtoGrd1(camPara.R, camPara.T, camPara.focalLen, ix, iy, gz, &gx, &gy);

		if (minx1>gx) minx1 = gx;
		if (maxx1<gx) maxx1 = gx;
		if (miny1>gy) miny1 = gy;
		if (maxy1<gy) maxy1 = gy;
	}
	camPara.l = minx1;
	camPara.r = maxx1;
	camPara.t = miny1;
	camPara.b = maxy1;
	int oht = (maxy1 - miny1) / outResolution;
	int owd = (maxx1 - minx1) / outResolution;

	//save the geoinfo
	geoinfo.ht = oht;
	geoinfo.wd = owd;
	geoinfo.dx = outResolution;
	geoinfo.dy = -geoinfo.dx;
	geoinfo.left = minx1;
	geoinfo.top  = miny1;
	
	return 0;
}


//generate orthoimage and save into the memory, added by xiedonghai, 2018.1.4
//meanHeight: dem elevation from the sea surface
int GenerateRGBDOM(string srcFile,
	double rawResolution, double outResolution, double demResolution,
	double minx, double miny, double maxx, double maxy,
	CameraPara camPara, double  meanHeight,
	vector<vector<double>> demData,
	stGeoInfo&  geoinfo, vector<vector<unsigned char>>& rgb, int& oht, int& owd
	)
{
	//if (!camPara.bIsAddedIntoNet)
	//	return -1;

	char* filename = const_cast<char*>(srcFile.c_str());

	int demHt = demData.size();
	int demWd = demData[0].size();

	stGeoInfo ginfo;
	GetGeoInformation(filename, ginfo);
	int ht = ginfo.ht;
	int wd = ginfo.wd;
	int nBand = ginfo.nband;
	
	
	double bx[4], by[4];
	bx[0] = -wd / 2;	by[0] = -ht / 2;
	bx[1] = wd / 2;		by[1] = -ht / 2;
	bx[2] = -wd / 2;	by[2] = ht / 2;
	bx[3] = wd / 2;		by[3] = ht / 2;
	double minx1 = 100000000;
	double maxx1 = -10000000;
	double miny1 = 100000000;
	double maxy1 = -10000000;
	double gz = meanHeight;

	printf("meanHeight: %lf \n", gz);

	for (int j = 0; j<4; j++)
	{
		double ix = bx[j];
		double iy = by[j];
		double gx, gy;
		ImgtoGrd1(camPara.R, camPara.T, camPara.focalLen, ix, iy, gz, &gx, &gy);

		if (minx1>gx) minx1 = gx;
		if (maxx1<gx) maxx1 = gx;
		if (miny1>gy) miny1 = gy;
		if (maxy1<gy) maxy1 = gy;
	}
	printf("[GenerateRGBDOM] ground area: %lf %lf %lf %lf \n", minx1, maxx1, miny1, maxy1);
	

	/*
	//determine the ortho area only based on the height, written by xdh, 2018.1.11
	double ratio = fabs(meanHeight - camPara.zs) / camPara.focalLen;
	double gwd = double(wd) * ratio;
	double ght = double(ht) * ratio;
	double minx1 = camPara.xs - gwd*0.5;
	double maxx1 = camPara.xs + gwd*0.5;
	double miny1 = camPara.ys - ght*0.5;
	double maxy1 = camPara.ys + ght*0.5;
	*/

	camPara.l = minx1;
	camPara.r = maxx1;
	camPara.t = miny1;
	camPara.b = maxy1;
	oht = (maxy1 - miny1) / outResolution;
	owd = (maxx1 - minx1) / outResolution;
	printf("ortho image size: %d  %d \n", oht, owd);


	//save the geoinfo
	geoinfo.ht = oht;
	geoinfo.wd = owd;
	geoinfo.dx = outResolution;
	geoinfo.dy = -geoinfo.dx;
	geoinfo.left = minx1;
	geoinfo.top  = maxy1; //revised by xdh, 2018.5.31


	//image buffer for r,g,b channels
	rgb.resize(3);
	for (int j = 0; j < 3; j++)
		rgb[j].resize(oht*owd);

	//the ortho-images of some photos may be very large because of rotation angle, 
	//just discard them, added by Xie Donghai, 2015.11.28
	int maxImageSize = sqrt(double(ht*ht + wd*wd));
	if (max(oht, owd)>maxImageSize)
	{
		printf("[GenerateRGBDOM]: the size of orthoimage is two big! \n");
		return -1;
	}
	
	//double tResolution = 1.0 / outResolution;
	double tDemResolution = 1.0 / demResolution;

	//generate the mapping of resampling
	int* srcIndex = (int*)malloc(owd*oht*sizeof(int));
	memset(srcIndex, -1, sizeof(owd*oht)*sizeof(int));

	for (int j = 0; j<oht; j++)
	{
		for (int i = 0; i<owd; i++)
		{
			double gy = j*outResolution + miny1;
			double gx = i*outResolution + minx1;

			//DEM image coordinate
			int demY = (gy - miny) * tDemResolution;
			int demX = (gx - minx) * tDemResolution;
			demY = max(0, min(demHt - 1, demY));
			demX = max(0, min(demWd - 1, demX));

			//3D point
			double gP[3];
			gP[0] = gx;
			gP[1] = gy;
            gP[2] = meanHeight; //demData[demY][demX];

			//3D to 2D projection
			double ix, iy;
			GrdToImg1(camPara.R, camPara.T, camPara.focalLen,
				camPara.k1, camPara.k2, gP[0], gP[1], gP[2], &ix, &iy, ht, wd);
			ix = int(ix);
			iy = int(iy);

			//resample
			if (ix>0 && ix<wd && iy>0 && iy<ht)
			{
				//int index1 = (oht-j-1)*owd + i;
				int index2 = (ht - iy - 1)*wd + ix;

				srcIndex[j*owd + i] = index2;
			}
		}
	}

	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8", "NO");
	//reading file
	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filename, GA_ReadOnly);
	if (poDataset == NULL)
	{
		printf("Open file using gdal failed ! \n");
		return -1;
	}

	double gx, gy;
	unsigned char* pBuffer = (unsigned char*)malloc(sizeof(unsigned char)*wd*ht);

	for (int k = 0; k<3; k++)
	{
		//reading the buffer of band
		GDALRasterBand *poBand = poDataset->GetRasterBand(k + 1);
		wd = poBand->GetXSize();
		ht = poBand->GetYSize();

		poBand->RasterIO(GF_Read, 0, 0, wd, ht, pBuffer, wd, ht, GDT_Byte, 0, 0);

		//clear the destiny buffer
		for (int j = 0; j<oht; j++)
		{
			for (int i = 0; i<owd; i++)
			{
				//re-sample
				if (srcIndex[j*owd + i]>0)
				{
					int index1 = (oht - j - 1)*owd + i;
					int index2 = srcIndex[j*owd + i];
					//pDstBuffer[index1] = pBuffer[index2];
					rgb[k][index1] = pBuffer[index2];
				}
			}
		}
	}
	
	free(pBuffer);
	free(srcIndex);
	//close the input file
	GDALClose((GDALDatasetH)poDataset);

	return 0;
}


/* generate one orthoimage as jpeg and header file, written by xiedonghai, 2016.1.5
*/
int  GenerateSingleZipOrtho(char*   filename,  // jpg 
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
						 char*   outImageFile, //.jpeg
						 char*   outGeoFile    //.geo
						 )
{
	if(camParas.f==0)
		return -1;

	stGeoInfo ginfo;
	GetGeoInformation(filename, ginfo);
	int ht = ginfo.ht;
	int wd = ginfo.wd;
	int nBand = ginfo.nband;

	///////// Single Image Color Correction, added by Xie Donghai, 2015.12.15//////////////////
    //VignettingCorrect(pImage);
	///////////////////////////////////////////////////////////////////////////////////////////

	
	double bx[4],by[4];
	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	double minx1 = 100000000;
	double maxx1 = -10000000;
	double miny1 = 100000000;
	double maxy1 = -10000000;

	double gz = meanHeight;
	for(int j=0; j<4; j++)
	{
		double ix = bx[j];
		double iy = by[j];
		double gx,gy;
		ImgtoGrd1(camParas.R, camParas.T, camParas.f, ix, iy, gz, &gx, &gy);

		if(minx1>gx) minx1=gx;
		if(maxx1<gx) maxx1=gx;
		if(miny1>gy) miny1=gy;
		if(maxy1<gy) maxy1=gy;
	}
	camParas.l = minx1;
	camParas.r = maxx1;
	camParas.t = miny1;
	camParas.b = maxy1;

	int oht = (maxy1 - miny1) / outResolution;
	int owd = (maxx1 - minx1) / outResolution;

	//the ortho-images of some photos may be very large because of rotation angle, 
	//just discard them, added by Xie Donghai, 2015.11.28
	int maxImageSize = sqrt( double(ht*ht+wd*wd) );
	if( max(oht,owd)>maxImageSize )
	{
		printf("[GenerateSingleZipOrtho]: the size of orthoimage is two big! \n");
		return -1;
	}

	FILE* fp = fopen(outGeoFile, "w");
	fprintf(fp, "%d \n", zoneNumber);
	fprintf(fp, "%d %d \n", oht, owd);
	fprintf(fp, "%lf %lf %lf \n", minx1, maxy1, outResolution);
	fclose(fp);


	double tResolution = 1.0 / outResolution;
	double tDemResolution = 1.0 / (rawResolution*demScale);

	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	//reading file
	GDALDataset *poDataset = (GDALDataset*)GDALOpen(filename, GA_ReadOnly); 
	if(poDataset==NULL)
	{
		printf("Open file using gdal failed ! \n");
		return 0;
	}
    
	//create the MEM to write the buffer
	char **papszOptions = NULL;
	char* drvName="MEM";
	GDALDriver *memDriver;
	memDriver = GetGDALDriverManager()->GetDriverByName(drvName);
	GDALDataset* ds_preview = memDriver->Create("default", owd, oht, nBand, GDT_Byte, papszOptions);
	
	unsigned char* pDstBuffer = (unsigned char*)malloc(oht*owd);	

	//generate the mapping of resampling
	int* srcIndex = (int*)malloc(owd*oht*sizeof(int));
	memset(srcIndex, -1, sizeof(owd*oht)*sizeof(int));

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
			if(ix>0 && ix<wd && iy>0 && iy<ht)
			{
				//int index1 = (oht-j-1)*owd + i;
				int index2 = (ht-iy-1)*wd + ix;
				
				srcIndex[j*owd+i] = index2;
			}
		}
	}
    

	double gx,gy;	
	for(int k=0; k<nBand; k++)
	{
		//reading the buffer of band
		GDALRasterBand *poBand  = poDataset->GetRasterBand(k+1);
		wd = poBand->GetXSize();
		ht = poBand->GetYSize();

		unsigned char* pBuffer=(unsigned char*)malloc(sizeof(unsigned char)*wd*ht);
		poBand->RasterIO(GF_Read,0,0,wd,ht,pBuffer,wd,ht,GDT_Byte,0,0);
		
		//clear the destiny buffer
		memset(pDstBuffer, 0, oht*owd);

		for(int j=0; j<oht; j++)
		{
			for(int i=0; i<owd; i++)
			{				
				//re-sample
				if( srcIndex[j*owd+i]>0 )
				{
					int index1 = (oht-j-1)*owd + i;
					int index2 = srcIndex[j*owd+i];
					pDstBuffer[index1] = pBuffer[index2];
				}
			}
		}
		
		GDALRasterBand *outBand  = ds_preview->GetRasterBand(k+1);
		outBand->RasterIO(GF_Write,0,0,owd,oht,pDstBuffer,owd,oht,GDT_Byte,0,0);	

		free(pBuffer);
	}
	

	//write the jpeg file
	char* jpgDrvName="JPEG";
	GDALDriver *jpgDriver=GetGDALDriverManager()->GetDriverByName(jpgDrvName);
	GDALDataset *dsjpg=jpgDriver->CreateCopy(outImageFile, ds_preview,0,NULL,NULL,NULL);
	GDALClose( (GDALDatasetH) ds_preview );
	GDALClose( (GDALDatasetH) dsjpg );

	free(pDstBuffer);
	free(srcIndex);

	//close the input file
	GDALClose( (GDALDatasetH) poDataset );
}


/* generate one orthoimage as GeoTiff, written by xiedonghai, 2016.9.1
*/
int  GenerateSingleOrthoTiff(char*   filename,      
							 double  rawResolution,
							 double  outResolution, 
							int	    demHt,
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
							char*   outImageFile
						 )
{
	if(camParas.f==0)
		return -1;

	stGeoInfo ginfo;
	GetGeoInformation(filename, ginfo);
	int ht = ginfo.ht;
	int wd = ginfo.wd;

	///////// Single Image Color Correction, added by Xie Donghai, 2015.12.15//////////////////
	//VignettingCorrect(pImage);
	///////////////////////////////////////////////////////////////////////////////////////////

	double bx[4],by[4];
	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	double minx1 = 100000000;
	double maxx1 = -10000000;
	double miny1 = 100000000;
	double maxy1 = -10000000;

	double gz = meanHeight;
	for(int j=0; j<4; j++)
	{
		double ix = bx[j];
		double iy = by[j];
		double gx,gy;
		ImgtoGrd1(camParas.R, camParas.T, camParas.f, ix, iy, gz, &gx, &gy);

		if(minx1>gx) minx1=gx;
		if(maxx1<gx) maxx1=gx;
		if(miny1>gy) miny1=gy;
		if(maxy1<gy) maxy1=gy;
	}
	camParas.l = minx1;
	camParas.r = maxx1;
	camParas.t = miny1;
	camParas.b = maxy1;

	//the ortho-images of some photos may be very large because of rotation angle, 
	//just discard them, added by Xie Donghai, 2015.11.28
	int oht = (maxy1 - miny1) / outResolution;
	int owd = (maxx1 - minx1) / outResolution;
	int maxImageSize = sqrt( double(ht*ht+wd*wd) );
	if( max(oht,owd)>maxImageSize )
	{
		return -1;
	}

	double tResolution = 1.0 / outResolution;
	double tDemResolution = 1.0 / (rawResolution*demScale);

	//get the type of data
	GDALDataType ntype = GetDataType(filename);
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
		OrthoResampleGeneral(pByteBuffer, nByte, ntype, filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);		
		break;
	case 2: //unsigned short
		printf("unsigned short ...\n");
		nByte = 2;
		OrthoResampleGeneral(pUShortBuffer, nByte, ntype,filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);
		break;
	case 3: //short
		printf("short ...\n");
		nByte = 2;
		OrthoResampleGeneral(pShortBuffer, nByte, ntype,filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);
		break;
	case 4: //unsigned int
		nByte = 4;
		OrthoResampleGeneral(pUIntBuffer, nByte, ntype,filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);
		break;
	case 5: //int 
		nByte = 4;
		OrthoResampleGeneral(pIntBuffer, nByte, ntype,filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);
		break;
	case 6: //float
		nByte = 4;
		OrthoResampleGeneral(pFloatBuffer, nByte, ntype,filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);
		break;
	case 7: //double
		nByte = 8;
		OrthoResampleGeneral(pDoubleBuffer, nByte, ntype,filename, 
			rawResolution, outResolution, demHt,demWd, demScale, iz, 
			minx, miny, maxx, maxy, camParas, meanHeight, zoneNumber, outImageFile);
		break;
	}

	return 0;
}


/* generate one orthoimage from input image and DEM
*/
int  GenerateSingleOrtho(char*   filename,
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
						 char*   outFile)
{
	if(camParas.f==0)
		return -1;

	IplImage* pImage = cvLoadImage(filename, 1);
	int ht = pImage->height;
	int wd = pImage->width;


	///////// Single Image Color Correction, added by Xie Donghai, 2015.12.15//////////////////
    //VignettingCorrect(pImage);
	///////////////////////////////////////////////////////////////////////////////////////////


	double bx[4],by[4];
	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	double minx1 = 100000000;
	double maxx1 = -10000000;
	double miny1 = 100000000;
	double maxy1 = -10000000;

	double gz = meanHeight;
	for(int j=0; j<4; j++)
	{
		double ix = bx[j];
		double iy = by[j];
		double gx,gy;
		ImgtoGrd1(camParas.R, camParas.T, camParas.f, ix, iy, gz, &gx, &gy);

		if(minx1>gx) minx1=gx;
		if(maxx1<gx) maxx1=gx;
		if(miny1>gy) miny1=gy;
		if(maxy1<gy) maxy1=gy;
	}
	camParas.l = minx1;
	camParas.r = maxx1;
	camParas.t = miny1;
	camParas.b = maxy1;

	int oht = (maxy1 - miny1) / outResolution;
	int owd = (maxx1 - minx1) / outResolution;
	
	//the ortho-images of some photos may be very large because of rotation angle, 
	//just discard them, added by Xie Donghai, 2015.11.28
	int maxImageSize = sqrt( double(ht*ht+wd*wd) );
	if( max(oht,owd)>maxImageSize )
	{
		cvReleaseImage(&pImage);
		return -1;
	}

	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );  

	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	printf("gdal creating file... \n");

	
	/*
	double ratio = CalculateRatioFor32bitOS(oht, owd, 3);
	outResolution /= ratio;
	oht = (maxy1 - miny1) / outResolution;
	owd = (maxx1 - minx1) / outResolution;
	*/
	
	
	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );	

	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx1;
	geoTransform[3] = maxy1;
	geoTransform[1] = outResolution;
	geoTransform[5] = -outResolution;
	poDataset->SetGeoTransform(geoTransform);
#ifndef _DEBUG
	poDataset->SetProjection(pszWKT);
#endif


	double tResolution = 1.0 / outResolution;
	double tDemResolution = 1.0 / (rawResolution*demScale);
	unsigned char* pOrtho1 = (unsigned char*)malloc(oht*owd);
	unsigned char* pOrtho2 = (unsigned char*)malloc(oht*owd);
	unsigned char* pOrtho3 = (unsigned char*)malloc(oht*owd);
	memset(pOrtho1, 0, oht*owd);
	memset(pOrtho2, 0, oht*owd);
	memset(pOrtho3, 0, oht*owd);


	double gx,gy;
	//mosaic
	//for( gy=camParas.t; gy<camParas.b; gy+=outResolution )
	for(int j=0; j<oht; j++)
	{
		//for( gx=camParas.l; gx<camParas.r; gx+=outResolution )
		for(int i=0; i<owd; i++)
		{
			gy = j*outResolution + miny1;
			gx = i*outResolution + minx1;

			//orthoimage coordinate
			//int j = (gy-miny1) * tResolution;
			//int i = (gx-minx1) * tResolution;				
			//j = max(0, min(oht-1, j));
			//i = max(0, min(owd-1, i));

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
			if(ix>0 && ix<wd && iy>0 && iy<ht)
			{
				int index1 = (oht-j-1)*owd+i;
				int index2 = (ht-iy-1)*pImage->widthStep + ix*3;

				pOrtho1[ index1 ] = pImage->imageData[ index2 ];
				pOrtho2[ index1 ] = pImage->imageData[ index2+1 ];
				pOrtho3[ index1 ] = pImage->imageData[ index2+2 ];
			}				
		}
	}
	cvReleaseImage(&pImage);		

	
	poBand = poDataset->GetRasterBand( 3 );
	poBand->RasterIO(GF_Write, 0, 0, owd, oht, pOrtho1, owd, oht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand( 2 );
	poBand->RasterIO(GF_Write, 0, 0, owd, oht, pOrtho2, owd, oht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand( 1 );
	poBand->RasterIO(GF_Write, 0, 0, owd, oht, pOrtho3, owd, oht, GDT_Byte, 0, 0);
	GDALClose( (GDALDatasetH) poDataset );	
	
	free(pOrtho1);
	free(pOrtho2);
	free(pOrtho3);
}


char* GenerateTifPath(char* imagefile)
{
	char tifFile[256];

	
	char* pdes = strrchr(imagefile, '\\');
	int   index = pdes - imagefile;

	//get the imagepath
	char imagepath[512];
	strcpy(imagepath, imagefile);
	strcpy(imagepath+index, "\0");
	printf("image path: %s \n", imagepath);

	//get the file title
	char title[512];
	strcpy(title, imagefile+index+1);
	printf("title: %s \n", title);

	//generate new file
	sprintf(tifFile, "%s\\product\\%s", imagepath, title);
	


	//strcpy(tifFile, imagefile);

	int len = strlen(tifFile);
	tifFile[len-3] = 't';
	tifFile[len-2] = 'i';
	tifFile[len-1] = 'f';

	return tifFile;
}

int bIsFileExit(char* filepath)
{
	FILE* fp = fopen(filepath, "r");
	if(fp!=NULL)
	{
		fclose(fp);
		return 1;
	}
	else
		return 0;
}


/*
   generate orthoimage, mosaic and fuse, written by xiedonghai, 2015.10.13
*/
int  GenerateMosaicImageWithFusion(char**  filenames,
						 double  rawResolution,
						 double  outResolution, 
						 int  demHt,
						 int  demWd,
						 double  demScale,
						 float* iz,
						 vector<stPOS>& camParas, stAbsPOS absPosParams,
						 double  minx, double maxx, double miny, double maxy,
						 double  meanHeight,
						 char* outFile)
{
	int zoneNumber = GetZoneNumber(camParas);

	
	//for progress bar
	int nValidate = 0;
	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f > 0)
			nValidate++;
	}
	double dstep = 80.0 / double(nValidate);
	//////////////////////////////////////////////////////////////////////////

	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f == 0)
			continue;

		/*
		char* orthofile = (char*)malloc(256);
		GenerateProductFile(filenames[i], "product", "tif", &orthofile);
		printf("%lf \n", camParas[i].f);
		printf("%s \n",  orthofile);
		GenerateSingleOrtho(filenames[i], rawResolution, outResolution,
			demHt, demWd, demScale, iz, minx, miny, maxx, maxy, camParas[i], meanHeight, zoneNumber, orthofile);	
		free(orthofile);
		*/

		char* orthofile = (char*)malloc(256);
		GenerateProductFile(filenames[i], "product", "jpeg", &orthofile);
		char* geofile = (char*)malloc(256);
		GenerateProductFile(filenames[i], "product", "geo", &geofile);
		
		//if( bIsFileExit(orthofile) && bIsFileExit(geofile) )
		//	continue;

		char* postfix = NULL;
		GetPostfix(filenames[i], &postfix);
		
		if( strcmp(postfix, "jpg")==0 || strcmp(postfix, "JPG")==0 )
		{
			//if input format is jpg, then we compress the orthoimage as jpg and save the geodata into the *.geo file
			GenerateSingleZipOrtho(filenames[i], rawResolution, outResolution,
				demHt, demWd, demScale, iz, minx, miny, maxx, maxy, camParas[i], meanHeight, zoneNumber, 
				orthofile, geofile);
		}
		else if( strcmp(postfix, "tif")==0 || strcmp(postfix, "TIF")==0 )
		{
			//if input format is tif, we save the orthoimage as tif file directly
			GenerateProductFile(filenames[i], "product", "tif", &orthofile);
			GenerateSingleOrthoTiff(filenames[i], rawResolution, outResolution,
				demHt, demWd, demScale, iz, minx, miny, maxx, maxy, camParas[i], meanHeight, zoneNumber, 
				orthofile);
		}
		else
		{

		}

		free(orthofile);
		free(geofile);

		WriteProgressValueToFile(dstep);
	}

	return 0;
}


/* generate mosaic based on ortho image 
   input
	ht,wd: original image dimension
*/
int  GenerateMosaicImage(int ht, int wd, 
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
						 char* outFile)
{

	printf("Mosaic without fusion.... \n");

	double resolution = outResolution;
	
	//generating orthoimage
	int oht = (maxy-miny) / resolution ;
	int owd = (maxx-minx) / resolution ;

	//calculate the area of each image
	double bx[4],by[4];
	bx[0] = -wd/2;	by[0] = -ht/2;
	bx[1] = wd/2;	by[1] = -ht/2;
	bx[2] = -wd/2;	by[2] = ht/2;
	bx[3] = wd/2;	by[3] = ht/2;

	double gz = meanHeight;
	for(int i=0; i<camParas.size(); i++)
	{
		double minx1 = 100000000;
		double maxx1 = -10000000;
		double miny1 = 100000000;
		double maxy1 = -10000000;

		if(camParas[i].f==0)
			continue;

		for(int j=0; j<4; j++)
		{
			double ix = bx[j];
			double iy = by[j];
			double gx,gy;
			ImgtoGrd1(camParas[i].R, camParas[i].T, camParas[i].f, ix, iy, gz, &gx, &gy);

			if(minx1>gx) minx1=gx;
			if(maxx1<gx) maxx1=gx;
			if(miny1>gy) miny1=gy;
			if(maxy1<gy) maxy1=gy;
		}
		camParas[i].l = minx1;
		camParas[i].r = maxx1;
		camParas[i].t = miny1;
		camParas[i].b = maxy1;
	}


	printf("Mosaic..... \n");
	printf("out ht & wd: %d %d \n", oht, owd);
	unsigned char* mosaicImg = (unsigned char*)malloc(oht*owd);
	if(mosaicImg==NULL)
	{
		printf(" [GenerateMosaicImage]: malloc for mosaic image failed.... \n");
		return -1;
	}


	GDALAllRegister();
	CPLSetConfigOption("GDAL_FILENAME_IS_UTF8","NO");

	GDALDriver* poDriver   = NULL;
	GDALDataset *poDataset = NULL;   
	GDALRasterBand *poBand = NULL;
	char **papszOptions    = NULL;
	int zoneNumber = GetZoneNumber(camParas);
	OGRSpatialReference oSRS;
	oSRS.SetUTM(zoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );  
    

	//GDALRasterBand *poBand;
	poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	if( poDriver == NULL )
	{
		printf("Gdal Open file failed ! \n");
		return 0;
	}
	printf("gdal creating file... \n");
	poDataset = poDriver->Create(outFile, owd, oht, 3, GDT_Byte, papszOptions );	
    	
	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = resolution;
	geoTransform[5] = -resolution;
	poDataset->SetGeoTransform(geoTransform);
	poDataset->SetProjection(pszWKT);
    

	double tResolution = 1.0 / resolution;
	double tDemResolution = 1.0 / (rawResolution*demScale);

	double dstep = 80.0 / (double)camParas.size();

	for(int nBandId=0; nBandId<3; nBandId++)
	{
		for(int k=0; k<camParas.size(); k++)
		{	
			WriteProgressValueToFile(dstep);

			printf(" Image: %d  %s \n", k, filenames[k]);

			if(camParas[k].f==0)
				continue;
			
			IplImage* pImage = cvLoadImage(filenames[k], 1);
			double gx,gy;
			
			//mosaic
			for( gy=camParas[k].t; gy<camParas[k].b; gy+=resolution )
			{
				for( gx=camParas[k].l; gx<camParas[k].r; gx+=resolution )
				{
					int j = (gy-miny) * tResolution;
					int i = (gx-minx) * tResolution;				
					j = max(0, min(oht-1, j));
					i = max(0, min(owd-1, i));

					int demY = (gy-miny) * tDemResolution;
					int demX = (gx-minx) * tDemResolution;				

					demY = max(0, min(demHt-1, demY) );
					demX = max(0, min(demWd-1, demX) );

					double gP[3];
					gP[0] = gx; 
					gP[1] = gy; 
					gP[2] = iz[demY*demWd+demX];
										
					double ix,iy;
					GrdToImg1(camParas[k].R,  camParas[k].T,  camParas[k].f, 
						camParas[k].k1, camParas[k].k2, gP[0], gP[1], gP[2], &ix, &iy, ht, wd);

					ix = int(ix);
					iy = int(iy);

					if(ix>0 && ix<wd && iy>0 && iy<ht)
					{
						//CvScalar color = cvGet2D(pImage, ht-iy-1, ix);
						//mosaicImg[ j*owd+i ] = color.val[nBandId];
						mosaicImg[ (oht-j-1)*owd+i ] = pImage->imageData[ (int)( (ht-iy-1)*pImage->widthStep + ix*3 + nBandId) ];
					}				
				}
			}

			cvReleaseImage(&pImage);
		}

		poBand = poDataset->GetRasterBand( 3-nBandId );
		poBand->RasterIO(GF_Write, 0, 0, owd, oht, mosaicImg, owd, oht, GDT_Byte, 0, 0);
	}
	

	GDALClose( (GDALDatasetH) poDataset );

	free(mosaicImg);

	//WriteProgressValueToFile()

	return 0;
}

/* reading Geo Tag information from jpeg image files
   return the number of images containing geo tag
*/
/*
int ReadingPosFromJpegImages(char* listfile, vector<stPOS>& camParas)
{
	int nImageofGeotags = 0;
	
	int nImage = GetFileRows(listfile);
	char imagefile[256];
	int  ti;
	int    ht,wd;
	double focalLen;
	double ccdWidth;
	double lat=0,lon=0,alt=0;
	int    GpsInfoPresent;

	FILE* fp =fopen(listfile, "r");
	if(fp==NULL)
		return -1;

	for(int i=0; i<nImage; i++)
	{

		fscanf(fp, "%s %d %lf ", imagefile, &ti, &focalLen);
		
		char* postfix;
		GetPostfix(imagefile, &postfix);

		GpsInfoPresent = 0;
		if( strcmp(postfix, "jpg")==0 || strcmp(postfix, "JPG")==0 )
		{
			GetInfoFromJpegExif(imagefile, &ht, &wd, 
				&focalLen, &ccdWidth, 
				&lat, &lon, &alt, &GpsInfoPresent);
		}
		
		if(GpsInfoPresent)
		{
			camParas[i].lon = lon;
			camParas[i].lat = lat;
			camParas[i].altitude = alt;
			nImageofGeotags++;
		}
		else
		{
			camParas[i].lon = 0;
			camParas[i].lat = 0;
			camParas[i].altitude = -100000;
			//camParas[i].f = 0;
		}
	}
	fclose(fp);

	return nImageofGeotags;
}
*/


int ReadingPosFile(char* posFile, char* imageListFile, vector<stPOS>& camParas)
{
	vector<stPOS> allPosData;

	int nImage = GetFileRows(posFile) - 1;

	if(nImage<1)
	{
		printf("pos file is not valid! \n");
		return -1;
	}

	char** titles = f2c(nImage, 256);

	char sline[512];
	FILE* fp = fopen(posFile, "r");
	if(fp==NULL)
	{
		printf("error to load pos file! \n");
		return -1;
	}
	fgets(sline, 512, fp);
	for(int i=0; i<nImage; i++)
	{
		stPOS cam;
		fscanf(fp, "%s %lf %lf %lf %lf %lf %lf ", titles[i], 
			&(cam.lat), &(cam.lon), &(cam.altitude),
			&(cam.yaw), &(cam.pitch), &(cam.roll) );
		allPosData.push_back(cam);
	}
	fclose(fp);
	
	int nBundleImage = GetFileRows(imageListFile);
	fp = fopen(imageListFile, "r");
	for(int i=0; i<nBundleImage; i++)
	{
		char filename[256];
		int tn,focalLen;
		fscanf(fp, "%s %d %lf", filename, &tn, &focalLen);

		for(int k=0; k<nImage; k++)
		{
			char* pos = strstr(filename, titles[k]);
			if(pos!=NULL)
			{
				camParas[i].lon = allPosData[k].lon;
				camParas[i].lat = allPosData[k].lat;
				camParas[i].altitude = allPosData[k].altitude;
			}
		}
	}
	fclose(fp);

	return 0;
}


/* absolute pose estimation, written by xiedonghai, 2015.10.2
   
*/
int AbsPosEstimation(stAbsPOS& absPosParams, vector<stPOS>& camParas, vector<stTrack>& tracks)
{		
	//collect valid cameras
	vector<int> validCameraIndex;
	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f != 0)
		{
			validCameraIndex.push_back(i);
		}
	}

	//compute the position of each camera 
	// form RX+T to R( X - (-inv(R)*T) )
	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];		
		
		double t1[3];
		double R1[9];
		memcpy(R1, camParas[ci].R, 9*sizeof(double) );
		invers_matrix(R1, 3);
		mult(R1, camParas[ci].T, t1, 3, 3, 1);
		camParas[ci].xs = -t1[0]; 
		camParas[ci].ys = -t1[1]; 
		camParas[ci].zs = -t1[2];
	}

	//convert from lon,lat to ground coordinate
	int zoneNumber =GetZoneNumber(camParas);
	printf("zonenumber: %d \n", zoneNumber);
	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];
		double lat = camParas[ci].lat;
		double lon = camParas[ci].lon;
		double gx,gy;
		LLtoUTM(23, lat, lon, gy, gx, zoneNumber);
		camParas[ci].gx = gx;
		camParas[ci].gy = gy;
		camParas[ci].gz = camParas[ci].altitude;
		printf("%lf %lf %lf \n", gx, gy, camParas[ci].altitude);
	}

	//calculate the step
	int numCtrlPt = 2;
	int nstep = validCameraIndex.size() / (numCtrlPt+1);
	vector<int> ctrlPtIndex;
	ctrlPtIndex.resize(numCtrlPt+2);
	ctrlPtIndex[0] = validCameraIndex[0];
	ctrlPtIndex[1] = validCameraIndex[validCameraIndex.size()-1];
	for(int i=0; i<numCtrlPt; i++)
	{
		ctrlPtIndex[i+2] = validCameraIndex[nstep*(i+1)];
	}
	//validCameraIndex = ctrlPtIndex;


	//calculate the scale
	double sumScale = 0;
	int    ns = 0;
	for(int j=0; j<ctrlPtIndex.size(); j++)
	{
		int cj = ctrlPtIndex[j];

		double xj = camParas[cj].xs;
		double yj = camParas[cj].ys;
		double zj = camParas[cj].zs;
		double gxj = camParas[cj].gx;
		double gyj = camParas[cj].gy;
		double gzj = camParas[cj].gz;

		for(int i=j+1; i<ctrlPtIndex.size(); i++)
		{
			int ci = ctrlPtIndex[i];

			double xi = camParas[ci].xs;
			double yi = camParas[ci].ys;
			double zi = camParas[ci].zs;
			double gxi = camParas[ci].gx;
			double gyi = camParas[ci].gy;
			double gzi = camParas[ci].gz;
		
			double gLen = sqrt(  (gxj-gxi)*(gxj-gxi) + (gyj-gyi)*(gyj-gyi) + (gzj-gzi)*(gzj-gzi)  );
			double fLen = sqrt(  (xj-xi)*(xj-xi) + (yj-yi)*(yj-yi) + (zj-zi)*(zj-zi)  );

			sumScale += gLen / fLen;

			ns ++;
		}
	}
	absPosParams.scale = sumScale / (double)(ns);

	//calculate the R using LS optimization
	double A[3][9];
	double At[9][3];
	double AtA[9][9];
	double L[3];
	double AtL[9];
	double sAtA[9][9];
	double sAtL[9];
	memset(*sAtA, 0, 81*sizeof(double));
	memset(sAtL, 0, 9*sizeof(double));

	for(int j=0; j<ctrlPtIndex.size(); j++)
	{
		int cj = ctrlPtIndex[j];

		double xj = camParas[cj].xs;
		double yj = camParas[cj].ys;
		double zj = camParas[cj].zs;
		double gxj = camParas[cj].gx;
		double gyj = camParas[cj].gy;
		double gzj = camParas[cj].gz;

		for(int i=j+1; i<ctrlPtIndex.size(); i++)
		{
			int ci = ctrlPtIndex[i];

			double xi = camParas[ci].xs;
			double yi = camParas[ci].ys;
			double zi = camParas[ci].zs;
			double gxi = camParas[ci].gx;
			double gyi = camParas[ci].gy;
			double gzi = camParas[ci].gz;

			double dx = absPosParams.scale*(xj-xi);
			double dy = absPosParams.scale*(yj-yi);
			double dz = absPosParams.scale*(zj-zi);
			double dgx = gxj-gxi;
			double dgy = gyj-gyi;
			double dgz = gzj-gzi;

			memset(*A, 0, 27*sizeof(double));
			memset(*At, 0, 27*sizeof(double));
			memset(*AtA, 0, 81*sizeof(double));
			memset(L, 0, 3*sizeof(double));
			memset(AtL, 0, 9*sizeof(double));

			A[0][0] = dx; A[0][1] = dy; A[0][2] = dz;
			A[1][3] = dx; A[1][4] = dy; A[1][5] = dz;
			A[2][6] = dx; A[2][7] = dy; A[2][8] = dz;
			L[0] = dgx;   L[1] = dgy;   L[2] = dgz;
			
			transpose(*A, *At, 3, 9);
			mult(*At, *A, *AtA, 9, 3, 9);
			mult(*At, L, AtL, 9, 3, 1);
	
			for(int m=0; m<9; m++)
			{
				sAtL[m] += AtL[m];
				for(int n=0; n<9; n++)
				{
					sAtA[m][n] += AtA[m][n];
				}
			}
		}
	}
	invers_matrix(*sAtA, 9);
	mult(*sAtA, sAtL, absPosParams.R, 9, 9, 1);

	//calculate the T
	double sT[3];
	memset(sT, 0, sizeof(double)*3);

	for(int j=0; j<ctrlPtIndex.size(); j++)
	{
		int cj = ctrlPtIndex[j];
	
		double gp[3];
		gp[0] = camParas[cj].gx;
		gp[1] = camParas[cj].gy;
		gp[2] = camParas[cj].gz;
		
		double fp[3];
		fp[0] = camParas[cj].xs;
		fp[1] = camParas[cj].ys;
		fp[2] = camParas[cj].zs;
	
		double frp[3];
		mult(absPosParams.R, fp, frp, 3, 3, 1);

		for(int k=0; k<3; k++)
			sT[k] +=  ( gp[k] - frp[k]*absPosParams.scale);
	}

	for(int k=0; k<3; k++)
		absPosParams.T[k] = sT[k] / (double)(ctrlPtIndex.size());
	
	//transform each camera parameters from current coordinate to new coordinate
	double Rg[9];
	memcpy(Rg, absPosParams.R, 9*sizeof(double) );
	invers_matrix(Rg, 3);
	double Tg[3];
	memcpy(Tg, absPosParams.T, 3*sizeof(double) );

	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];

		//new R
		double newR[9];
		mult(camParas[ci].R, Rg, newR, 3, 3, 3);

		//new T
		double newT[3];
		mult(newR, Tg, newT, 3, 3, 1);
		for(int k=0; k<3; k++)
			newT[k] = absPosParams.scale*camParas[ci].T[k] - newT[k];

		memcpy( camParas[ci].R, newR, sizeof(double)*9 );
		memcpy( camParas[ci].T, newT, sizeof(double)*3 );
		
		//convert form RX+T to R( X - (-inv(R)*T) )
		double t1[3];
		double R1[9];
		memcpy(R1, camParas[ci].R, 9*sizeof(double) );
		invers_matrix(R1, 3);
		mult(R1, camParas[ci].T, t1, 3, 3, 1);
		camParas[ci].xs = -t1[0]; 
		camParas[ci].ys = -t1[1]; 
		camParas[ci].zs = -t1[2];
	}

	//transform all tracks 
	for(int i=0; i<tracks.size(); i++)
	{
		double fP[3];
		fP[0] = tracks[i].x;
		fP[1] = tracks[i].y;
		fP[2] = tracks[i].z;

		double tp[3];
		mult(absPosParams.R, fP, tp, 3, 3, 1);

		for(int k=0; k<3; k++)
			fP[k] = absPosParams.scale*tp[k] + absPosParams.T[k];

		tracks[i].x = fP[0];
		tracks[i].y = fP[1];
		tracks[i].z = fP[2];
	}
	return 0;
}

//
////////////////////////////////////////////////////////////////////////////
//int  ReadCamFile(char* camfile, vector<stPOS>& camParas)
//{
//	FILE* fp = fopen(camfile,"rt");
//	if(fp==NULL)
//	{
//		printf("error to load camera point file ! \n");
//		return -1;
//	}
//
//	int j=0;
//	while( !feof(fp) )
//	{
//		stPOS cam;
//		int nImageIndex = 0;
//		if(fscanf(fp, "%d", &nImageIndex)==EOF) break;
//		cam.nImageIndex = nImageIndex;		// point on image index, image on j
//		//m_rev[nImageIndex] = j;
//		fscanf(fp, "%lf %lf %lf", &cam.f, &cam.k1, &cam.k2);
//		for(int i=0; i<9; i++)	fscanf(fp, "%lf", &(cam.R[i]));
//		for(int i=0; i<3; i++)	fscanf(fp, "%lf", &(cam.T[i]));
//		j++;
//		camParas.push_back(cam);
//	}
//	fclose(fp);
//
//	//sort according the image index
//	for(int j=0; j<camParas.size(); j++)
//		for(int i=j+1; i<camParas.size(); i++)
//		{
//			if(camParas[j].nImageIndex>camParas[i].nImageIndex)
//			{
//				stPOS tp;
//				tp = camParas[j];
//				camParas[j] = camParas[i];
//				camParas[i] = tp;
//			}
//		}
//
//	return 0;
//}
//
//int  ReadPtFile(char* ptfile, vector<stTrack>& tracks)
//{
//	int num_pts;
//	FILE* fp = fopen(ptfile,"rt");
//	if(fp==NULL)
//	{
//		printf("error to read point file ! \n");
//		return -1;
//	}
//
//	fscanf(fp, "%d\n",&num_pts);
//
//	for (int i=0;i<num_pts;i++)
//	{
//		stTrack tp;
//		int ptnum;
//
//		if(fscanf(fp, "%lf %lf %lf %d", &tp.x, &tp.y, &tp.z, &ptnum )==EOF)
//		{
//			num_pts = i;
//			printf("points num: %d\n", num_pts);
//			break;
//		}
//
//		if (ptnum<=0)
//		{
//			continue;
//		}
//		double t1,t2;
//		int    ptindex;
//		for (int j=0;j<ptnum;j++)
//		{
//			fscanf(fp, " %d %lf %lf", &ptindex, &t1, &t2);
//		}
//
//		tracks.push_back(tp);
//	}
//	fclose(fp);
//
//
//	return 0;
//}
//

//
//int  MosaicOnBundleWithWYOut(char* imageListFile, char* camfile,  char* ptfile, 
//							 int nLevel, char* outFile)
//{	
//	//reading camera's parameters
//	vector<stPOS> camParas;
//    if( ReadCamFile(camfile, camParas) < 0 )
//	{
//		printf("Loading cam file failed! \n");
//		return -1;
//	}
//	
//	//reading tracks
//	vector<stTrack> tracks;
//	ReadPtFile(ptfile, tracks);
//	
//	//reading image list file path
//	double focus;
//	int    tn;
//	char   tline[256];
//	int    nImageIndex;
//	int    wd,ht;
//	double scale, td1,td2;
//	int nImage = ( GetFileRows(imageListFile)-1) / 2;
//	if(nImage<1) return -1;
//	char** filenamesAll = f2c(nImage, 256);
//	FILE* fp = fopen(imageListFile, "r");
//	if(fp==NULL)
//	{
//		printf("error reading image list file! \n");
//		return -1;
//	}
//	fgets(tline, 256, fp);
//	for(int i=0; i<nImage; i++)
//	{
//		fscanf(fp, "%d %d %d  %lf %lf %lf ", &nImageIndex, &wd, &ht, &scale, &td1, &td2);
//		fscanf(fp, "\n%[^\n]", filenamesAll[i]);
//	}
//	fclose(fp);
//
//	//generate image list
//	int nImageValid = camParas.size();
//	char** filenames = f2c(nImageValid, 256);
//	for(int i=0; i<nImageValid; i++)
//	{   
//		int nImageId = camParas[i].nImageIndex; 
//		sprintf(filenames[i], "%s", filenamesAll[nImageId]);
//		printf("%s \n", filenames[i]);
//	}
//
//	//find the ground plane 
//	//(just like absolute pose estimation and transform the exterior parameters and tracks to new coordinate )
//	stAbsPOS absPosParams;
//	absPosParams.scale = 1;
//	memset(absPosParams.T, 0, 3*sizeof(double));
//	memset(absPosParams.R, 0, 9*sizeof(double));
//	absPosParams.R[0] = absPosParams.R[4] = absPosParams.R[8] = 1;
//   	{
//		FindGroundPlane(absPosParams, camParas, tracks);
//	}
//
//	//progress write
//	FILE *f_pro = fopen("process_4.txt","wt");	
//	fprintf(f_pro,"10\n");
//	fclose(f_pro);
//	
//	//calculate the mean height of plane
//	double planeHeiMean = 0;
//	int nValidCameras = 0;
//	for(int i=0; i<camParas.size(); i++)
//	{
//		if(camParas[i].f == 0 )
//			continue;
//		planeHeiMean += camParas[i].zs;
//		nValidCameras ++;
//	}
//	planeHeiMean /= double( nValidCameras );
//	printf("plane mean height: %lf \n", planeHeiMean);
//
//	//calculate the area and height of points after absolute pose estimation
//	double minx,maxx,miny,maxy;
//	double meanHeight;
//	CalculateGroudArea(tracks, minx, maxx, miny, maxy, meanHeight);
//	double relativeHei = fabs(meanHeight - planeHeiMean);
//    
//  
//	//calculate the resolution after absolute pose estimation
//	double rawResolution;
//	rawResolution = CalculateResolution(ht, wd, meanHeight, camParas);
//	
//	
//	//get the resolution
//	double outResolution = rawResolution*(nLevel*2 + 1); 
//	int oht = (maxy-miny) / outResolution ;
//	int owd = (maxx-minx) / outResolution ;
//	//get the maximum dimension for 32bit OS
//	double ratio = CalculateRatioFor32bitOS(oht, owd, 1);
//	oht *= ratio;
//	owd *= ratio;
//	outResolution = outResolution / ratio;
//    
//	//progress write
//	f_pro = fopen("process_4.txt","wt");	
//	fprintf(f_pro,"20\n");
//	fclose(f_pro);
//
//
//	//DEM interpolation
//	float* iz = NULL;
//	int demHt, demWd;
//	ratio = CalculateRatioFor32bitOS(oht, owd, 4);
//	int demScale = (outResolution/rawResolution)*(1.0/ratio)*2;
//	int np = tracks.size();
//	double* px = (double*)malloc(np*sizeof(double));
//	double* py = (double*)malloc(np*sizeof(double));
//	double* pz = (double*)malloc(np*sizeof(double));
//	int index = 0;
//	double heiRatio = 0.25;
//	for(int i=0; i<np; i++)
//	{
//		//remove error 3d point
//		double dif = fabs(tracks[i].z - meanHeight);
//        if(dif > fabs(relativeHei*heiRatio) )
//			continue;
//
//		px[index] = tracks[i].x;
//		py[index] = tracks[i].y;
//		pz[index] = tracks[i].z;
//		index ++;
//	}
//	MBAInterpolation(px, py, pz, index, demScale, rawResolution, demHt, demWd, &iz);
//	printf("dem ht*wd: %d %d \n", demHt, demWd);
//	free(px);
//	free(py);
//	free(pz);
//
//	//progress write
//	f_pro = fopen("process_4.txt","wt");	
//	fprintf(f_pro,"30\n");
//	fclose(f_pro);
//
//
//	/*
//	stGeoInfo geoinfo;
//	memset( &geoinfo, 0, sizeof(stGeoInfo) );
//	geoinfo.left = minx;
//	geoinfo.top = miny;
//	geoinfo.dx = rawResolution*demScale;
//	geoinfo.dy = rawResolution*demScale;
//	GdalWriteFloat("d:\\dem.tif", iz, demHt, demWd, geoinfo);
//	*/
//	
//	//generating orthoimage
//	//GenerateMosaicImage(ht, wd, filenames, rawResolution, outResolution, demHt, demWd, demScale, 
//	//	iz, camParas, absPosParams, minx, maxx, miny, maxy, meanHeight, outFile);
//	
//	//progress write
//	f_pro = fopen("process_4.txt","wt");	
//	fprintf(f_pro,"100\n");
//	fclose(f_pro);
//    
//	
//    free(iz);	
//	
//	return 0;
//}
////////////////////////////////////////////////////////////////////////////
//

/* mosaic uav images according the bundler output , written by xie donghai, 2015.9.30
   inputs
	imageListFile: the file containing image absolute path
	bundleFile:    the output of bundler
	outFile:       mosaic image file
	nLevel:        0-high, 1-medium, 2-low
*/
int MosaicOnBundleWithDEM(char* imageListFile, char* bundleFile, 
						  char* smoothedFile,
						  char* posFile, 
						  int nIsGeoCorrection, char* outFile, int nLevel)
{
	printf("[MosaicOnBundleWithDEM] ... \n");

	//reading image list file path
	int nImage = GetFileRows(imageListFile);
	if(nImage<1)
	{
		printf(" %s : image list file is empty, exit! \n", imageListFile);
		return 1;
	}

	char** filenames = f2c(nImage, 256);
	double initFocus;
	int tn;
	FILE* fp = fopen(imageListFile, "r");
	if(fp==NULL)
	{
		printf("error reading image list file! \n");
		return -1;
	}
	for(int i=0; i<nImage; i++)
	{
		fscanf(fp, "%s %d %lf ", filenames[i], &tn, &initFocus);
	}
	fclose(fp);
	
	if(nImage<2)
	{
		printf("image number is not enough ! \n");
		return -1;
	}

	/*
	//get the height and width of image
	int ht,wd;
	IplImage* pImage = cvLoadImage(filenames[0], 0);
	if(pImage==NULL)
	{
		printf("Loading %s failed ! \n", filenames[0]);
		return -1;
	}
	ht = pImage->height;
	wd = pImage->width;
	cvReleaseImage(&pImage);
	*/

	
	//using gdal to get the image information
	stGeoInfo geoinfo;
	GetGeoInformation(filenames[0], geoinfo);
	int ht = geoinfo.ht;
	int wd = geoinfo.wd;
	printf("ht:%d wd:%d \n", ht, wd);
    

	//reading bundle file
	vector<stPOS> camParas;
	int res = ReadBundlerOutFileBak(bundleFile, camParas);
	if(res<0)
	{
		printf("error to read bundle out file! \n");
		return -1;
	}


	//reading smoothed point cloud
	vector<stTrack> tracks;
	ReadSmoothedPointCloud(smoothedFile, tracks);	
	WriteProgressValueToFile(5.0);
	
	//getchar();
	
	//reading POS parameters 
	bool bIsHavePos = false;
	if(nIsGeoCorrection)
	{
		if(posFile != NULL && strlen(posFile)>5 ) //from txt file
		{
			int res = ReadingPosFile(posFile, imageListFile, camParas);
			if(res==0)
			{
				bIsHavePos = true;
			}
		}
		else //from  geo tag from images
		{	
			int nPOSTag = 0;
						
			//nPOSTag = ReadingPosFromJpegImages(imageListFile, camParas);
						
			if( nPOSTag>3 )
				bIsHavePos = true;
		}
	}
	
	//filter out the wrong cameras
    for(int i=0; i<camParas.size(); i++)
	{
		if( camParas[i].f == 0 )
			continue;

		double dif = fabs( camParas[i].f - initFocus );
		if( dif > (initFocus*3) )
			camParas[i].f = 0;
	}
	
	//calculate the valid number of camera
	int nValidCamera = 0;
	for(int i=0; i<camParas.size(); i++)
	{
		if( camParas[i].f != 0 )
			nValidCamera++;
	}
	if(nValidCamera<3)
	{
		printf("the number of valid cameras is less than 3, exit! \n");
		return -1;
	}

	//absolute orientation
	stAbsPOS absPosParams;
	absPosParams.scale = 1;
	memset(absPosParams.T, 0, 3*sizeof(double));
	memset(absPosParams.R, 0, 9*sizeof(double));
	absPosParams.R[0] = absPosParams.R[4] = absPosParams.R[8] = 1;
    if(bIsHavePos)
	{
		//AbsPosEstimation(absPosParams, camParas, tracks);
		double err = AbsOriOrthogonal(absPosParams, camParas, tracks);
		printf("absolute orientation error: %lf \n", err);
	}
	else
	{
		printf("without absolute orientation ... \n");
		//find the ground plane 
		//(just like absolute pose estimation and transform the exterior parameters and tracks to new coordinate )
		FindGroundPlane(absPosParams, camParas, tracks);
	}
	WriteProgressValueToFile(5.0);
	

	//calculate the mean height of plane
	double planeHeiMean = 0;
	int nValidCameras = 0;
	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f == 0 )
			continue;
		planeHeiMean += camParas[i].zs;
		nValidCameras ++;
	}
	planeHeiMean /= double( nValidCameras );
	printf("plane mean height: %lf \n", planeHeiMean);


	//calculate the area and height of points after absolute pose estimation
	double minx,maxx,miny,maxy;
	double meanHeight;
	CalculateGroudArea(tracks, minx, maxx, miny, maxy, meanHeight);
	double relativeHei = fabs(meanHeight - planeHeiMean);
	printf("relative Height: %lf \n", relativeHei);

	//calculate the resolution after absolute pose estimation
	double rawResolution;
	rawResolution = CalculateResolution(ht, wd, meanHeight, camParas);
	printf("raw Resolution: %lf \n", rawResolution);

			
	//set the orthoimage resolution
	double outResolution = rawResolution*1.25;   //sqrt(2.0); //*(nLevel*2 + 1); 	
	int oht = (maxy-miny) / outResolution ;
	int owd = (maxx-minx) / outResolution ;
	
	
	//get the maximum dimension for 32bit OS
	double ratio = CalculateRatioFor32bitOS(oht, owd, 1);
	oht *= ratio;
	owd *= ratio;
	outResolution = outResolution / ratio;

	//DEM interpolation
	float* iz = NULL;
	int demHt, demWd;
	ratio = CalculateRatioFor32bitOS(oht, owd, 4);

	int demScale = (outResolution/rawResolution)*(1.0/ratio)*2;
	int np = tracks.size();
	double* px = (double*)malloc(np*sizeof(double));
	double* py = (double*)malloc(np*sizeof(double));
	double* pz = (double*)malloc(np*sizeof(double));
	int index = 0;
	double heiRatio = 0.25;
	for(int i=0; i<np; i++)
	{
		//remove error 3d point
		double dif = fabs(tracks[i].z - meanHeight);
        if(dif > fabs(relativeHei*heiRatio) )
			continue;

		px[index] = tracks[i].x;
		py[index] = tracks[i].y;
		pz[index] = tracks[i].z;
		index ++;
	}
	int level = 10;
	MBAInterpolation(px, py, pz, index, demScale, rawResolution, level, demHt, demWd, &iz);
	printf("dem ht*wd: %d %d \n", demHt, demWd);
	free(px);
	free(py);
	free(pz);
	WriteProgressValueToFile(10.0);

	
//#ifdef _DEBUG
//	//SaveBmpGeneral(iz, demHt, demWd, "d:\\dem.bmp");
//	stGeoInfo geoinfo;	
//	memset(&geoinfo, 0, sizeof(stGeoInfo));
//	geoinfo.left = minx;
//	geoinfo.top = maxy;
//	geoinfo.dx = rawResolution*demScale;
//	geoinfo.dy = rawResolution*demScale;
//	GdalWriteFloat("d:\\dem.tif", iz, demHt, demWd, geoinfo);
//#endif
	

	/*
	stGeoInfo geoinfo;
	memset( &geoinfo, 0, sizeof(stGeoInfo) );
	geoinfo.left = minx;
	geoinfo.top = miny;
	geoinfo.dx = rawResolution*demScale;
	geoinfo.dy = rawResolution*demScale;
	GdalWriteFloat("d:\\dem.tif", iz, demHt, demWd, geoinfo);
	*/
	
    //mosaic without fusion
	//GenerateMosaicImage(ht, wd, filenames, rawResolution, outResolution, demHt, demWd, demScale, 
	//	iz, camParas, absPosParams, minx, maxx, miny, maxy, meanHeight, outFile);
	

	printf("generating orthoimage.... \n");
	//mosaic with fusion
	GenerateMosaicImageWithFusion(filenames, rawResolution, outResolution, demHt, demWd, demScale, 
		iz, camParas, absPosParams, minx, maxx, miny, maxy, meanHeight, outFile );

    free(iz);

	printf("Finished ! \n");	
	return 0;
}

void FitPlane1(double* px, double* py, double* pz, int np, double* R)
{	
	int i,j,k;

	// 拟合平面，旋转到与XZ平面平行
	double *errorinfo = new double[np];
	double *errorpoints = new double[np];
	double *tmpotx = new double[np];
	double *tmpoty = new double[np];
	double *tmpotz = new double[np];
	double para[100];
	int paranum = 100;
	int infonum = 0;
	FILE* hp=NULL;
	//char  tmppath[256] = "d:\\R.txt";

	FitObject(2,px,py,pz,np,
		para,paranum,
		errorinfo,errorpoints,
		tmpotx,tmpoty,tmpotz,infonum);

	/*
	//写文件	
	hp = fopen(tmppath,"wt");
	for (i=0;i<4;i++)
	{
		fprintf(hp, "%d %lf %lf %lf\n", i, tmpotx[i],tmpoty[i],tmpotz[i]);
	}
	fprintf(hp, "%lf %lf %lf %lf\n",para[0],para[1],para[2],para[3]);
	*/

	int pairarray[100];
	pairarray[0] = 31;
	int objectinfo[3] = {1,0,0};
	Normalize(para,3);
	double tx[1] = {para[0]};
	double ty[1] = {para[1]};
	double tz[1] = {para[2]};
	ComputeParameter(tx, ty, tz, pairarray, objectinfo, 1, para, paranum);
	//fprintf(hp, "%lf %lf %lf\n  %lf %lf %lf\n %lf %lf %lf\n", 
	//	para[0],para[1],para[2],para[3],para[4],para[5],para[6],para[7],para[8]);
	//fclose(hp);

	//rotation matrix from original data to horizontal plane
	for(i=0; i<9; i++)
		R[i] = para[i];	

	delete []errorinfo;
	delete []errorpoints;
	delete []tmpotx;
	delete []tmpoty;
	delete []tmpotz;
}




/*
generate orthoimage, mosaic and fuse, written by xiedonghai, 2015.10.13
*/
/*
int  GenerateOrthoImage(vector<string> filenames,
	double  rawResolution, double  outResolution,
	int  demHt,	int  demWd,	double  demScale,
	float* iz,	vector<CameraPara>& camParas, stAbsPOS absPosParams,
	double  minx, double maxx, double miny, double maxy,
	double  meanHeight,	string outpath)
{
	int zoneNumber = GetZoneNumber(camParas);
	
	//for progress bar
	int nValidate = 0;
	for (int i = 0; i<camParas.size(); i++)
	{
		if (camParas[i].focalLen > 0)
			nValidate++;
	}
	double dstep = 80.0 / double(nValidate);
	//////////////////////////////////////////////////////////////////////////

	for (int i = 0; i<camParas.size(); i++)
	{
		if (camParas[i].focalLen == 0)
			continue;

		char* orthofile = (char*)malloc(256);
		GenerateProductFile(filenames[i], "product", "jpeg", &orthofile);
		char* geofile = (char*)malloc(256);
		GenerateProductFile(filenames[i], "product", "geo", &geofile);

		//if( bIsFileExit(orthofile) && bIsFileExit(geofile) )
		//	continue;

		char* postfix = NULL;
		GetPostfix(filenames[i], &postfix);

		if (strcmp(postfix, "jpg") == 0 || strcmp(postfix, "JPG") == 0)
		{
			//if input format is jpg, then we compress the orthoimage as jpg and save the geodata into the *.geo file
			GenerateSingleZipOrtho(filenames[i], rawResolution, outResolution,
				demHt, demWd, demScale, iz, minx, miny, maxx, maxy, camParas[i], meanHeight, zoneNumber,
				orthofile, geofile);
		}
		else if (strcmp(postfix, "tif") == 0 || strcmp(postfix, "TIF") == 0)
		{
			//if input format is tif, we save the orthoimage as tif file directly
			GenerateProductFile(filenames[i], "product", "tif", &orthofile);
			GenerateSingleOrthoTiff(filenames[i], rawResolution, outResolution,
				demHt, demWd, demScale, iz, minx, miny, maxx, maxy, camParas[i], meanHeight, zoneNumber,
				orthofile);
		}
		else
		{

		}

		free(orthofile);
		free(geofile);

		WriteProgressValueToFile(dstep);
	}

	return 0;
}
*/