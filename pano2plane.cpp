
#include "pano2plane.h"

#include "absOri.hpp"
#include "baselib.h"

#include "panorama.hpp"

#include "Matrix.h"



/*
  generate perspective image of one direction
*/
IplImage*  PanoToPlane(IplImage* panoImage,double* direction, double vangle, double hangle)
{
	double R[9];
	double xT[3] = {0,0,1};

	double nd = sqrt( direction[0]*direction[0] + direction[1]*direction[1] + direction[2]*direction[2] );

	//calculate the rotation between spherical coordinate to projection coordinate
	POINT3D zp;
	zp.x = -direction[0] / nd;
	zp.y = -direction[1] / nd;
	zp.z = -direction[2] / nd;
	POINT3D xp,yp;
	GenerateProjectAxis(zp, xp, yp);
	vector<POINT3D> srcPts;
	vector<POINT3D> dstPts;
	srcPts.resize(3);
	dstPts.resize(3);
	srcPts[0].x = 1; srcPts[0].y = 0; srcPts[0].z = 0;
	srcPts[1].x = 0; srcPts[1].y = 1; srcPts[1].z = 0;
	srcPts[2].x = 0; srcPts[2].y = 0; srcPts[2].z = 1;  
	dstPts[0] = xp;
	dstPts[1] = yp;
	dstPts[2] = zp;

	//orthogonal Procrustes algorithm
	RotationAlign(srcPts, dstPts, R); //from projection image space to sphere space

	int ht = panoImage->height;
	int wd = panoImage->width;
	int scanWd = panoImage->widthStep;
	double radius = double(wd)/(2*PI);
	//printf("radius: %lf \n", radius);

	//calculate the projection plane size
	double radianAngle = 1.0 / 180.0 * PI;
	double projFocus = radius;    //radius*0.8; 
	//printf("focal length: %.4lf \n", radius);
	//printf("proj focalLen:   %.4lf \n", projFocus);

	int projHt = projFocus*tan(vangle*radianAngle*0.5)*2;
	int projWd = projFocus*tan(hangle*radianAngle*0.5)*2;

	//printf("proj width:%d   height:%d \n", projWd, projHt);

	//image re-projection
	int nBand = panoImage->nChannels;
	IplImage* planeImage = cvCreateImage( cvSize(projWd, projHt), 8, nBand);
	int projScanWd = planeImage->widthStep;

	double grd[3];
	grd[2] = -projFocus;  //z axis

	for(int y=-projHt*0.5; y<projHt*0.5; y++)
	{
		for(int x=-projWd*0.5; x<projWd*0.5; x++)
		{  
			//calculate the image coordinates
			int pj = -y + projHt*0.5;
			int pi =  x + projWd*0.5;
			pj = max( 0, min(projHt-1, pj) );
			pi = max( 0, min(projWd-1, pi) );

			//3D point coordinates
			grd[0] = x;
			grd[1] = y;

			//from projection image space to sphere space
			double rg[3];
			mult(R, grd, rg, 3, 3, 1);

			//for sphere 3D to panorama image space
			double ix, iy;
			GrdToSphere_center( rg[0], rg[1], rg[2], radius, ix, iy);	
			//printf("%lf %lf \n", ix, iy);
			ix = ix + wd*0.5;
			iy = ht*0.5 - iy;

			if(ix>=wd) ix=wd-1;
			if(iy>=ht) iy=ht-1;
			if(ix<0) ix = 0;
			if(iy<0) iy = 0;

			int nx = ix;
			int ny = iy;

			if(nBand==3)
			{
				planeImage->imageData[pj*projScanWd + pi*3 ]    = panoImage->imageData[ny*scanWd+nx*3];
				planeImage->imageData[pj*projScanWd + pi*3 + 1] = panoImage->imageData[ny*scanWd+nx*3 + 1];
				planeImage->imageData[pj*projScanWd + pi*3 + 2] = panoImage->imageData[ny*scanWd+nx*3 + 2];
			}
			else if(nBand==1)
			{
				planeImage->imageData[pj*projScanWd + pi ]    = panoImage->imageData[ny*scanWd+nx];
			}
		}
	}

	//printf("pano to plane projection Finished! \n");

	return planeImage;
}

/*

outputs:
	pR:  from spherical 3D to image 3D

*/
IplImage*  PanoToPlane(IplImage* panoImage, double  vangle, double hangle, 
	double* direction, double focalLenRatio, 
	double& focalLen, int& outHt, int& outWd, double* pR)
{

	double R[9];
	double xT[3] = {0,0,1};

	//calculate the rotation between spherical coordinate to projection coordinate
	POINT3D zp;
	zp.x = direction[0];
	zp.y = direction[1];
	zp.z = direction[2];
	POINT3D xp,yp;
	GenerateProjectAxis(zp, xp, yp);
	vector<POINT3D> srcPts;
	vector<POINT3D> dstPts;
	srcPts.resize(3);
	dstPts.resize(3);
	srcPts[0].x = 1; srcPts[0].y = 0; srcPts[0].z = 0;
	srcPts[1].x = 0; srcPts[1].y = 1; srcPts[1].z = 0;
	srcPts[2].x = 0; srcPts[2].y = 0; srcPts[2].z = 1;  
	dstPts[0] = xp;
	dstPts[1] = yp;
	dstPts[2] = zp;
	
	//orthogonal Procrustes algorithm
	RotationAlign(srcPts, dstPts, R); //from projection image space to sphere space
	
	
	int ht = panoImage->height;
	int wd = panoImage->width;
	int scanWd = panoImage->widthStep;
	double radius = double(wd)/(2*PI);
	//printf("radius: %lf \n", radius);

	//calculate the projection plane size
	double radianAngle = 1.0 / 180.0 * PI;
	double projFocus = radius*focalLenRatio;    //radius*0.8; 
	//printf("focal length: %.4lf \n", radius);
	//printf("proj focalLen:   %.4lf \n", projFocus);

	int projHt = projFocus*tan(vangle*radianAngle*0.5)*2;
	int projWd = projFocus*tan(hangle*radianAngle*0.5)*2;
	
	//printf("proj width:%d   height:%d \n", projWd, projHt);

	//image re-projection
	int nBand = panoImage->nChannels;
	IplImage* planeImage = cvCreateImage( cvSize(projWd, projHt), 8, nBand);
	int projScanWd = planeImage->widthStep;

	double grd[3];
	grd[2] = -projFocus;  //z axis

	for(int y=-projHt*0.5; y<projHt*0.5; y++)
	{
		for(int x=-projWd*0.5; x<projWd*0.5; x++)
		{  
			//calculate the image coordinates
			int pj = -y + projHt*0.5;
			int pi =  x + projWd*0.5;
			pj = max( 0, min(projHt-1, pj) );
			pi = max( 0, min(projWd-1, pi) );

			//pj = projHt - pj - 1;

			//3D point coordinates
			grd[0] = x;
			grd[1] = y;

			//from projection image space to sphere space
			double rg[3];
			mult(R, grd, rg, 3, 3, 1);

			//for sphere 3D to panorama image space
			double ix, iy;
			GrdToSphere_center( rg[0], rg[1], rg[2], radius, ix, iy);	
			//printf("%lf %lf \n", ix, iy);
			ix = ix + wd*0.5;
			iy = ht*0.5 - iy;

			if(ix>=wd) ix=wd-1;
			if(iy>=ht) iy=ht-1;
			if(ix<0) ix = 0;
			if(iy<0) iy = 0;

			int nx = ix;
			int ny = iy;

			if(nBand==3)
			{
				planeImage->imageData[pj*projScanWd + pi*3 ]    = panoImage->imageData[ny*scanWd+nx*3];
				planeImage->imageData[pj*projScanWd + pi*3 + 1] = panoImage->imageData[ny*scanWd+nx*3 + 1];
				planeImage->imageData[pj*projScanWd + pi*3 + 2] = panoImage->imageData[ny*scanWd+nx*3 + 2];
			}
			else if(nBand==1)
			{
				planeImage->imageData[pj*projScanWd + pi ]    = panoImage->imageData[ny*scanWd+nx];
			}
		}
	}

	//generate the projection rotation
	memcpy(pR, R, sizeof(double)*9);
	invers_matrix(pR, 3);

	//get the inner parameters of projection camera
	focalLen = projFocus;
	outHt = projHt;
	outWd = projWd;

	printf("pano to plane projection Finished! \n");

	return planeImage;
}


int PanoToPlanes(IplImage* panoImage, double anglestep,
	double vangle, double hangle, double fratio,
	double* R, double* T, 
	vector<IplImage*>& projImages,
	vector<CameraPara>& camParas)
{
	//projection from panorama to plane 
	int nAngleStep  = anglestep;
	int nProjectNum = 360 / nAngleStep;
	for(int i=0; i<nProjectNum; i++)
	{

		double lat = 90 / 180.0 * PI; //from top to bottom: 0-180
		double lon = i*nAngleStep / 180.0 * PI; //from left to right: 0-360

		//calculate the direction according to (lon,lat), must be the opposite
		double direction[3];
		direction[0] = -sin(lon)*sin(lat);  //x
		direction[1] = -cos(lon)*sin(lat);  //y
		direction[2] = -cos(lat);           //z


		//printf("direction: %lf %lf %lf \n", direction[0], direction[1], direction[2]);

		double Rp[9]; //rotation matrix for plane 
		//generate the plane projection image and save it
		double focalLen;
		int outHt, outWd;
		IplImage* planeImage = PanoToPlane( panoImage, vangle, hangle, direction, fratio, 
			focalLen, outHt, outWd, Rp);

		/*for(int kj=0; kj<3; kj++)
		{
			for(int ki=0; ki<3; ki++)
			{
				printf("%6.4lf ", Rp[kj*3+ki]);
			}
			printf("\n");
		}*/

		projImages.push_back(planeImage);

		//transform from spherical space to image space,
		//Xs = Rg.Xg + Tg, Xp = Rp.Xs ---> Xp = Rp.Rg.Xg + Rp.Tg, define Rpg = Rp.Rg, Tpg=Rp.Tg
		double Rpg[9];
		double Tpg[3];
		mult(Rp, R, Rpg, 3, 3, 3);
		mult(Rp, T, Tpg, 3, 3, 1);
		
		//save the R and T into the array
		CameraPara cam;
		cam.focalLen = focalLen;
		cam.k1 = 0;
		cam.k2 = 0;
		cam.rows = outHt;
		cam.cols = outWd;
		memcpy(cam.R, Rpg, 9*sizeof(double));
		memcpy(cam.T, Tpg, 3*sizeof(double));

		camParas.push_back(cam);
	}

	return 0;
}


/* panorama to one plane projection
inputs:
   srcImageFile: panoram file 
   outImageFile: output file
   vangle:    vertical fov angle
   hangle:    horizonal fov angle
   direction: the plane image normal direction
   focalLenRatio:  focalLen / panorama radius
output:
	 focalLen: 
	 outHt,outWd:
	 pR: rotation from the spherical coordinate to projection plane coordinate
   
*/
int PanoToPlane(char* srcImageFile, char* outImageFile, 
				double vangle, double hangle,
				double* direction, double focalLenRatio, 
				double& focalLen, int& outHt, int& outWd, double* pR 
				)
{

	double R[9];
	double xT[3] = {0,0,1};

	//calculate the rotation between spherical coordinate to projection coordinate
	//first method: R.xT = direction
	//CalculateAlignRotation1(xT, direction, R);
	//second method: orthogonal Procrustes algorithm

	POINT3D zp;
	zp.x = direction[0];
	zp.y = direction[1];
	zp.z = direction[2];
	POINT3D xp,yp;
	GenerateProjectAxis(zp, xp, yp);
	vector<POINT3D> srcPts;
	vector<POINT3D> dstPts;
	srcPts.resize(3);
	dstPts.resize(3);
	srcPts[0].x = 1; srcPts[0].y = 0; srcPts[0].z = 0;
	srcPts[1].x = 0; srcPts[1].y = 1; srcPts[1].z = 0;
	srcPts[2].x = 0; srcPts[2].y = 0; srcPts[2].z = 1;  
	dstPts[0] = xp;
	dstPts[1] = yp;		
	dstPts[2] = zp;
	RotationAlign(srcPts, dstPts, R); //from projection image space to sphere space

	/*
	char rotatefile[256];
	strcpy(rotatefile, outImageFile);
	strcat(rotatefile, ".txt");
	FILE* fp = fopen(rotatefile, "w");
	for(int j=0; j<3; j++)
	{
		for(int i=0; i<3; i++)
		{
			fprintf(fp, "%lf ", R[j*3+i]);		
		}
		fprintf(fp, "\n");
	}
	fclose(fp);
	*/


	//char   filename[256];	
	//strcpy(filename, argv[1]);
	//double vangle = atof(argv[2]); //fov vertical angle
	//double hangle = atof(argv[3]); //fov horizonal angle
	//double widthAngle = PI*(1.0/3.0); // 60 degree
	//printf("input image: %s \n", filename);
	//printf("input angle: %lf %lf \n", vangle,hangle);

	IplImage* pImage = cvLoadImage(srcImageFile, 1);
	//cvNamedWindow("panoimage");
	//cvShowImage("panoimage", pImage);

	int ht = pImage->height;
	int wd = pImage->width;
	int scanWd = pImage->widthStep;
	double radius = double(wd)/(2*PI);
	//printf("radius: %lf \n", radius);

	//calculate the projection plane size
	double radianAngle = 1.0 / 180.0 * PI;
	double projFocus = radius*focalLenRatio;    //radius*0.8; 
	printf("focal length: %.4lf \n", radius);
	printf("proj focalLen:   %.4lf \n", projFocus);

	int projHt = projFocus*tan(vangle*radianAngle*0.5)*2;
	int projWd = projFocus*tan(hangle*radianAngle*0.5)*2;


	printf("proj width:%d   height:%d \n", projWd, projHt);

	//image reprojection
	IplImage* planeImage = cvCreateImage( cvSize(projWd, projHt), 8, 3);
	int projScanWd = planeImage->widthStep;

	double grd[3];
	grd[2] = -projFocus;  //z axis

	for(int y=-projHt*0.5; y<projHt*0.5; y++)
	{
		for(int x=-projWd*0.5; x<projWd*0.5; x++)
		{  
			//calculate the image coordinates
			int pj = -y + projHt*0.5;
			int pi =  x + projWd*0.5;
			pj = max( 0, min(projHt-1, pj) );
			pi = max( 0, min(projWd-1, pi) );

			//pj = projHt - pj - 1;

			//3D point coordinates
			grd[0] = x;
			grd[1] = y;

			//from projection image space to sphere space
			double rg[3];
			mult(R, grd, rg, 3, 3, 1);

			//for sphere 3D to panoram image spance
			double ix, iy;
			GrdToSphere_center( rg[0], rg[1], rg[2], radius, ix, iy);	
			//printf("%lf %lf \n", ix, iy);

			ix = ix + wd*0.5;
			iy = ht*0.5 - iy;

			//ix = min(wd-1, ix);
			//iy = min(ht-1, iy);
			if(ix>=wd) ix=wd-1;
			if(iy>=ht) iy=ht-1;
			if(ix<0) ix = 0;
			if(iy<0) iy = 0;

			int nx = ix;
			int ny = iy;
			planeImage->imageData[pj*projScanWd + pi*3 ]    = pImage->imageData[ny*scanWd+nx*3];
			planeImage->imageData[pj*projScanWd + pi*3 + 1] = pImage->imageData[ny*scanWd+nx*3 + 1];
			planeImage->imageData[pj*projScanWd + pi*3 + 2] = pImage->imageData[ny*scanWd+nx*3 + 2];
		}
	}

	char file[256];
	cvSaveImage(outImageFile, planeImage);

	//cvNamedWindow("projimage");
	//cvShowImage("projimage", planeImage);

	cvReleaseImage(&planeImage);
	cvReleaseImage(&pImage);

	//generate the projection rotation
	memcpy(pR, R, sizeof(double)*9);
	invers_matrix(pR, 3);

	//get the inner parameters of projection camera
	focalLen = projFocus;
	outHt = projHt;
	outWd = projWd;

	printf("pano to plane projection Finished! \n");

	return 0;
}


int PanoToPlanes(char* srcFile, double anglestep, double vangle, 
	             double hangle, double fratio, char* outpath)
{
	char title[256];
	strcpy(title, srcFile);
	strcpy(title+strlen(title)-4, "\0");
	
	//projection from panorama to plane 
	int nAngleStep  = anglestep;
	int nProjectNum = 360 / nAngleStep;
	for(int i=0; i<nProjectNum; i++)
	{
		double lat = 90 / 180.0 * PI; //from top to bottom: 0-180
		double lon = i*nAngleStep / 180.0 * PI; //from left to right: 0-360

		//calculate the direction according to (lon,lat), must be the opposite
		double direction[3];
		direction[0] = -sin(lon)*sin(lat);  //x
		direction[1] = -cos(lon)*sin(lat);  //y
		direction[2] = -cos(lat);           //z

		printf("direction: %lf %lf %lf \n", direction[0], direction[1], direction[2]);

		char outfile[256];
		//sprintf(outfile, "%s_%d_%d.jpg", title, int(lat/PI*180+0.5), int(lon/PI*180+0.5));
		sprintf(outfile, "%s-%.2d.jpg", title, i);
		printf("plane file: %s \n", outfile);
		
		double Rp[9]; //rotation matrix for plane 
		//generate the plane projection image and save it
		double focalLen;
		int outHt, outWd;
		PanoToPlane( srcFile, outfile, vangle, hangle, direction, fratio, 
			focalLen, outHt, outWd, Rp);
	}

	return 0;
}

/* projection from sphere space to several planes
inputs:
  nImageIndex: 0,1,2.... N, the image index
	srcFile: the panorama file
	R:  the rotation matrix of the panorama 
	T:  the translation vector of the panorama
outputs:
	projected plane image series
*/
int PanoToPlanes( int nImageIndex, char* srcFile, double anglestep,
				  double vangle, double hangle, double fratio,
				  double* R, double* T,
				  vector<CameraPara>& camParas,
				  char* outpath)
{

	char title[256];
	strcpy(title, srcFile);
	strcpy(title+strlen(title)-4, "\0");
	
	//projection from panorama to plane 
	int nAngleStep  = anglestep;
	int nProjectNum = 360 / nAngleStep;
	//double vangle   = 60;
	//double hangle   = 60;
	//double ratio = 1;

	//nProjectNum = 1;

	for(int i=0; i<nProjectNum; i++)
	{

		double lat = 90 / 180.0 * PI; //from top to bottom: 0-180
		double lon = i*nAngleStep / 180.0 * PI; //from left to right: 0-360

		//printf("input image: %s \n", filename);
		//printf("input angle: %lf %lf \n", vangle,hangle);

		//calculate the direction according to (lon,lat), must be the opposite
		double direction[3];
		direction[0] = -sin(lon)*sin(lat);  //x
		direction[1] = -cos(lon)*sin(lat);  //y
		direction[2] = -cos(lat);           //z


		printf("direction: %lf %lf %lf \n", direction[0], direction[1], direction[2]);


		char outfile[256];
		//sprintf(outfile, "%s_%d_%d.jpg", title, int(lat/PI*180+0.5), int(lon/PI*180+0.5));
		sprintf(outfile, "%s\\visualize\\%.8d.jpg", outpath, (nImageIndex*nProjectNum+i) );
		//sprintf(outfile, "%s\\visualize\\%.8d.jpg", outpath, nImageIndex);
		printf("plane file: %s \n", outfile);


		double Rp[9]; //rotation matrix for plane 
		//generate the plane projection image and save it
		double focalLen;
		int outHt, outWd;
		PanoToPlane( srcFile, outfile, vangle, hangle, direction, fratio, 
			focalLen, outHt, outWd, Rp);

		//transform from spherical space to image space,
		//Xs = Rg.Xg + Tg, Xp = Rp.Xs ---> Xp = Rp.Rg.Xg + Rp.Tg, define Rpg = Rp.Rg, Tpg=Rp.Tg
		double Rpg[9];
		double Tpg[3];
		mult(Rp, R, Rpg, 3, 3, 3);
		mult(Rp, T, Tpg, 3, 3, 1);


		//generate the projection matrix for PMVS and save it	  
		char projFile[256];
		//sprintf(projFile, "%s_%d_%d.txt", title, int(lat/PI*180+0.5), int(lon/PI*180+0.5));
		sprintf(projFile, "%s\\txt\\%.8d.txt", outpath, (nImageIndex*nProjectNum+i) );
		//sprintf(projFile, "%s\\txt\\%.8d.txt", outpath, nImageIndex);

		double K[9] = 
		{ -focalLen, 0.0,       0.5 * outWd - 0.5,
		  0.0,       focalLen,  0.5 * outHt - 0.5,
		  0.0,       0.0,       1.0 };

		double Ptmp[12] = 
		{ Rpg[0], Rpg[1], Rpg[2], Tpg[0],
		  Rpg[3], Rpg[4], Rpg[5], Tpg[1],
		  Rpg[6], Rpg[7], Rpg[8], Tpg[2] };

		double P[12];
		dll_matrix_product(3, 3, 3, 4, K, Ptmp, P);
		dll_matrix_scale(3, 4, P, -1.0, P);

		FILE* f = fopen(projFile, "w");
		fprintf(f, "CONTOUR\n");
		fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[0], P[1], P[2],  P[3]);
		fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[4], P[5], P[6],  P[7]);
		fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[8], P[9], P[10], P[11]);
		fprintf(f, "%lf %d %d \n", focalLen, outHt, outWd);
		//append the R and T for other applications
		fprintf(f, "%lf %lf %lf \n", Tpg[0], Tpg[1], Tpg[2]);
		fprintf(f, "%lf %lf %lf \n", Rpg[0], Rpg[1], Rpg[2]);
		fprintf(f, "%lf %lf %lf \n", Rpg[3], Rpg[4], Rpg[5]);
		fprintf(f, "%lf %lf %lf \n", Rpg[6], Rpg[7], Rpg[8]);
		fclose(f);

		//save the R and T into the array
		CameraPara cam;
		cam.focalLen = focalLen;
		cam.k1 = 0;
		cam.k2 = 0;
		cam.rows = outHt;
		cam.cols = outWd;
		memcpy(cam.R, Rpg, 9*sizeof(double));
		memcpy(cam.T, Tpg, 3*sizeof(double));

		camParas.push_back(cam);
	}

	return 0;
}
