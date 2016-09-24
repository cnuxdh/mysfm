
#include"stdio.h"
#include"string.h"
#include"stdlib.h"
#include"math.h"


//opencv, these should be put previously, otherwise errors appear
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"


#include "panorama.hpp"
#include "absOri.hpp"
#include "defines.hpp"
//#include "relativepose.hpp"


//corelib
#include "corelib/matrix.h"


//matrix lib
#include "matrix/matrix.h"



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
	printf("proj focus: %lf \n", projFocus);
	
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
             
			//3D point coordinates
			grd[0] = x;
			grd[1] = y;

			//from projection image space to sphere space
			double rg[3];
			mult(R, grd, rg, 3, 3, 1);

			//for sphere 3D to panoram image spance
			double ix, iy;
			GrdToSphere( rg[0], rg[1], rg[2], radius, ix, iy);	
			//printf("%lf %lf \n", ix, iy);

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
									vector<CameraPara>& camParas)
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
		sprintf(outfile, "visualize/%.8d.jpg", (nImageIndex*nProjectNum+i) );
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
		sprintf(projFile, "txt/%.8d.txt", (nImageIndex*nProjectNum+i) );
		
	  double K[9] = 
            { -focalLen, 0.0, 0.5 * outWd - 0.5,
              0.0, focalLen,  0.5 * outHt - 0.5,
              0.0, 0.0, 1.0 };

    double Ptmp[12] = 
        { Rpg[0], Rpg[1], Rpg[2], Tpg[0],
          Rpg[3], Rpg[4], Rpg[5], Tpg[1],
          Rpg[6], Rpg[7], Rpg[8], Tpg[2] };
    
    double P[12];
    matrix_product(3, 3, 3, 4, K, Ptmp, P);
    matrix_scale(3, 4, P, -1.0, P);

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
	  cam.focus = focalLen;
	  cam.k1 = 0;
	  cam.k2 = 0;
	  cam.rows = outHt;
	  cam.cols = outWd;
	  memcpy(cam.R, Rpg, 9*sizeof(double));
	  memcpy(cam.t, Tpg, 3*sizeof(double));
	  
	  camParas.push_back(cam);
	}
	
	return 0;
}







//convert spherical image coordinate to 3D points
//x,y: original point is the center of image
int SphereTo3D_center(double x, double y, double radius, double& gx, double& gy, double& gz)
{
	int wd = 2*PI*radius;
    int ht = wd*0.5;

	double lon = (x - wd*0.5) / radius;
	double lat = (ht*0.5 - y) / radius;

	gx = radius*sin(lon)*cos(lat);
	gy = radius*cos(lon)*cos(lat);
	gz = radius*sin(lat);
	
	return 0;
}


//convert spherical image coordinate to 3D points
//x,y: image col and row
int SphereTo3D(double x, double y, double radius, double& gx, double& gy, double& gz)
{
	
	double lon = x / radius;
	double lat = y / radius;       
	
	//from lon,lat to 3D point
	gx = radius*sin(lon)*sin(lat);
	gy = radius*cos(lon)*sin(lat);
	gz = radius*cos(lat);
	
	
	return 0;
}


//from 3D point to spherical coordinate
int GrdToSphere_center(double gx, double gy, double gz, double radius, double& ix, double& iy)
{
	double sita, fai; //sita-lontitude, fai-latitude

	sita = atan2( gx, gy );
	//sita = atan2( gy, gx );
	//if(sita<0) sita += 2*PI;

	fai  = atan2(gz, sqrt(gx*gx+gy*gy) );
	//fai  = atan2( gz, sqrt(gx*gx+gy*gy) );
	//if(fai<0)  fai += PI;

	ix = radius*sita + 2*PI*radius ;
	iy = -radius*fai + PI*radius ;

	return 0;
}

//from 3D point to spherical coordinate
int GrdToSphere(double gx, double gy, double gz, double radius, double& ix, double& iy)
{
	double sita, fai; //sita-lontitude, fai-latitude
    
	sita = atan2( gx, gy );
	//sita = atan2( gy, gx );
    if(sita<0) sita += 2*PI;

	fai  = atan2( sqrt(gx*gx+gy*gy), gz );
	//fai  = atan2( gz, sqrt(gx*gx+gy*gy) );
	if(fai<0)  fai += PI;
 
	ix = radius*sita;
	iy = radius*fai;

	return 0;
}


//panoramic essential matrix
int CalculatePanoramicEpipolar(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, int width, double* F)
{
	int size = lPts.size();

	//normalized 3D points of sphere
	vector<Point3DDouble> grdLeft;
	vector<Point3DDouble> grdRight;
	grdLeft.resize(size);
	grdRight.resize(size);
	
	//from image coordinates to spherical 
	//calculate the radius of spherical projection 
	double radius = (double)(width) / 2*PI;

	//spherical coordinates to 3D 
	for(int i=0; i<size; i++)
	{
		double x = lPts[i].p[0];
		double y = lPts[i].p[1];
		double gx,gy,gz;
		SphereTo3D(x, y, radius, gx, gy, gz);
		grdLeft[i].p[0] = gx;
		grdLeft[i].p[1] = gy;
		grdLeft[i].p[2] = gy;
		printf("left: %lf %lf %lf \n", gx, gy, gz);

		x = rPts[i].p[0];
		y = rPts[i].p[1];
		SphereTo3D(x, y, radius, gx, gy, gz);
		grdRight[i].p[0] = gx;
		grdRight[i].p[1] = gy;
		grdRight[i].p[2] = gy;
		printf("right: %lf %lf %lf \n", gx, gy, gz);
	}
	return 0;
}

/*
//panoramic relative pose estimation
int PanoramicRelativePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, int width, CameraPara& cam)
{
	int size = lPts.size();
	//the vector points
	vector<Point2DDouble> lVectorPts;
	vector<Point2DDouble> rVectorPts;
	lVectorPts.resize(size);
	rVectorPts.resize(size);
	//calculate the radius of spherical projection 
	double radius = (double)(width) / 2*PI;
	//spherical coordinates to 3D 
	for(int i=0; i<size; i++)
	{
		double x = lPts[i].p[0];
		double y = lPts[i].p[1];
		double gx,gy,gz;
		SphereTo3D(x, y, radius, gx, gy, gz);
		gx = gx / gz;
		gy = gy / gz;
		lVectorPts[i].p[0] = gx;
		lVectorPts[i].p[1] = gy;

		x = rPts[i].p[0];
		y = rPts[i].p[1];
		SphereTo3D(x, y, radius, gx, gy, gz);
		gx = gx / gz;
		gy = gy / gz;
		rVectorPts[i].p[0] = gx;
		rVectorPts[i].p[1] = gy;	
	}
	//relative pose estimation
	CRelativePoseBase* pRP = new CEstimatePose5Point();
	CameraPara cam1;
	CameraPara cam2;	
	cam1.focus = 1;
	cam2.focus = 1;
	pRP->EstimatePose(lVectorPts, rVectorPts, cam1, cam2);
	return 0;
}
*/

void NormalizeVector(double* v, int dim)
{
   double s = 0;
   for(int i=0; i<dim; i++)
	   s += v[i]*v[i];
	s = sqrt(s);

   for(int i=0; i<dim; i++)
	   v[i] /= s;
}

int CalculateAlignRotation(double* a, double* b, double* R)
{
	NormalizeVector(a, 3);
	NormalizeVector(b, 3);

	R[0]=b[0]*a[0]; R[1]=b[0]*a[1]; R[2]=b[0]*a[2]; 
	R[3]=b[1]*a[0]; R[4]=b[1]*a[1]; R[5]=b[1]*a[2]; 
	R[6]=b[2]*a[0]; R[7]=b[2]*a[1]; R[8]=b[2]*a[2]; 	

	return 0;
}

int CalculateAlignRotation1(double* a, double* b, double* R)
{
	double oa[3];
	double ob[3];

	memcpy(oa, a, sizeof(double)*3);
	memcpy(ob, b, sizeof(double)*3);

	NormalizeVector(a, 3);
	NormalizeVector(b, 3);

	double theta = acos(a[0]*b[0] + a[1]*b[1] + a[2]*b[2]);
    if(theta<0.001)
	{
		memset(R, 0, sizeof(double)*9);
		R[0]=1; R[4]=1; R[8]=1;
		return 0;
	}
	
	// to be revised 
	if(  fabs(PI-theta)<0.001)
	{
		memset(R, 0, sizeof(double)*9);
		
		//rotation by Y axis
		R[0]=-1; R[4]=1; R[8]=1;

		return 0;
	}	
    
	//cross product
	double cv[3];
	cv[0] = a[1]*b[2] - b[1]*a[2];
	cv[1] = a[2]*b[0] - b[2]*a[0];
	cv[2] = a[0]*b[1] - b[0]*a[1];
	NormalizeVector(cv, 3);

	//skew symmetric matrix
    double A[9];
	A[0]=0;		 A[1]=-cv[2];  A[2]=cv[1];
	A[3]=cv[2];  A[4]=0;       A[5]=-cv[0];
	A[6]=-cv[1]; A[7]=cv[0];   A[8]=0;

	double A2[9];
	mult(A, A, A2, 3, 3, 3);

	memset(R, 0, sizeof(double)*9);
	R[0]=1; R[4]=1; R[8]=1;

	//exponential map
    for(int i=0; i<9; i++)
		R[i] += sin(theta)*A[i] + (1-cos(theta))*A2[i];
	
	memcpy(a, oa, sizeof(double)*3);
	memcpy(b, ob, sizeof(double)*3);

	return 0;
}

/* rotate srcVector around fixedVector
input:  angle(degree)
return: the rotated vector
*/
Point3DDouble RotateByVector(Point3DDouble srcVector, Point3DDouble fixedVector, double angle)
{
	Point3DDouble dstVector;

	double R[9];
	double radianAngle = angle*DPI;

	double c = cos(radianAngle);
	double s = sin(radianAngle);

	double x = fixedVector.p[0];
	double y = fixedVector.p[1];
	double z = fixedVector.p[2];

	double norm = sqrt(x*x+y*y+z*z);
	x /= norm;
	y /= norm;
	z /= norm;
	
	R[0] = x*x*(1-c)+c; 	R[1] = x*y*(1-c)-z*s;	R[2] = x*z*(1-c)+y*s;
	R[3] = y*x*(1-c)+z*s;	R[4] = y*y*(1-c)+c;		R[5] = y*z*(1-c)-x*s;
	R[6] = x*z*(1-c)-y*s;	R[7] = y*z*(1-c)+x*s;	R[8] = z*z*(1-c)+c;

	double sp[3];
	sp[0] = srcVector.p[0];
	sp[1] = srcVector.p[1];
	sp[2] = srcVector.p[2];

	double res[3];
	mult(R, sp, res, 3, 3, 1);

	dstVector.p[0] = res[0];
	dstVector.p[1] = res[1];
	dstVector.p[2] = res[2];

	return dstVector;
}



/* generate the normals of epipolar planes    
*/
vector<Point3DDouble> GenerateEpipolarPlaneNormals(double xs, double ys, double zs, int num)
{
	vector<Point3DDouble> normals;

	//the initial point vertical to the (xs, ys, zs)
	Point3DDouble initialPt;
	initialPt.p[0] = -ys;
	initialPt.p[1] = xs;
	initialPt.p[2] = 0;
    
	Point3DDouble fixedPt;
	fixedPt.p[0] = xs;
	fixedPt.p[1] = ys;
	fixedPt.p[2] = zs;
	 
	//rotation 
	double interval = 180.0 / num;
	for(int i=0; i<num; i++)
	{
		double angle = i*interval;

		Point3DDouble planeNorm = RotateByVector(initialPt, fixedPt, angle );
	
		normals.push_back(planeNorm);
	}

	return normals;
}

DLL_EXPORT int GeneratePanoEpipolarImageHeading(double* R, double* T, char* leftFile, char* rightFile)
{

	IplImage* pLeftImage = cvLoadImage(leftFile);
	IplImage* pRightImage = cvLoadImage(rightFile);
	int ht = pLeftImage->height;
	int wd = pLeftImage->width;
	int scanwd = pLeftImage->widthStep;

	IplImage* pEpipolarLeft  = cvCloneImage(pLeftImage);
	IplImage* pEpipolarRight = cvCloneImage(pRightImage);

	//calculate the position of second camera relative to the first camera, from RX+T to R(X-T')
	double sT[3];
	double iR[9];
	memcpy(iR, R, sizeof(double)*9);
	invers_matrix(iR, 3);
	mult(iR, T, sT, 3, 3, 1);
	sT[0] = -sT[0];
	sT[1] = -sT[1];
	sT[2] = -sT[2];
	printf("position: %lf %lf %lf \n", sT[0], sT[1], sT[2]);

	double radius = wd / (2*PI);

	double nR[9];
	double a[3];
	double b[3];
	a[0]=0;		a[1]=0;		a[2]=1;
	b[0]=-T[0];	b[1]=T[2];	b[2]=T[1];
	CalculateAlignRotation1(a, b, nR);
	//invers_matrix(nR, 3);
	
	//right epipolar image
	for(int j=0; j<ht; j++)
	{
		printf("%d \n", j);

		for(int i=0; i<wd; i++)
		{
			double ngx,ngy,ngz;
			SphereTo3D( (double)(i), (double)(j), radius, ngx, ngy, ngz );
			
			//transform to the old coordinate
			double gt[3];
			gt[0]=ngx;	gt[1]=ngy;	gt[2]=ngz;
			double res[3];
			mult(nR, gt, res, 3, 3, 1);

			//change projection mode to original
			double gx,gy,gz;
			gx = -res[0];	gy=res[2];	gz=res[1];

			//convert from 3D to spherical coordinate
			double dx,dy;
			GrdToSphere(gx, gy, gz, radius, dx, dy);
			int ix = min(wd-1, int(dx+0.5));
			int iy = min(ht-1, int(dy+0.5));

			pEpipolarRight->imageData[j*scanwd + 3*i]     = pRightImage->imageData[iy*scanwd + 3*ix];
			pEpipolarRight->imageData[j*scanwd + 3*i + 1] = pRightImage->imageData[iy*scanwd + 3*ix + 1];
			pEpipolarRight->imageData[j*scanwd + 3*i + 2] = pRightImage->imageData[iy*scanwd + 3*ix + 2];
		}
	}

	//left epipolar image
	for(int j=0; j<ht; j++)
	{
		printf("%d \n", j);

		for(int i=0; i<wd; i++)
		{
			double ngx,ngy,ngz;
			SphereTo3D((double)(i), (double)(j), radius, ngx, ngy, ngz);
			double gt[3];
			gt[0]=ngx;	gt[1]=ngy;	gt[2]=ngz;
			
			//transform to the old coordinate
			double res[3];
			mult(nR, gt, res, 3, 3, 1);

			//change projection mode to original
			double og[3];
			og[0]=-res[0];	og[1]=res[2];	og[2]=res[1];

			//from right to left
			double lt[3];
			mult(iR, og, lt, 3, 3, 1);
			double gx,gy,gz;
			gx=lt[0];	gy=lt[1];	gz=lt[2];
            
			//convert from 3D to spherical coordinate
			double dx,dy;
			GrdToSphere(gx, gy, gz, radius, dx, dy);
			int ix = min(wd-1, int(dx+0.5));
			int iy = min(ht-1, int(dy+0.5));

			pEpipolarLeft->imageData[j*scanwd + 3*i]     = pLeftImage->imageData[iy*scanwd + 3*ix];
			pEpipolarLeft->imageData[j*scanwd + 3*i + 1] = pLeftImage->imageData[iy*scanwd + 3*ix + 1];
			pEpipolarLeft->imageData[j*scanwd + 3*i + 2] = pLeftImage->imageData[iy*scanwd + 3*ix + 2];
		}
	}

	cvSaveImage("d:\\epipolarLeftHeading.jpg",  pEpipolarLeft);
	cvSaveImage("d:\\epipolarRightHeading.jpg", pEpipolarRight);

	cvReleaseImage(&pLeftImage);
	cvReleaseImage(&pRightImage);

	printf(" \n GeneratePanoEpipolarImage(..) is finished ! \n");

	return 0;
}


//int GeneratePanoEpipolarImage(double* R, double* T, IplImage* pLeftImage, IplImage* pRightImage)
DLL_EXPORT int GeneratePanoEpipolarImage(double* R, double* T, char* leftFile, char* rightFile)
{
	IplImage* pLeftImage = cvLoadImage(leftFile);
	IplImage* pRightImage = cvLoadImage(rightFile);
    

	int ht = pLeftImage->height;
	int wd = pLeftImage->width;
	int scanwd = pLeftImage->widthStep;

	IplImage* pEpipolarLeft  = cvCloneImage(pLeftImage);
	IplImage* pEpipolarRight = cvCloneImage(pRightImage);

	double radius = wd / (2*PI);

	//calculate the position of second camera relative to the first camera, from RX+T to R(X-T')
	double sT[3];
	double iR[9];
	memcpy(iR, R, sizeof(double)*9);

	invers_matrix(iR, 3);
	mult(iR, T, sT, 3, 3, 1);
	sT[0] = -sT[0];
	sT[1] = -sT[1];
	sT[2] = -sT[2];
	printf("position: %lf %lf %lf \n", sT[0], sT[1], sT[2]);
	
     
	//generate the normals 
	/*
	double a = sT[0];
	double b = sT[0];
	double c = sT[0];
	double dIntervals = 180.0 / double(ht);
	vector<double> normalLat;
	vector<double> normalLon;
	vector<Point3DDouble> normals;
	for(int i=0; i<ht; i++)
	{
		//double lat = i*dIntervals;
		double lon = i*dIntervals*DPI;
		double lat = atan( - c / (a*sin(lon)+b*cos(lon)) );

		//sphere to 3D 
		Point3DDouble n;
		n.p[0] = sin(lon)*sin(lat);
		n.p[1] = cos(lon)*sin(lat);
		n.p[2] = cos(lat);
		normals.push_back(n);
	}
	*/

	//generate epipolar lines	
	double dLonInterval = 360.0 / (double)(wd);

	//left image
	vector<Point3DDouble> normals = GenerateEpipolarPlaneNormals(sT[0], sT[1], sT[2], ht); 
	
	for(int i=0; i<normals.size(); i++)
	{
		printf(" %d ", i);

		double a = normals[i].p[0];
		double b = normals[i].p[1];
		double c = normals[i].p[2];

		if( c==0 )
			continue;

		for(int k=0; k<wd; k++)
		{
			double lon = k*dLonInterval*DPI;
			double lat = atan( - c / (a*sin(lon)+b*cos(lon)) );
			if(lat<0) lat += PI;
			
			//image coordinates
			//printf("lon: %lf  lat: %lf  \n", lon, lat);            
			int x = radius*lon;
			int y = radius*lat;

			//calculate the 3D point coordinate
			int nx = k;
			int ny = i;
			/*
			double gx = sin(lon)*sin(lat);
			double gy = cos(lon)*sin(lat);
			double gz = cos(lat);
			//rotate the 3D point
			double newlon,newlat;
			//GrdToSphere(-gx, gz, gy, 1, newlon, newlat); //heading as the Z axis
			GrdToSphere(gx, gy, gz, 1, newlon, newlat); //vertical to heading as the Z axis
			nx = newlon*radius;
			ny = newlat*radius;
			*/
            
			pEpipolarLeft->imageData[ny*scanwd + 3*nx]     = pLeftImage->imageData[y*scanwd + 3*x];
			pEpipolarLeft->imageData[ny*scanwd + 3*nx + 1] = pLeftImage->imageData[y*scanwd + 3*x + 1];
			pEpipolarLeft->imageData[ny*scanwd + 3*nx + 2] = pLeftImage->imageData[y*scanwd + 3*x + 2];
		}
	}	

	//right image
	for(int i=0; i<normals.size(); i++)
	{
		printf(" %d ", i);
		double tN[3];
		tN[0] = normals[i].p[0];
		tN[1] = normals[i].p[1];
		tN[2] = normals[i].p[2];

		double rN[3];
		mult(R, tN, rN, 3, 3, 1);

		double a = rN[0];
		double b = rN[1];
		double c = rN[2];

		if( c==0 )
			continue;

		for(int k=0; k<wd; k++)
		{
			double lon = k*dLonInterval*DPI;
			double lat = atan( - c / (a*sin(lon)+b*cos(lon)) );
            

			if(lat<0) lat += PI;
			//printf("lon: %lf  lat: %lf  \n", lon, lat);            
			int x = radius*lon;
			int y = radius*lat;

			//calculate the 3D point coordinate
			int nx = k;
			int ny = i;
			/*
			double gx = sin(lon)*sin(lat);
			double gy = cos(lon)*sin(lat);
			double gz = cos(lat);
			//rotate the 3D point
			double newlon,newlat;
			GrdToSphere(gx, gy, gz, 1, newlon, newlat);	
			nx = newlon*radius;
			ny = newlat*radius;
			*/

			pEpipolarRight->imageData[ny*scanwd + 3*nx]     = pRightImage->imageData[y*scanwd + 3*x];
			pEpipolarRight->imageData[ny*scanwd + 3*nx + 1] = pRightImage->imageData[y*scanwd + 3*x + 1];
			pEpipolarRight->imageData[ny*scanwd + 3*nx + 2] = pRightImage->imageData[y*scanwd + 3*x + 2];
		}
	}

	cvSaveImage("d:\\epipolarLeft.jpg", pEpipolarLeft);
	cvSaveImage("d:\\epipolarRight.jpg", pEpipolarRight);

	cvReleaseImage(&pLeftImage);
	cvReleaseImage(&pRightImage);
	cvReleaseImage(&pEpipolarLeft);
	cvReleaseImage(&pEpipolarRight);

	printf(" \n GeneratePanoEpipolarImage(..) is finished ! \n");


	return 0;
}


//3d to cylinder projection point
int GrdToCylinder(double gx, double gy, double gz, double radius,
				  double& ix, double& iy)
{
	//
	int wd = 2*PI*radius;
	int ht = wd*0.5;
	
	//normalize the 3D point to the radius
	double ratio = radius / sqrt(gx*gx+gy*gy);
	
	gx *= ratio;
	gy *= ratio;
	gz *= ratio;
	
	//
	iy = ht*0.5 - gz;
    
	//calculate the azimuth angle
	double azi = atan2(gy, gx);
	if(azi<0) azi += 2*PI;
    
	ix = azi*radius;

	return 0;
}


int CylinderToGrd(double ix, double iy, double radius,
				  double& gx, double& gy, double& gz)
{

	int wd = radius*2*PI;
	int ht = wd*0.5;	
	//
	gz = ht*0.5 - iy;
	
	double angle = ix / (double)(wd) * 2 * PI;
	gy = radius*cos(angle);
	gx = radius*sin(angle);
    
	return 0;
}


// convert the sphere panoram image to cylinder projection
DLL_EXPORT int SphereToCilinder(char* infile, char* outfile)
{
	IplImage* inputImage = cvLoadImage(infile);

	int ht = inputImage->height;
	int wd = inputImage->width;
	int scanwd = inputImage->widthStep;

	double radius = double(wd) / (2*PI);

	IplImage* outImage = cvCreateImage( cvSize(wd, ht), 8, 3) ;

	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			double gx,gy,gz;
			//SphereTo3D(i, j, radius, gx, gy, gz);
			CylinderToGrd(i, j, radius, gx, gy, gz);
			

			double ix,iy;
			//GrdToCylinder(gx, gy, gz, radius, ix, iy);
			GrdToSphere(gx, gy, gz, radius, ix, iy);

            if( iy<0 || iy>=ht )
				continue;

			int cx = (int)(ix);
			int cy = (int)(iy);

			int r = (unsigned char)(inputImage->imageData[cy*scanwd+cx*3]);
			int g = (unsigned char)(inputImage->imageData[cy*scanwd+cx*3+1]);
			int b = (unsigned char)(inputImage->imageData[cy*scanwd+cx*3+2]);

			outImage->imageData[j*scanwd+i*3]   = r;
			outImage->imageData[j*scanwd+i*3+1] = g;
			outImage->imageData[j*scanwd+i*3+2] = b;
		}

	cvSaveImage(outfile, outImage);

	cvReleaseImage(&inputImage);
	cvReleaseImage(&outImage);

	return 0;
}