
#include"stdio.h"
#include"string.h"
#include"stdlib.h"
#include"math.h"

#include "panorama.hpp"
#include "defines.hpp"
//corelib
#include "matrix.h"


//convert spherical image coordinate to 3D points
//x,y: original point is the center of image
int SphereTo3D_center(double x, double y, double radius, double& gx, double& gy, double& gz)
{
	int wd = 2*PI*radius;
    int ht = wd*0.5;

	double lon = x/radius; //(x - wd*0.5) / radius;
	double lat = y/radius; //(ht*0.5 - y) / radius;

	gx = radius*sin(lon)*cos(lat);
	gy = radius*cos(lon)*cos(lat);
	gz = radius*sin(lat);
	
	return 0;
}

//convert spherical image coordinate to 3D points
//x,y: image col and row, original point is the top-left corner of the image
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

	ix = radius*sita; //radius*sita + 2*PI*radius ;
	iy = radius*fai;  //-radius*fai + PI*radius ;

	return 0;
}

//from 3D point to spherical coordinate, the image origin is the top-left 
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

/* generate the points lying on the panorama epipolar line  
   input
		normal: the normal of epipolar plane
   return
		3d point vectors of epipolar line 
*/
vector<Point3DDouble> GenerateEpipolarPlaneVectors(Point3DDouble normal,int num)
{
	vector<Point3DDouble> epipolarVecs;

	double xs = normal.p[0];
	double ys = normal.p[1];
	double zs = normal.p[2];

	//the initial point vertical to the (xs, ys, zs)
	Point3DDouble initialPt;
	initialPt.p[0] = -zs;   //-ys;
	initialPt.p[1] = 0;  //xs;
	initialPt.p[2] = xs;  

	Point3DDouble fixedPt;
	fixedPt.p[0] = xs;
	fixedPt.p[1] = ys;
	fixedPt.p[2] = zs;

	//rotation 
	double interval = 360.0 / num;
	for(int i=0; i<num; i++)
	{
		double angle = -i*interval;

		Point3DDouble vec = RotateByVector(initialPt, fixedPt, angle );

		epipolarVecs.push_back(vec);
	}

	return epipolarVecs;
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

DLL_EXPORT int GeneratePanoEpipolarImageHeading(double* R, double* T, IplImage* pLeftImage, IplImage* pRightImage)
{
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
		//printf("%d \n", j);

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
			gx = res[0];	gy=-res[2];	gz=res[1];

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
		//printf("%d \n", j);

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
			og[0]=res[0];	og[1]=-res[2];	og[2]=res[1];

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

	cvSaveImage("c:\\temp\\epipolarLeftHeading.jpg",  pEpipolarLeft);
	cvSaveImage("c:\\temp\\epipolarRightHeading.jpg", pEpipolarRight);

	cvReleaseImage(&pEpipolarLeft);
	cvReleaseImage(&pEpipolarRight);

	printf(" \n GeneratePanoEpipolarImage(..) is finished ! \n");

	return 0;
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
		//printf("%d \n", j);

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
			gx = res[0];	gy=-res[2];	gz=res[1];

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
		//printf("%d \n", j);

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
			og[0]=res[0];	og[1]=-res[2];	og[2]=res[1];

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

	cvSaveImage("c:\\temp\\epipolarLeftHeading.jpg",  pEpipolarLeft);
	cvSaveImage("c:\\temp\\epipolarRightHeading.jpg", pEpipolarRight);

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

	cvSaveImage("c:\\temp\\epipolarLeft.jpg", pEpipolarLeft);
	cvSaveImage("c:\\temp\\epipolarRight.jpg", pEpipolarRight);

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



/*
CIntegratedPanoMatch::CIntegratedPanoMatch()
{
	m_IsEssentialMatrixReady = false;
}

CIntegratedPanoMatch::~CIntegratedPanoMatch()
{

}

int CIntegratedPanoMatch::InitEpipolarConstraint(CameraPara left, CameraPara right )
{	
	//calculate the relative pose according to the POS
	m_leftPanoCam = left;
	m_rightPanoCam = right;
	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, 
		m_rightPanoCam.R, m_rightPanoCam.t, 
		m_leftPanoCam.bIsExplicit, 
		m_R, m_T, m_EM);
	
	m_IsEssentialMatrixReady = true;
	return 1;
}

int CIntegratedPanoMatch::Match(IplImage* pLeftImage, IplImage* pRightImage, vector<ImagePair>& matches)
{
	//
	////if we do not know the relative pose between two images
	//if(!m_IsEssentialMatrixReady)
	//{
	//	IplImage* pLeftSmall  = ResizeImage(pLeftImage,  PANORAMA_HT_SMALL);
	//	IplImage* pRightSmall = ResizeImage(pLeftImage,  PANORAMA_HT_SMALL);

	//	//relative pose estimation directly using panorama
	//	CFeatureBase* pFeatDetect = new CSIFTFloat();
	//	ImgFeature lFeats;
	//	pFeatDetect->Detect(pLeftSmall, lFeats);
	//	ImgFeature rFeats;
	//	pFeatDetect->Detect(pLeftSmall,rFeats);
	//	delete pFeatDetect;
	//	cvReleaseImage(&pLeftSmall);
	//	cvReleaseImage(&pRightSmall);

	//	CMatchBase* pMatch = new CPanoMatch();
	//	PairMatchRes mr;
	//	pMatch->Match(lFeats, rFeats, mr);
	//	int matchNum = mr.matchs.size();
	//	vector<Point2DDouble> lpts,rpts;
	//	lpts.resize(matchNum);
	//	rpts.resize(matchNum);
	//	for(int i=0; i<matchNum; i++)
	//	{
	//		int li = mr.matchs[i].l;
	//		int ri = mr.matchs[i].r;
	//		//default origin of feature point is at the center of image
	//		lpts[i] = lFeats.GetCenteredPt(li);
	//		rpts[i] = rFeats.GetCenteredPt(ri);
	//	}
	//	delete pMatch;

	//	//relative pose estimation
	//	CRelativePoseBase* pRP = new CEstimatePose5PointPano();
	//	m_leftPanoCam.rows  = pLeftSmall->height;
	//	m_leftPanoCam.cols  = pLeftSmall->width;
	//	m_rightPanoCam.rows = pRightSmall->height;
	//	m_rightPanoCam.cols = pRightSmall->width;
	//	pRP->EstimatePose(lpts, rpts, m_leftPanoCam, m_rightPanoCam );  
	//	delete pRP;
	//}
	//
		
	//resize the image
	IplImage* resizeLeft = NULL;
	IplImage* resizeRight = NULL;
	ResizeImage(pLeftImage,  PANORAMA_HT);
	ResizeImage(pRightImage, PANORAMA_HT);


	//split the panorama image into several perspective images 
	double R[9] = {1,0,0,0,1,0,0,0,1};
	double T[3] = {0,0,0};
	vector<IplImage*>  lProjImages;
	vector<CameraPara> lCamParas;
	PanoToPlanes(resizeLeft,  60, 90, 90, 1, R, T, lProjImages, lCamParas);
	vector<IplImage*>  rProjImages;
	vector<CameraPara> rCamParas;
	PanoToPlanes(resizeRight, 60, 90, 90, 1, R, T, rProjImages, rCamParas);


	//detect feature points
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	vector<ImgFeature> lImageFeats;
	lImageFeats.resize(lProjImages.size());
	for(int i=0; i<lProjImages.size(); i++)
	{
		pFeatDetect->Detect(lProjImages[i], lImageFeats[i]);
	} 
	vector<ImgFeature> rImageFeats;
	lImageFeats.resize(rProjImages.size());
	for(int i=0; i<rProjImages.size(); i++)
	{
		pFeatDetect->Detect(rProjImages[i], rImageFeats[i]);
	} 
	
	//matching one by one between left and right perspective images

	for(int j=0; j<lImageFeats.size(); j++)
	{
		double zaxis[3] = {0, 0, -1};
		double Rl[9];
		memcpy(Rl, lCamParas[j].R, sizeof(double)*9);
		invers_matrix(Rl, 3);
		double imageDirLeft[3];
		mult(Rl, zaxis, imageDirLeft, 3, 3, 1);
		Point3DDouble lp;
		lp.p[0] = imageDirLeft[0];
		lp.p[1] = imageDirLeft[1];
		lp.p[2] = imageDirLeft[2];

		int bestIndex = 0;
		int minAngle = 180;
		for(int i=0; i<rImageFeats.size(); i++)
		{
			double Rr[9];
			memcpy(Rr, rCamParas[i].R, sizeof(double)*9);
			invers_matrix(Rr, 3);
			double imageDirRight[3];
			mult(Rr, zaxis, imageDirRight, 3, 3, 1);
			Point3DDouble rp;
			rp.p[0] = imageDirRight[0];
			rp.p[1] = imageDirRight[1];
			rp.p[2] = imageDirRight[2];

			//calculate the included angle
			double angle = angleOfVector(lp, rp);
					
			if(angle<minAngle)
			{
				minAngle = angle;
				bestIndex = i;
			}
		}
	
	}


	return 0;
}


int CIntegratedPanoMatch::Match(IplImage* pLeftImage, IplImage* pRightImage, 
	CameraPara left, CameraPara right, vector<ImagePair>& matches)
{
	//construct the relative pose using POS data
	m_leftPanoCam = left;
	m_rightPanoCam = right;
	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, 
		m_rightPanoCam.R, m_rightPanoCam.t, 
		m_leftPanoCam.bIsExplicit, 
		m_R, m_T, m_EM);
	
	//



	return 0;
}*/