


#ifndef  CV_PANORAMA_H
#define  CV_PANORAMA_H

#include "export.hpp"
#include "defines.hpp"




#include <vector>
using namespace std;


DLL_EXPORT int SphereTo3D_center(double x, double y, double radius, double& gx, double& gy, double& gz);
DLL_EXPORT int SphereTo3D(double x, double y, double radius, double& gx, double& gy, double& gz);

DLL_EXPORT int GrdToSphere_center(double gx, double gy, double gz, double radius, double& ix, double& iy);
DLL_EXPORT int GrdToSphere(double gx, double gy, double gz, double radius, double& ix, double& iy);


//for cylinder projection
//from 3d space to cylinder projection
int GrdToCylinder(double gx, double gy, double gz, double radius, double& ix, double& iy);
//from cylinder projection to 3d space
int CylinderToGrd(double ix, double iy, double radius, double& gx, double& gy, double& gz);



/*  convert spherical panorama image into cylinder projection
*/
DLL_EXPORT int SphereToCilinder(char* infile, char* outfile);





/* panorama to plane projection
   srcImageFile: panoram file 
   outImageFile: output file
   vangle:    vertical fov angle
   hangle:    horizonal fov angle
   direction: the plane image normal direction
   focalLenRatio:  focalLen / panaram radius
   
output:
	 pR: rotation from the spherical coordinate to projection plane coordinate
	 
*/
DLL_EXPORT int PanoToPlane(char* srcImageFile, char* outImageFile, 
				double vangle, double hangle, double* direction, double focalLenRatio, 
				double& focalLen, int& outHt, int& outWd,
				double* pR);


//panoram image is projected to several images
int PanoToPlanes(int nImageIndex, char* srcFile, double anglestep,
									double vangle, double hangle, double fratio,
									double* R, double* T);



/* relative pose estimation for spherical panoramic images, based on 5-point algorithm, 
   written by xdh, 2015.7.8
input:
	lPts,rPts: image points in panoramic image
	radius:    the width of panoramic image
output:
	cam:       save the R and T
*/
//DLL_EXPORT int PanoramicRelativePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, int width, CameraPara& cam);



/* simple method: calculate the rotation matrix align 3D vector a to b  
                Ra = b
   but this method may be wrong!!!!
*/
DLL_EXPORT int CalculateAlignRotation(double* a, double* b, double* R);



/* calculate the rotation matrix align 3D vector a to b  
                Ra = b
*/
DLL_EXPORT int CalculateAlignRotation1(double* a, double* b, double* R);



/*
input: relative estimation results, projection model x = RX + T
*/
DLL_EXPORT int GeneratePanoEpipolarImage(double* R, double* T, char* leftFile, char* rightFile);


DLL_EXPORT int GeneratePanoEpipolarImageHeading(double* R, double* T, char* leftFile, char* rightFile);


#endif

