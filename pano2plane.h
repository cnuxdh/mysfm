#ifndef PANO_TO_PLANE_H
#define PANO_TO_PLANE_H


#include "export.hpp"
#include "defines.hpp"

//opencv, these should be put previously, otherwise errors appear
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"





IplImage*  PanoToPlane(IplImage* panoImage,double* direction, double vangle, double hangle);


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
						double  vangle, double hangle, double* direction, double focalLenRatio, 
						double& focalLen, int& outHt, int& outWd,
						double* pR);

DLL_EXPORT int PanoToPlanes(char* srcFile, double anglestep, double vangle, 
							double hangle, double fratio, char* outpath);

//panoram image is projected to several images and save the projection matrix 
DLL_EXPORT int PanoToPlanes(int nImageIndex, char* srcFile, double anglestep,
							double vangle, double hangle, double fratio,
							double* R, double* T, vector<CameraPara>& camParas,
							char* outpath);


DLL_EXPORT IplImage*  PanoToPlane(IplImage* panoImage, double  vangle, double hangle, 
	double* direction, double focalLenRatio, 
	double& focalLen, int& outHt, int& outWd, double* pR);

//panoram image is projected to several images and save the images and projection matrix 
DLL_EXPORT int PanoToPlanes(IplImage* panoImage, double anglestep,
	double vangle, double hangle, double fratio,
	double* R, double* T, 
	vector<IplImage*>& projImages,
	vector<CameraPara>& camParas);


#endif