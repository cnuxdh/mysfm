
#ifndef MYUI_HPP
#define MYUI_HPP


#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"

using namespace cv;

#include"export.hpp"

DLL_EXPORT int SaveAsPNG(char* outfile, Mat& image);

DLL_EXPORT int CalculateDstRect(int ht, int wd, double a, double b, double c, double d,
	int& minx, int& miny, int& maxx, int& maxy);

DLL_EXPORT int ImageSimilarityWarp(Mat& srcImage, Mat& dstImage,
	double a, double b, double c, double d,
	int& minx, int& miny, int& maxx, int& maxy, Mat& mask);



DLL_EXPORT int WarpVideoFrame(char* outfile, Mat& srcImage,
	double a, double b, double c, double d);

//DLL_EXPORT int WarpVideoFrame(Mat& srcImage,
//	double a, double b, double c, double d,
//	Mat& dstImage, int& minx, int& miny, int& maxx, int& maxy);
//


#endif
