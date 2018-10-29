
#include"myui.hpp"

//cvlib
#include "geotiff.h"


int SaveAsPNG(char* outfile, Mat& image)
{
	int ht = image.rows;
	int wd = image.cols;

	unsigned char* pr = new unsigned char[ht*wd];
	unsigned char* pg = new unsigned char[ht*wd];
	unsigned char* pb = new unsigned char[ht*wd];
	unsigned char* alpha = new unsigned char[ht*wd];
	memset(alpha, 0, ht*wd);

	for (int i = 0; i<ht; i++)
	{
		uchar* data = image.ptr<uchar>(i);
		for (int j = 0; j<wd; j++)
		{
			pr[i*wd + j] = data[j * 3];
			pg[i*wd + j] = data[j * 3 + 1];
			pb[i*wd + j] = data[j * 3 + 2];

			int ts = data[j * 3] + data[j * 3 + 1] + data[j * 3 + 2];
			if (ts>0)
				alpha[i*wd + j] = 255;
		}
	}

	GdalWritePNGImage(outfile, pb, pg, pr, alpha, ht, wd);

	free(pr);
	free(pg);
	free(pb);
	free(alpha);

	return 0;
}


int CalculateDstRect(int ht, int wd, double a, double b, double c, double d,
	int& minx, int& miny, int& maxx, int& maxy)
{

	double px[4] = { -wd*0.5, wd*0.5, wd*0.5, -wd*0.5 };
	double py[4] = { ht*0.5, ht*0.5, -ht*0.5, -ht*0.5 };

	minx = 1000000;
	maxx = -1000000;
	miny = 1000000;
	maxy = -1000000;

	for (int i = 0; i < 4; i++)
	{
		int tx = px[i] * a - py[i] * b + c;
		int ty = px[i] * b + py[i] * a + d;

		minx = min(tx, minx);
		maxx = max(tx, maxx);
		miny = min(ty, miny);
		maxy = max(ty, maxy);
	}
		
	return 0;
}



int ImageSimilarityWarp(Mat& srcImage, Mat& dstImage, 
	double a, double b, double c, double d,
	int& minx, int& miny, int& maxx, int& maxy, Mat& mask)
{
	double ht = srcImage.rows;
	double wd = srcImage.cols;

	double px[4] = { -wd*0.5, wd*0.5, wd*0.5, -wd*0.5 };
	double py[4] = { ht*0.5, ht*0.5, -ht*0.5, -ht*0.5 };

	int dstHt = dstImage.rows;
	int dstWd = dstImage.cols;

	int nChannel = dstImage.channels();

	//inverse resample
	double lamda = 1.0 / (a*a+b*b);
	for (int j = 0; j < dstHt; j++)
	{
		uchar* dstdata = dstImage.ptr<uchar>(j);
		uchar* maskdata = mask.ptr<uchar>(j);

		for (int i = 0; i < dstWd; i++)
		{
			double y = (dstHt*0.5-j) + miny + dstHt*0.5;
			double x = (i-dstWd*0.5) + minx + dstWd*0.5;

			int ix = (a*(x - c) + b*(y - d))*lamda;
			int iy = (-b*(x - c) + a*(y - d))*lamda;

			ix = ix + wd*0.5;
			iy = ht*0.5 - iy;

			//ix = max(0, min(ix, int(wd)));
			//iy = max(0, min(iy, int(ht)));

			if (ix >= 0 && ix < wd && iy >= 0 && iy < ht)
			{
				uchar* srcdata = srcImage.ptr<uchar>(iy);

				dstdata[i * nChannel] = srcdata[ix * 3];
				dstdata[i * nChannel + 1] = srcdata[ix * 3 + 1];
				dstdata[i * nChannel + 2] = srcdata[ix * 3 + 2];
				if (nChannel == 4) 
					dstdata[i * nChannel + 3] = 255;

				//mask.at<unsigned char>(j, i) = 1;
				maskdata[i] = 1;
			}		
		}
	}
	
	return 0;
}



int WarpVideoFrame(char* outfile, Mat& srcImage, 
	double a, double b, double c, double d)
{
	int minx, miny, maxx, maxy;
	int ht = srcImage.rows;
	int wd = srcImage.cols;
	CalculateDstRect(ht, wd, a, b, c, d, minx, miny, maxx, maxy);
	
	int dstHt = maxy - miny;
	int dstWd = maxx - minx;
	Mat dstImage(dstHt, dstWd, CV_8UC3, Scalar(0, 0, 0));
	Mat mask(dstHt, dstWd, CV_8UC1, Scalar(0));
	ImageSimilarityWarp(srcImage, dstImage, a, b, c, d, minx, miny, maxx, maxy, mask);
	
	//imwrite("c:\\temp\\warp.jpg", dstImage);
	SaveAsPNG(outfile, dstImage);

	return 0;
}


/*
int WarpVideoFrame(Mat& srcImage,
	double a, double b, double c, double d, 
	Mat& dstImage, int& minx, int& miny, int& maxx, int& maxy)
{
	//int minx, miny, maxx, maxy;
	int ht = srcImage.rows;
	int wd = srcImage.cols;
	CalculateDstRect(ht, wd, a, b, c, d, minx, miny, maxx, maxy);

	int dstHt = maxy - miny;
	int dstWd = maxx - minx;
	dstImage(dstHt, dstWd, CV_8UC4, Scalar(0, 0, 0, 0));
	ImageSimilarityWarp(srcImage, dstImage, a, b, c, d, minx, miny, maxx, maxy);
	

	return 0;
}
*/