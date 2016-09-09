

#include "opticalflow.hpp"


//opencv lib
#include "opencv2/video/tracking.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/core.hpp"


#include <cstdio>
#include <iostream>

using namespace cv;
using namespace std;

#define APP_NAME "simpleflow_demo : "


bool readOpticalFlowFromFile(FILE* file, Mat& flow) 
{
	char header[5];
	if (fread(header, 1, 4, file) < 4 && (string)header != "PIEH") {
		return false;
	}

	int cols, rows;
	if (fread(&cols, sizeof(int), 1, file) != 1||
		fread(&rows, sizeof(int), 1, file) != 1) {
			return false;
	}

	flow = Mat::zeros(rows, cols, CV_32FC2);

	for (int i = 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			Vec2f flow_at_point;
			if (fread(&(flow_at_point[0]), sizeof(float), 1, file) != 1 ||
				fread(&(flow_at_point[1]), sizeof(float), 1, file) != 1) {
					return false;
			}
			flow.at<Vec2f>(i, j) = flow_at_point;
		}
	}

	return true;
}


// binary file format for flow data specified here:
// http://vision.middlebury.edu/flow/data/
void writeOpticalFlowToFile(const Mat& flow, FILE* file) 
{
	int cols = flow.cols;
	int rows = flow.rows;

	fprintf(file, "PIEH");

	if (fwrite(&cols, sizeof(int), 1, file) != 1 ||
		fwrite(&rows, sizeof(int), 1, file) != 1) {
			printf(APP_NAME "writeOpticalFlowToFile : problem writing header\n");
			exit(1);
	}

	for (int i= 0; i < rows; ++i) {
		for (int j = 0; j < cols; ++j) {
			Vec2f flow_at_point = flow.at<Vec2f>(i, j);

			if (fwrite(&(flow_at_point[0]), sizeof(float), 1, file) != 1 ||
				fwrite(&(flow_at_point[1]), sizeof(float), 1, file) != 1) {
					printf(APP_NAME "writeOpticalFlowToFile : problem writing data\n");
					exit(1);
			}
		}
	}
}

static Scalar randomColor(RNG& rng) 
{   
	int icolor = (unsigned)rng;   
	return Scalar(icolor&255, (icolor>>8)&255, (icolor>>16)&255); 
}

void HarrisDetect(char* filename, vector<int>& featx, vector<int>& featy, int number)
{
	IplImage* img_copy= cvLoadImage(filename); //cvCloneImage(img);  
	IplImage* img_gray=cvCreateImage(cvGetSize(img_copy),IPL_DEPTH_8U,1);  
	IplImage* eig_img=cvCreateImage(cvGetSize(img_copy),IPL_DEPTH_32F,1);  
	IplImage* temp_img=cvCloneImage(eig_img);  

	if(img_copy->nChannels==3)
	{
		cvCvtColor(img_copy,img_gray,CV_BGR2GRAY);  
	}
	else
	{
		cvReleaseImage(&img_gray);
		img_gray = cvCloneImage(img_copy);
	}

	const int MAX_CORNERS=number;   //定义角点个数最大值  
	CvPoint2D32f* corners=new CvPoint2D32f[MAX_CORNERS];//分配保存角点的空间  
	int corner_count=MAX_CORNERS;  
	double quality_level=0.01;      //or 0.01  
	double min_distance=10;  
	cvGoodFeaturesToTrack(img_gray,
		eig_img,temp_img,corners,&corner_count,
		quality_level,min_distance);  

	for(int i=0;i<corner_count;i++)  
	{  
		int x = corners[i].x;
		int y = corners[i].y;
		featx.push_back(x);
		featy.push_back(y);
	}

	cvReleaseImage(&img_copy);
	cvReleaseImage(&img_gray);
	cvReleaseImage(&eig_img);
	cvReleaseImage(&temp_img);
	delete[] corners;

}


void HarrisDetect(IplImage* img, vector<int>& featx, vector<int>& featy, int number)
{
	IplImage* img_copy=cvCloneImage(img);  
	IplImage* img_gray=cvCreateImage(cvGetSize(img),IPL_DEPTH_8U,1);  
	IplImage* eig_img=cvCreateImage(cvGetSize(img),IPL_DEPTH_32F,1);  
	IplImage* temp_img=cvCloneImage(eig_img);  

	cvCvtColor(img,img_gray,CV_BGR2GRAY);  
	const int MAX_CORNERS=number;   //定义角点个数最大值  
	CvPoint2D32f* corners=new CvPoint2D32f[MAX_CORNERS];//分配保存角点的空间  
	int corner_count=MAX_CORNERS;  
	double quality_level=0.01;      //or 0.01  
	double min_distance=10;  
	cvGoodFeaturesToTrack(img_gray,
		eig_img,temp_img,corners,&corner_count,
		quality_level,min_distance);  

	for(int i=0;i<corner_count;i++)  
	{  
		int x = corners[i].x;
		int y = corners[i].y;
		featx.push_back(x);
		featy.push_back(y);
	}

	cvReleaseImage(&img_copy);
	cvReleaseImage(&img_gray);
	cvReleaseImage(&eig_img);
	cvReleaseImage(&temp_img);
	delete[] corners;
}


int SimpleFlow(char* lfile, char* rfile, char* outfile, 
			   int& ht, int& wd, 
			   vector<double>& tx,  
			   vector<double>& ty)
{
	//IplImage* pImageLeft = cvLoadImage(lfile);
	//printf("%d %d \n", pImageLeft->height, pImageLeft->width);
	//IplImage* pImageRight = cvLoadImage(rfile);
	//printf("%d %d \n", pImageLeft->height, pImageLeft->width);

	Mat frame1 = imread(lfile); //(pImageLeft,0);
	Mat frame2 = imread(rfile);

	if (frame1.empty()) {
		printf(APP_NAME "Image #1 : %s cannot be read\n", lfile);
		exit(1);
	}

	if (frame2.empty()) {
		printf(APP_NAME "Image #2 : %s cannot be read\n", rfile);
		exit(1);
	}

	if (frame1.rows != frame2.rows && frame1.cols != frame2.cols) {
		printf(APP_NAME "Images should be of equal sizes\n");
		exit(1);
	}

	if (frame1.type() != 16 || frame2.type() != 16) {
		printf(APP_NAME "Images should be of equal type CV_8UC3\n");
		exit(1);
	}

	printf(APP_NAME "Read two images of size [rows = %d, cols = %d]\n",
		frame1.rows, frame1.cols);

	Mat flow;

	printf("simpleflow .... \n");

	float start = (float)getTickCount();
	calcOpticalFlowSF(frame1, frame2, flow,	2, 2, 2);
	/*
		4.1, 
		25.5, 
		18, 
		55.0, 
		25.5, 
		0.35, 
		18, 
		55.0, 
		25.5, 
		10);
		*/

	printf(APP_NAME "calcOpticalFlowSF : %lf sec\n", (getTickCount() - start) / getTickFrequency());

	FILE* file = fopen(outfile, "wb");
	if (file == NULL) {
		printf(APP_NAME "Unable to open file '%s' for writing\n", outfile);
		exit(1);
	}
	printf(APP_NAME "Writing to file\n");
	writeOpticalFlowToFile(flow, file);
	fclose(file);

	//output to the array
	wd = flow.cols;
	ht = flow.rows;
	tx.resize(wd*ht);
	ty.resize(wd*ht);
	RNG rng(0xFFFFFFFF);
	for (int i= 0; i < ht; ++i) 
	{
		for (int j = 0; j < wd; ++j) 
		{
			Vec2f flow_at_point = flow.at<Vec2f>(i, j);
			tx[i*wd+j] = flow_at_point[0];
			ty[i*wd+j] = flow_at_point[1];
              
#if 0
			if(i%30==0 && j%30==0)
			{
				Point cl,cr;
				cl.x = j;
				cl.y = i;
				cr.x = j+tx[i*wd+j];
				cr.y = i+ty[i*wd+j];
				//circle(frame1, cl, 2, randomColor(rng));
				//circle(frame2, cr, 2, randomColor(rng));
				circle(frame1, cl, 2, CV_RGB(255,0,0), 2);
				circle(frame2, cr, 2, CV_RGB(255,0,0), 2);
			}
#endif
		}
	}
	
	//harris feature detect 
	vector<int> featx;
	vector<int> featy;
	HarrisDetect(lfile, featx, featy, 512);

	for(int i=0; i<featx.size(); i++)
	{
		Point cl,cr;
		cl.x = featx[i];
		cl.y = featy[i];
		cr.x = cl.x+tx[cl.y*wd+cl.x];
		cr.y = cl.y+ty[cl.y*wd+cl.x];
		circle(frame1, cl, 2, CV_RGB(255,0,0), 2);
		//arrowedLine(frame1, cl, cr, CV_RGB(0,255,0));
		circle(frame2, cr, 2, CV_RGB(255,0,0), 2);
	}

#if 1
	imwrite("d:\\panaleft.jpg", frame1);
	imwrite("d:\\panaright.jpg", frame2);
#endif

	printf("simple flow finished! \n");

	return 0;
}

