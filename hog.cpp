
#include "defines.hpp"
#include "hog.hpp"

#include <vector>
using namespace std;


#ifdef OPENCV_1X 
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#else
//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc_c.h"
#endif

//calculate the HOG of 2*2 cells for image 
int CalculateHOG(IplImage* pImage, vector<double>& hog)
{
	//degree
	const int angleStep = 20;
    int histDim = 360 / angleStep;

	IplImage* pSrc = cvCloneImage(pImage);
	
	int nchannels = pSrc->nChannels;
	if(nchannels==3)
	{
		IplImage* pGray = cvCreateImage(cvGetSize(pSrc), IPL_DEPTH_8U, 1);//创建目标图像  
		cvCvtColor(pSrc, pGray, CV_BGR2GRAY);//cvCvtColor(src,des,CV_BGR2GRAY)  
		cvReleaseImage(&pSrc);
		pSrc = cvCloneImage(pGray);
		cvReleaseImage(&pGray);
	}
	
	int ht = pSrc->height;
	int wd = pSrc->width;
	int scanWd = pSrc->widthStep;

    //smooth the image
	cvSmooth(pSrc, pSrc);
    
	int x[3] = {0, wd*0.5, wd-1};
	int y[3] = {0, ht*0.5, ht-1};

	double dSum = 0;
	for(int m=0; m<2; m++)
	{
		for(int n=0; n<2; n++)
		{
			int t = y[m];
			int b = y[m+1];
			int l = x[n];
			int r = x[n+1];

			//calculate dx, dy	
			vector<double> cellHist;
			cellHist.resize(histDim);

			for(int j=t; j<=b; j++)
			{
				for(int i=l; i<=r; i++)
				{
					int p0 = (unsigned char)(pSrc->imageData[j*scanWd+i]);
					int pl = (unsigned char)(pSrc->imageData[j*scanWd+i+1]);
					int pb = (unsigned char)(pSrc->imageData[(j+1)*scanWd+i]);
					
					double dx = pl - p0;
					double dy = pb - p0;

					//ignore the area without textures
					if( fabs(dx)<4 && fabs(dy)<4 )
						continue;
					
					double gradAngle = atan2(dy, dx);
					
					if(gradAngle<0) gradAngle += 2*PI;
					gradAngle = gradAngle / PI * 180;
					double dIndex = gradAngle / angleStep;

					int nLeftIndex = int(dIndex);
					int nRightIndex = (nLeftIndex+1) % histDim;

					cellHist[nLeftIndex]  += 1-(dIndex-nLeftIndex);
					cellHist[nRightIndex] += (dIndex-nLeftIndex);

					dSum ++;
				}
			}

			for(int k=0; k<cellHist.size(); k++)
			{
				hog.push_back(cellHist[k]);
			}
		}
	}

	//normalize the hist
	for(int i=0; i<hog.size(); i++)
	{
		hog[i] /= dSum;
	}

	cvReleaseImage(&pSrc);

	return 0;
}
