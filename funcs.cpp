
#include "funcs.hpp"


double distanceVec(Point2DDouble p1, Point2DDouble p2)
{
	double dis = 0;

	double dx = p1.p[0] - p2.p[0];
	double dy = p1.p[1] - p2.p[1];

	dis = sqrt( dx*dx + dy*dy );

	return dis;
}

double norm(Point3DDouble p)
{
	return sqrt( p.p[0]*p.p[0] + p.p[1]*p.p[1] + p.p[2]*p.p[2] );
}

double dot(Point3DDouble p1, Point3DDouble p2)
{
	return ( p1.p[0]*p2.p[0] + p1.p[1]*p2.p[1] + p1.p[2]*p2.p[2] );
}

double angleOfVector(Point3DDouble p1, Point3DDouble p2)
{
	double n1 = norm(p1);
	double n2 = norm(p2);
	double angle = acos( dot(p1,p2)/(n1*n2) );

	angle = angle / PI * 180;

	return angle;
}

//mosaic the image for feature matching
IplImage* VerticalMosaic(IplImage* pLeft, IplImage* pRight)
{	
	int nLeftWd = pLeft->width;
	int nLeftHt = pLeft->height;
	int nScanWdLeft = pLeft->widthStep;

	int nRightWd = pRight->width;
	int nRightHt = pRight->height;
	int nScanWdRight = pRight->widthStep;

	int oht = nLeftHt + nRightHt;
	int owd = max(nLeftWd, nRightWd);

	IplImage* pMosaicImage = cvCreateImage( cvSize(owd,oht), 8, 3);
	int oscanWd = pMosaicImage->widthStep;


	for(int j=0; j<nLeftHt; j++)
		for(int i=0; i<nLeftWd; i++)
		{
			pMosaicImage->imageData[j*oscanWd + i*3]   = pLeft->imageData[j*nScanWdLeft + i*3];
			pMosaicImage->imageData[j*oscanWd + i*3+1] = pLeft->imageData[j*nScanWdLeft + i*3+1];
			pMosaicImage->imageData[j*oscanWd + i*3+2] = pLeft->imageData[j*nScanWdLeft + i*3+2];
		}

		for(int j=0; j<nRightHt; j++)
			for(int i=0; i<nRightWd; i++)
			{
				pMosaicImage->imageData[(j+nLeftHt)*oscanWd + i*3]   = pRight->imageData[j*nScanWdRight + i*3];
				pMosaicImage->imageData[(j+nLeftHt)*oscanWd + i*3+1] = pRight->imageData[j*nScanWdRight + i*3+1];
				pMosaicImage->imageData[(j+nLeftHt)*oscanWd + i*3+2] = pRight->imageData[j*nScanWdRight + i*3+2];
			}

			return pMosaicImage;
}



int DrawMatches(char* filename, IplImage* pLeft, IplImage* pRight, vector<Point2DDouble> lpts, vector<Point2DDouble> rpts)
{
	IplImage* pMosaic = VerticalMosaic(pLeft, pRight);

	int ht = pLeft->height;
	int wd = pLeft->width;

	srand(0);

	for(int i=0; i<lpts.size(); i++)
	{				
		int r,g,b;

		{
			r = (double)(rand())/(double)(RAND_MAX)*255;
			g = (double)(rand())/(double)(RAND_MAX)*255;
			b = (double)(rand())/(double)(RAND_MAX)*255;

			CvPoint pl;
			pl.x = lpts[i].p[0];
			pl.y = lpts[i].p[1];
			CvPoint pr;
			pr.x = rpts[i].p[0];
			pr.y = rpts[i].p[1];	

			//remove the match without motion
			//double dis = sqrt( (double)( (pl.x-pr.x)*(pl.x-pr.x) + (pl.y-pr.y)*(pl.y-pr.y)) );
			//if(dis<5) continue;

			//cvDrawCircle(pLeft, pl, 2, CV_RGB(r,g,b),2);
			//cvDrawCircle(rImage, pr, 2, CV_RGB(r,g,b),2);
			cvDrawLine(pMosaic, pl, cvPoint( pr.x, pr.y+ht ), CV_RGB(r,g,b));
		}				
	}	

	cvSaveImage(filename, pMosaic);
	cvReleaseImage(&pMosaic);

	return 0;
}




