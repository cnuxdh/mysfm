
#include "interactReg.hpp"

#include "register.hpp"
#include "relativepose.hpp"
#include "sift.hpp"
#include "panorama.hpp"
#include "funcs.hpp"

#include "commondata.h"
#include "CommonFuncs.h"
#include "matrix.h"

#include "delaunay.h"




int DrawCross1(CvPoint p, int len, IplImage* pImage)
{
	cvLine( pImage, cvPoint( p.x-len,p.y), cvPoint( p.x+len,p.y),CV_RGB(255,0,0),3,8,0);
	cvLine( pImage, cvPoint( p.x,p.y-len), cvPoint( p.x,p.y+len),CV_RGB(255,0,0),3,8,0);

	return 0;
}

//generate the histogram data around the pixel(x,y) from image
int GetHistData(int x, int y, IplImage* pImage, int* hist, int nGrayStep, int blockSize)
{
	if(pImage->nChannels==3)
	{
		printf("Only gray image is supported !");
		return -1;
	}

	int ht = pImage->height;
	int wd = pImage->width;
	int scanWd = pImage->widthStep;

	int tx = x - blockSize*0.5;
	int ty = y - blockSize*0.5;
	if(tx<0) tx=0;
	if(ty<0) ty=0;

	int bx = x + blockSize*0.5;
	int by = y + blockSize*0.5;
	if(bx>=wd) bx = wd-1;
	if(by>=ht) by = ht-1;
	
	for(int j=ty; j<=by; j++)
		for(int i=tx; i<=bx; i++)
		{
			int pixelValue = (unsigned char)(pImage->imageData[j*scanWd+i]);
			pixelValue = pixelValue/nGrayStep;
			hist[pixelValue] ++;
		}
}


CIPanoReg::CIPanoReg()
{
	m_pLeft = NULL;
	m_pRight = NULL;
}
CIPanoReg::~CIPanoReg()
{
	cvReleaseImage(&m_pLeft);
	cvReleaseImage(&m_pRight);

	for(int i=0; i<m_pLeftPlaneImages.size(); i++)
		cvReleaseImage(&m_pLeftPlaneImages[i]);

	for(int i=0; i<m_pRightPlaneImages.size(); i++)
		cvReleaseImage(&m_pRightPlaneImages[i]);
}

int CIPanoReg::Init(char* pLeftFile, char* pRightFile)
{
	//loading gray image
	m_pLeft = cvLoadImage(pLeftFile, 0);
	m_pRight = cvLoadImage(pRightFile, 0);
	
	//detect feature points
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	ImgFeature lFeats;
	pFeatDetect->Detect(pLeftFile, lFeats, 480);
	ImgFeature rFeats;
	pFeatDetect->Detect(pRightFile,rFeats, 480);
	delete pFeatDetect;

	//matching
	CMatchBase* pMatch = new CPanoMatch();
	PairMatchRes mr;
	pMatch->Match(lFeats, rFeats, mr);
	int matchNum = mr.matchs.size();
	vector<Point2DDouble> lpts,rpts;
	lpts.resize(matchNum);
	rpts.resize(matchNum);
	for(int i=0; i<matchNum; i++)
	{

		int li = mr.matchs[i].l;
		int ri = mr.matchs[i].r;

		//default origin of feature point is at the center of image
		lpts[i] = lFeats.GetCenteredPt(li);
		rpts[i] = rFeats.GetCenteredPt(ri);
	}

	//relative pose estimation
	CRelativePoseBase* pRP = new CEstimatePose5PointPano();
	//CameraPara leftCam, rightCam;
	m_leftPanoCam.rows  = m_pLeft->height;
	m_leftPanoCam.cols  = m_pLeft->width;
	m_rightPanoCam.rows = m_pRight->height;
	m_rightPanoCam.cols = m_pRight->width;
	pRP->EstimatePose(lpts, rpts, m_leftPanoCam, m_rightPanoCam );  


	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, m_rightPanoCam.R, m_rightPanoCam.t, true, m_EM);

	
	//generate perspective images
	m_pLeftPlaneImages.clear();
	m_leftPlaneCams.clear();
	PanoToPlanes(m_pLeft, 60, 90, 90,  1, m_leftPanoCam.R, m_leftPanoCam.t, m_pLeftPlaneImages, m_leftPlaneCams);
	
	m_pRightPlaneImages.clear();
	m_rightPlaneCams.clear();
	PanoToPlanes(m_pRight, 60, 90, 90, 1, m_rightPanoCam.R, m_rightPanoCam.t, m_pRightPlaneImages, m_rightPlaneCams);

	//for debug
	if(0)
	{
		char file[256];
		for(int i=0; i<m_pLeftPlaneImages.size(); i++)
		{
			sprintf(file, "c:\\temp\\persperctive_left_%d.jpg", i);
			printf("%s \n", file);
			cvSaveImage(file, m_pLeftPlaneImages[i]);

			sprintf(file, "c:\\temp\\persperctive_right_%d.jpg", i);
			printf("%s \n", file);
			cvSaveImage(file, m_pRightPlaneImages[i]);
		}
	}

	return 0;
}

//input point: the origin is at the top-left of image
//nImageIndex: 0-left, 1-right
int CIPanoReg::PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex)
{
	if(m_pLeft==NULL || m_pRight==NULL)
		return -1;

	double radius = (double)(m_pLeft->width) / (2*PI);
	
	double gx,gy,gz;
	SphereTo3D( srcPt.p[0], srcPt.p[1], radius, gx, gy, gz );
	
	//calculate the epipolar constraint
	if(nImageIndex==0) //  lies in the left image
	{
		double leftPt[3];
		leftPt[0] = gx;
		leftPt[1] = gy;
		leftPt[2] = gz;

		double gn = sqrt(gx*gx+gy*gy+gz*gz);
		gx /= gn;
		gy /= gn;
		gz /= gn;

		//locate the image containing the srcPt
		int nTargetImage = 0;
		double fMinAngle = 1000000;
		for(int k=0; k<m_pLeftPlaneImages.size(); k++)
		{
			double vz[3]={0,0,1};
			double camDir[3];
			//calculate the view direction of camera
			double R[9];
			memcpy(R, m_leftPlaneCams[k].R, sizeof(double)*9);
			invers_matrix(R, 3);

			//from plane image space to spherical 3D coordinate
			mult(R, vz, camDir, 3, 3, 1);
			//printf("%lf %lf %lf \n", camDir[0], camDir[1], camDir[2]);

			double nc = sqrt( camDir[0]*camDir[0] + camDir[1]*camDir[1] + camDir[2]*camDir[2] );
			for(int ki=0; ki<3; ki++)
				camDir[ki] /= nc;

			double angle = acos( -(camDir[0]*gx + camDir[1]*gy + camDir[2]*gz) );
			if(fMinAngle>angle)
			{
				fMinAngle = angle;
				nTargetImage = k;
			}
		}		
		printf("Target: %d \n", nTargetImage);
		
		//calculate the point coordinate of perspective image 
		double ipt[3];
		mult( m_leftPlaneCams[nTargetImage].R, leftPt, ipt, 3, 3, 1 );
		double focalLen = m_leftPlaneCams[nTargetImage].focus;
		int sx = ipt[0] / ipt[2] * (-focalLen);
		int sy = ipt[1] / ipt[2] * (-focalLen);
		int ht = m_pLeftPlaneImages[nTargetImage]->height;
		int wd = m_pLeftPlaneImages[nTargetImage]->width;
		sx = sx + wd*0.5;
		sy = ht*0.5 - sy;

		if(1)
		{
			CvPoint pt;
			pt.x = sx;
			pt.y = sy;
			IplImage* pImage = cvCloneImage(m_pLeftPlaneImages[nTargetImage]);
			DrawCross1(pt, 21, pImage);
			cvSaveImage("c:\\temp\\target.jpg", pImage);
			cvReleaseImage(&pImage);
		}
		
		//generate the histogram of srcPt
		const int histDim = 64;
		int srcHist[histDim];
		int nGrayStep = 256 / histDim;
		int blockSize = m_pLeftPlaneImages[nTargetImage]->width * 0.1;
		memset(srcHist, 0, sizeof(int)*histDim);
		GetHistData(sx, sy, m_pLeftPlaneImages[nTargetImage], srcHist, nGrayStep, blockSize);


		//calculate the normal of epipolar plane
		double pn[3];
		mult(m_EM, leftPt, pn, 3, 3, 1);
		
		//calculate the function of epipolar line for each perspective image
		int dstHist[histDim];
		int minDist = 1000000;
		
		for(int k=0; k<m_pRightPlaneImages.size(); k++)
		{
			printf("image: %d \n", k);

			int ht = m_pRightPlaneImages[k]->height;
			int wd = m_pRightPlaneImages[k]->width;
			double pR[9];
			memcpy(pR, m_rightPlaneCams[k].R, sizeof(double)*9);
			invers_matrix(pR, 3);

			/*
			//judge the visibility
			Point3DDouble baselineDir;
			baselineDir.p[0] = m_leftPanoCam.t[0] - m_rightPanoCam.t[0];
			baselineDir.p[1] = m_leftPanoCam.t[1] - m_rightPanoCam.t[1];
			baselineDir.p[2] = m_leftPanoCam.t[2] - m_rightPanoCam.t[2];

			Point3DDouble leftBundleDir;
			leftBundleDir.p[0] = gx;
			leftBundleDir.p[1] = gy;
			leftBundleDir.p[2] = gz;

			Point3DDouble rightBundleDir;
			double vz[3]={0,0,1};
			double camDir[3];
			//from plane image space to spherical 3D coordinate
			mult(pR, vz, camDir, 3, 3, 1);
			camDir[0] = -camDir[0];
			camDir[1] = -camDir[1];
			camDir[2] = -camDir[2];

			//from right pano to the left pano
			double rPanoR[9];
			memcpy(rPanoR, m_rightPanoCam.R, sizeof(double)*9);
			invers_matrix(rPanoR, 3);
			double ct[3];
			mult(rPanoR, camDir, ct, 3, 3, 1);
			rightBundleDir.p[0] = ct[0];
			rightBundleDir.p[1] = ct[1];
			rightBundleDir.p[2] = ct[2];
		    
			double angle1 = angleOfVector(baselineDir, leftBundleDir);
			double angle2 = angleOfVector(baselineDir,   rightBundleDir);
			double angle3 = angleOfVector(leftBundleDir, rightBundleDir);
			if( (angle2-angle1)>0 || (angle3-angle1)>0 )
				continue;
			*/

			double focalLen = m_rightPlaneCams[k].focus;
			
			double lA =   pR[0]*pn[0] + pR[3]*pn[1] + pR[6]*pn[2];
			double lB =   pR[1]*pn[0] + pR[4]*pn[1] + pR[7]*pn[2];
			double lC = -(pR[2]*pn[0] + pR[5]*pn[1] + pR[8]*pn[2])*focalLen;

			//calculate the intersection point
			double x1 = -wd*0.5;
			double y1 = -(lA*x1+lC) / lB; 

			double x2 =  wd*0.5;
			double y2 = -(lA*x2+lC) / lB;

			if(1)
			{
				IplImage* pDisp = cvCloneImage(m_pRightPlaneImages[k]);
				cvLine(pDisp, cvPoint(x1+wd*0.5,ht*0.5-y1), cvPoint(x2+wd*0.5,ht*0.5-y2), CV_RGB(255,0,0));
				char file[256];
				sprintf(file, "c:\\temp\\epipolar_%d.jpg", k);
				cvSaveImage(file, pDisp);
				cvReleaseImage(&pDisp);
			}

			//coarse search along the epipolar line using histogram matching
			//printf("coarse search....\n");
			for(int i=0; i<wd; i++)
			{
				int cx = i  - wd*0.5;
				int cy = (double)(lA*cx+lC) / (double)(-lB+MINIMAL_VALUE);
				//printf("%d %d, ", x,y);

				int x = cx + wd*0.5;
				int y = ht*0.5 - cy;

				if( x<0 || x>=wd || y<0 || y>=ht)
					continue;
								
				memset(dstHist, 0, sizeof(int)*histDim);
				GetHistData(x, y, m_pRightPlaneImages[k], dstHist, nGrayStep, blockSize);

				int ndist = 0;
				for(int ki=0; ki<histDim; ki++)
				{
					ndist = abs( srcHist[ki] - dstHist[ki] );
				}
				if(ndist<minDist)
				{
					minDist = ndist;
					double ipt[3];
					ipt[0] = cx;
					ipt[1] = cy;
					ipt[2] = -focalLen;
					double spt[3];
					mult(pR, ipt, spt, 3, 3, 1);

					double ix,iy;
					GrdToSphere(spt[0], spt[1], spt[2], radius, ix, iy);
					dstPt.p[0] = ix;
					dstPt.p[1] = iy;
				}				
			}
			//printf("\n");
		}
	}
	else if(nImageIndex==1)
	{

	}

	return 0;
}






CIPanoRegTri::CIPanoRegTri()
{
	m_pLeft = NULL;
	m_pRight = NULL;
}
CIPanoRegTri::~CIPanoRegTri()
{
	cvReleaseImage(&m_pLeft);
	cvReleaseImage(&m_pRight);

	for(int i=0; i<m_pLeftPlaneImages.size(); i++)
		cvReleaseImage(&m_pLeftPlaneImages[i]);

	for(int i=0; i<m_pRightPlaneImages.size(); i++)
		cvReleaseImage(&m_pRightPlaneImages[i]);
}

int CIPanoRegTri::Init(char* pLeftFile, char* pRightFile)
{
	//loading gray image
	m_pLeft  = cvLoadImage(pLeftFile,  1);
	m_pRight = cvLoadImage(pRightFile, 1);

	int ht = m_pLeft->height;
	int wd = m_pLeft->width;

	//detect feature points
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	ImgFeature lFeats;
	pFeatDetect->Detect(pLeftFile, lFeats, 480);
	ImgFeature rFeats;
	pFeatDetect->Detect(pRightFile,rFeats, 480);
	delete pFeatDetect;

	//matching
	CMatchBase* pMatch = new CPanoMatch();
	PairMatchRes mr;
	pMatch->Match(lFeats, rFeats, mr);
	int matchNum = mr.matchs.size();
	vector<Point2DDouble> lpts,rpts;

	for(int i=0; i<matchNum; i++)
	{
		int li = mr.matchs[i].l;
		int ri = mr.matchs[i].r;

		/*Point2DDouble lcp = lFeats.GetTopLeftPt(li);
		Point2DDouble rcp = rFeats.GetTopLeftPt(ri);
		lpts_tl.push_back( lcp );
		rpts_tl.push_back( rcp );*/
		
		//if(lcp.p[1]>ht*0.6)
		//	continue;

		//default origin of feature point is at the center of image
		lpts.push_back( lFeats.GetCenteredPt(li) );
		rpts.push_back( rFeats.GetCenteredPt(ri) );
	}
	printf("feature point number: %d \n", lpts.size());

	//relative pose estimation
	CRelativePoseBase* pRP = new CEstimatePose5PointPano();
	//CameraPara leftCam, rightCam;
	m_leftPanoCam.rows  = m_pLeft->height;
	m_leftPanoCam.cols  = m_pLeft->width;
	m_rightPanoCam.rows = m_pRight->height;
	m_rightPanoCam.cols = m_pRight->width;
	pRP->EstimatePose(lpts, rpts, m_leftPanoCam, m_rightPanoCam );  
	printf("inlier points: %d \n", lpts.size());
	
	//GeneratePanoEpipolarImageHeading(m_rightPanoCam.R, m_rightPanoCam.t, m_pLeft, m_pRight);
	
	vector<Point2DDouble> lpts_tl,rpts_tl;
	for(int i=0; i<lpts.size(); i++)
	{
		Point2DDouble pt;
		pt.p[0] = lpts[i].p[0] + wd*0.5;
		pt.p[1] = ht*0.5 - lpts[i].p[1];
		lpts_tl.push_back(pt);

		pt.p[0] = rpts[i].p[0] + wd*0.5;
		pt.p[1] = ht*0.5 - rpts[i].p[1];
		rpts_tl.push_back(pt);
	}
	DrawMatches("c:\\temp\\match.jpg", m_pLeft, m_pRight, lpts_tl, rpts_tl);


	/*
	//calculate the Essential matrix
	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, m_rightPanoCam.R, m_rightPanoCam.t, true, m_EM);
	
	//generate triangles
	m_pTin = new CTINClass("aaa");
	int npt = lpts.size();
	double* px = new double[npt];
	double* py = new double[npt];
	for(int i=0; i<npt; i++)
	{
		px[i] = lpts[i].p[0] + wd*0.5;
		py[i] = ht*0.5 - lpts[i].p[1];
	}
	GenerateDelaunayTriangle(m_pTin, px, py, npt);

	IplImage* pDisp = cvCloneImage(m_pLeft);
	DrawDelaunayTriangle(pDisp, m_pTin);
	cvSaveImage("c:\\temp\\triangles.jpg", pDisp);
	cvReleaseImage(&pDisp);

	m_lpts = lpts;
	m_rpts = rpts;

	//generate perspective images
	m_pLeftPlaneImages.clear();
	m_leftPlaneCams.clear();
	PanoToPlanes(m_pLeft, 60, 90, 90,  1, m_leftPanoCam.R, m_leftPanoCam.t, m_pLeftPlaneImages, m_leftPlaneCams);

	m_pRightPlaneImages.clear();
	m_rightPlaneCams.clear();
	PanoToPlanes(m_pRight, 60, 90, 90, 1, m_rightPanoCam.R, m_rightPanoCam.t, m_pRightPlaneImages, m_rightPlaneCams);

	//for debug
	if(0)
	{
		char file[256];
		for(int i=0; i<m_pLeftPlaneImages.size(); i++)
		{
			sprintf(file, "c:\\temp\\persperctive_left_%d.jpg", i);
			printf("%s \n", file);
			cvSaveImage(file, m_pLeftPlaneImages[i]);

			sprintf(file, "c:\\temp\\persperctive_right_%d.jpg", i);
			printf("%s \n", file);
			cvSaveImage(file, m_pRightPlaneImages[i]);
		}
	}
	*/

	return 0;
}

int CIPanoRegTri::PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex)
{
	if(nImageIndex==0)
	{
		int ht = m_pLeft->height;
		int wd = m_pLeft->width;

	

		//############  coarse searching based on triangles ###############
		dstPt.p[0] = 0;
		dstPt.p[1] = 0;
		//find the triangle containing the srcpt
		int pIndex[3];
		int res = SelectTriangle(srcPt.p[0], srcPt.p[1], m_pTin, pIndex);

		if(res>0)
		{
			Point2DDouble lp[3];
			Point2DDouble rp[3];
			for(int i=0; i<3; i++)
			{
				lp[i] = m_lpts[pIndex[i]];
				rp[i] = m_rpts[pIndex[i]];

				//from center to top-left
				lp[i].p[0] = lp[i].p[0] + wd*0.5;
				lp[i].p[1] = ht*0.5 - lp[i].p[1];
				rp[i].p[0] = rp[i].p[0] + wd*0.5;
				rp[i].p[1] = ht*0.5 - rp[i].p[1];
			}

			double dis[3];
			double sumdis = 0;
			for(int i=0; i<3; i++)
			{
				dis[i] = distanceVec(srcPt, lp[i]);
				sumdis += dis[i];
			}
			for(int i=0; i<3; i++)
			{
				dis[i] /= sumdis;
			}
			for(int i=0; i<3; i++)
			{
				dstPt.p[0] += rp[i].p[0]*dis[i];
				dstPt.p[1] += rp[i].p[1]*dis[i];
			}
		}
		else
		{

		}

		//#################  precise searching #######################



	}
	return 0;
}




CIPanoRegDirect::CIPanoRegDirect()
{

}

CIPanoRegDirect::~CIPanoRegDirect()
{

}

int CIPanoRegDirect::Init(char* pLeftFile, char* pRightFile)
{

	//loading gray image
	m_pLeft = cvLoadImage(pLeftFile, 0);
	m_pRight = cvLoadImage(pRightFile, 0);

	//detect feature points
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	ImgFeature lFeats;
	pFeatDetect->Detect(pLeftFile, lFeats, 480);
	ImgFeature rFeats;
	pFeatDetect->Detect(pRightFile,rFeats, 480);
	delete pFeatDetect;

	//matching
	CMatchBase* pMatch = new CPanoMatch();
	PairMatchRes mr;
	pMatch->Match(lFeats, rFeats, mr);
	int matchNum = mr.matchs.size();
	vector<Point2DDouble> lpts,rpts;
	lpts.resize(matchNum);
	rpts.resize(matchNum);
	for(int i=0; i<matchNum; i++)
	{

		int li = mr.matchs[i].l;
		int ri = mr.matchs[i].r;

		//default origin of feature point is at the center of image
		lpts[i] = lFeats.GetCenteredPt(li);
		rpts[i] = rFeats.GetCenteredPt(ri);
	}

	//relative pose estimation
	CRelativePoseBase* pRP = new CEstimatePose5PointPano();
	//CameraPara leftCam, rightCam;
	m_leftPanoCam.rows  = m_pLeft->height;
	m_leftPanoCam.cols  = m_pLeft->width;
	m_rightPanoCam.rows = m_pRight->height;
	m_rightPanoCam.cols = m_pRight->width;
	pRP->EstimatePose(lpts, rpts, m_leftPanoCam, m_rightPanoCam );  


	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, m_rightPanoCam.R, m_rightPanoCam.t, true, m_EM);


	return 0;
}

int CIPanoRegDirect::PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex)
{

	//generate the normal of epipolar plane



	return 0;
}