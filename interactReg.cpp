
#include "interactReg.hpp"

#include "register.hpp"
#include "relativepose.hpp"
#include "sift.hpp"
#include "panorama.hpp"
#include "funcs.hpp"
#include "hog.hpp"


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

	return 0;
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
	m_pLeft  = cvLoadImage(pLeftFile, 0);
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


	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, m_rightPanoCam.R, m_rightPanoCam.t, true, 
		m_R, m_T, m_EM);
		
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
	m_pLeft = NULL;
	m_pRight = NULL;
	
}

CIPanoRegDirect::~CIPanoRegDirect()
{

}

int CIPanoRegDirect::Init(char* pLeftFile, char* pRightFile)
{

	//loading gray image
	m_pLeft = cvLoadImage(pLeftFile);
	m_pRight = cvLoadImage(pRightFile);

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


	//calculate the essential matrix
	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, m_rightPanoCam.R, m_rightPanoCam.t, 
		pRP->IsExplicit(), m_R,m_T, m_EM);
	
	return 0;
}


//for pos
int CIPanoRegDirect::Init(IplImage* pLeft, IplImage* pRight, CameraPara leftCam, CameraPara rightCam)
{
	if(m_pLeft!=NULL)
		cvReleaseImage(&m_pLeft);
	m_pLeft = cvCloneImage(pLeft);

	if(m_pRight!=NULL)
		cvReleaseImage(&m_pRight);
	m_pRight = cvCloneImage(pRight);

	m_leftPanoCam.rows  = m_pLeft->height;
	m_leftPanoCam.cols  = m_pLeft->width;
	m_rightPanoCam.rows = m_pRight->height;
	m_rightPanoCam.cols = m_pRight->width;

	m_leftPanoCam  = leftCam;
	m_rightPanoCam = rightCam;

	//calculate the essential matrix
	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.t, m_rightPanoCam.R, m_rightPanoCam.t, m_leftPanoCam.bIsExplicit, 
		m_R, m_T, m_EM);
	
	//generate depth image

	return 0;
}

int CIPanoRegDirect::MatchingPoint(Point2DDouble srcPt, Point2DDouble& dstPt)
{
	int ht = m_pLeft->height;
	int wd = m_pLeft->width;
	double cx = srcPt.p[0] - wd*0.5;
	double cy = ht*0.5 - srcPt.p[1];

	//from pano 2D to sphere 3D 
	double radius = (double)(m_pLeft->width) / (2*PI);
	double lp[3];
	SphereTo3D_center(cx, cy, radius, lp[0], lp[1], lp[2]);

	//generate perspective image
	IplImage* lPtPlane = PanoToPlane(m_pLeft, lp, 12, 12);
	cvSaveImage("c:\\temp\\left_pt_plane.jpg", lPtPlane);
	vector<double> lHist;
	//CalculateColorHist(lPtPlane, 8, lHist);
	CalculateHOG(lPtPlane, lHist);
	cvReleaseImage(&lPtPlane);
	
	Point3DDouble leftBundleDir;
	leftBundleDir.p[0] = lp[0];
	leftBundleDir.p[1] = lp[1];
	leftBundleDir.p[2] = lp[2];

	double lt[3];
	double iR[9];
	memcpy(iR, m_R, sizeof(double)*9);
	invers_matrix(iR, 3);
	mult(iR, m_T, lt, 3, 3, 1);
	Point3DDouble baselineDir;
	baselineDir.p[0] = lt[0];
	baselineDir.p[1] = lt[1];
	baselineDir.p[2] = lt[2];
	//baselineDir.p[0] = m_leftPanoCam.t[0] - m_rightPanoCam.t[0];
	//baselineDir.p[1] = m_leftPanoCam.t[1] - m_rightPanoCam.t[1];
	//baselineDir.p[2] = m_leftPanoCam.t[2] - m_rightPanoCam.t[2];
	//printf("baseline : %lf %lf %lf \n", baselineDir.p[0], baselineDir.p[1], baselineDir.p[2]);

	//calculate the epipolar normal
	double n[3];
	mult(m_EM, lp, n, 3, 3, 1);
	Point3DDouble enormal;
	enormal.p[0] = n[0];
	enormal.p[1] = n[1];
	enormal.p[2] = n[2];

	//generate the epipolar points
	vector<Point3DDouble> ptVecs =  GenerateEpipolarPlaneVectors(enormal, 360);


	//double rPanoR[9];
	//memcpy(rPanoR, m_R, sizeof(double)*9);
	//invers_matrix(rPanoR, 3);

	radius = (double)(m_pRight->width) / (2*PI);
	IplImage* pDisp = cvCloneImage(m_pRight);
	double minDif = 100000000;
	double maxSim = -100000000;
	int    nIndex = 0;
	vector<double> vecSim;
	for(int i=0; i<ptVecs.size(); i++)
	{
		//visibility decision							
		//from right pano to the left pano
		double rp[3];
		rp[0] = ptVecs[i].p[0];
		rp[1] = ptVecs[i].p[1];
		rp[2] = ptVecs[i].p[2];

		double ct[3];
		mult(iR, rp, ct, 3, 3, 1);

		Point3DDouble rightBundleDir;
		rightBundleDir.p[0] = ct[0];
		rightBundleDir.p[1] = ct[1];
		rightBundleDir.p[2] = ct[2];				

		//printf("%lf %lf %lf \n", )

		Point3DDouble rightTop,rightBottom;
		CalculateViewField(baselineDir, leftBundleDir, 100, 2, rightBottom, rightTop);

		double angle1,angle2,angle3;
		if(0)
		{
			angle1 = angleOfVector(baselineDir,   leftBundleDir);
			angle2 = angleOfVector(baselineDir,   rightBundleDir);
			angle3 = angleOfVector(leftBundleDir, rightBundleDir);
			//printf("%lf %lf %lf  %lf \n", angle1, angle2, angle3, angle1-angle2-angle3);
		}
		else
		{
			angle1 = angleOfVector(rightBottom,   rightTop);
			angle2 = angleOfVector(rightBottom,   rightBundleDir);
			angle3 = angleOfVector(rightTop,      rightBundleDir);
		}

		double ix,iy;
		GrdToSphere_center(ptVecs[i].p[0], ptVecs[i].p[1], ptVecs[i].p[2], radius, ix, iy);
		CvPoint ip;
		ip.x = ix+wd*0.5;
		ip.y = ht*0.5-iy;
		cvDrawCircle(pDisp, ip, 1, CV_RGB(255,0,0), 2);


		if( fabs(angle1-angle2-angle3)>5 )
			continue;

		cvDrawCircle(pDisp, ip, 1, CV_RGB(0,255,0), 2);

		//
		IplImage* rPtPlane = PanoToPlane(m_pRight, rp, 12, 12);
		//cvSaveImage("c:\\temp\\left_pt_plane.jpg", lPtPlane);
		vector<double> rHist;
		//CalculateColorHist(rPtPlane, 8, rHist);
		CalculateHOG(rPtPlane, rHist);


		//calculate the similarity
		double sim = Cov(lHist, rHist);
		vecSim.push_back(sim);

		if(0)
		{
			char filename[256];
			sprintf(filename, "c:\\temp\\patch_%d.jpg", vecSim.size());
			cvSaveImage(filename, rPtPlane);
		}

		/*
		for(int ki=0; ki<lHist.size(); ki++)
		{
		sim += fabs( lHist[ki]-rHist[ki] );
		}
		if(sim<minDif)
		{
		minDif = sim;
		dstPt.p[0] = ip.x;
		dstPt.p[1] = ip.y;
		//cvSaveImage("c:\\temp\\dstPtImg.jpg", rPtPlane);
		}
		*/
		if(sim>maxSim)
		{
			maxSim = sim;
			dstPt.p[0] = ip.x;
			dstPt.p[1] = ip.y;
			nIndex = vecSim.size();

			if(1)
			{
				cvSaveImage("c:\\temp\\dstPtImg.jpg", rPtPlane);
			}
		}

		cvReleaseImage(&rPtPlane);


		/*
		GrdToSphere_center(ptVecs[i+1].p[0], ptVecs[i+1].p[1], ptVecs[i+1].p[2], radius, ix, iy);
		CvPoint ip2;
		ip2.x = ix+wd*0.5;
		ip2.y = ht*0.5-iy;
		//cvDrawCircle(pDisp, ip, 1, CV_RGB(255,0,0), 1);
		cvLine(pDisp, ip1, ip2, CV_RGB(0,255,0),1);
		*/
	}

	if(1)
	{
		FILE* fp = fopen("c:\\temp\\sim.txt", "w");
		for(int k=0; k<vecSim.size(); k++)
		{
			fprintf(fp, "%lf \n", vecSim[k]);
		}
		fclose(fp);
	}

	printf("best: %d \n", nIndex);

	cvSaveImage("c:\\temp\\epipolarLine.jpg", pDisp);
	cvReleaseImage(&pDisp);

	return 0;
}

int CIPanoRegDirect::PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex)
{
	//generate the normal of epipolar plane
	double n[3];

	if(nImageIndex==0)
	{
		int ht = m_pLeft->height;
		int wd = m_pLeft->width;
		double cx = srcPt.p[0] - wd*0.5;
		double cy = ht*0.5 - srcPt.p[1];

		//from pano 2D to sphere 3D 
		double radius = (double)(m_pLeft->width) / (2*PI);
		double lp[3];
		SphereTo3D_center(cx, cy, radius, lp[0], lp[1], lp[2]);
		//printf("left bundle: %lf %lf %lf \n", lp[0], lp[1], lp[2]);


		if(srcPt.nType == 0) //for ground point
		{
			GroundPointReg(srcPt, dstPt);
		}
		else
		{
			//calculate the distance

		}
	}

	return 0;
}


int CIPanoRegDirect::GetEpipolarLinePts(Point2DDouble srcPt, int nImageIndex, vector<Point2DDouble>& epts)
{
	double n[3];

	epts.clear();

	if(nImageIndex==0)
	{
		int ht = m_pLeft->height;
		int wd = m_pLeft->width;
		double cx = srcPt.p[0] - wd*0.5;
		double cy = ht*0.5 - srcPt.p[1];

		//from pano 2D to sphere 3D 
		double radius = (double)(m_pLeft->width) / (2*PI);
		double lp[3];
		SphereTo3D_center(cx, cy, radius, lp[0], lp[1], lp[2]);

		Point3DDouble leftBundleDir;
		leftBundleDir.p[0] = lp[0];
		leftBundleDir.p[1] = lp[1];
		leftBundleDir.p[2] = lp[2];

		double lt[3];
		double iR[9];
		memcpy(iR, m_R, sizeof(double)*9);
		invers_matrix(iR, 3);
		mult(iR, m_T, lt, 3, 3, 1);
		Point3DDouble baselineDir;
		baselineDir.p[0] = lt[0];
		baselineDir.p[1] = lt[1];
		baselineDir.p[2] = lt[2];

		//calculate the epipolar normal
		mult(m_EM, lp, n, 3, 3, 1);
		Point3DDouble enormal;
		enormal.p[0] = n[0];
		enormal.p[1] = n[1];
		enormal.p[2] = n[2];

		//generate the points
		vector<Point3DDouble> ptVecs =  GenerateEpipolarPlaneVectors(enormal, 360);

		for(int i=0; i<ptVecs.size(); i++)
		{
			//visibility decision							
			//from right pano to the left pano
			double rp[3];
			rp[0] = ptVecs[i].p[0];
			rp[1] = ptVecs[i].p[1];
			rp[2] = ptVecs[i].p[2];

			double ct[3];
			mult(iR, rp, ct, 3, 3, 1);

			Point3DDouble rightBundleDir;
			rightBundleDir.p[0] = ct[0];
			rightBundleDir.p[1] = ct[1];
			rightBundleDir.p[2] = ct[2];				

			double angle1,angle2,angle3;
			angle1 = angleOfVector(baselineDir,   leftBundleDir);
			angle2 = angleOfVector(baselineDir,   rightBundleDir);
			angle3 = angleOfVector(leftBundleDir, rightBundleDir);
			
			if( fabs(angle1-angle2-angle3)>5 )
				continue;
			
			radius = (double)(m_pRight->width) / (2*PI);
			double ix,iy;
			GrdToSphere_center(ptVecs[i].p[0], ptVecs[i].p[1], ptVecs[i].p[2], radius, ix, iy);
			Point2DDouble ip;
			ip.p[0] = ix+wd*0.5;
			ip.p[1] = ht*0.5-iy;
			epts.push_back(ip);
		}
	}

	//order
	for(int j=0; j<epts.size(); j++)
		for(int i=j+1; i<epts.size(); i++)
		{
			if(epts[j].p[0]>epts[i].p[0])
			{
				Point2DDouble tp;
				tp = epts[j];
				epts[j]=epts[i];
				epts[i]=tp;
			}
		}
		
	return 0;
}

/* inputs:
		leftpt,rightpt: top-left image coordinate
*/
int CIPanoRegDirect::GroundPointReg(Point2DDouble leftPt, Point2DDouble& rightPt)
{
	int lwd = m_leftPanoCam.cols;
	int lht = m_leftPanoCam.rows;
	int relativeHei = m_leftPanoCam.camRelativeHei;

	//calculate the angle
	double angle = (leftPt.p[1] - lht*0.5) / (double)(lht) * PI ;

	//calculate the distance
	double distance = relativeHei / sin(angle);

	//calculate the 3D coordinate of ground point
	double sx,sy,sz;
	double cx = leftPt.p[0] - lwd*0.5;
	double cy = lht*0.5 - leftPt.p[1];
	double lradius = (double)(m_leftPanoCam.cols) / (2*PI);
	SphereTo3D_center(cx, cy, lradius, sx, sy, sz);
	double ng = sqrt(sx*sx + sy*sy + sz*sz);
	sx = sx / ng;
	sy = sy / ng;
	sz = sz / ng;
	
	//calculate the ground point
	double rp[3];
	rp[0] = sx*distance;
	rp[1] = sy*distance;
	rp[2] = sz*distance;
	double lR[9];
	memcpy(lR, m_leftPanoCam.R, sizeof(double)*9);
	invers_matrix(lR, 3);
	double res[3];
	mult(lR, rp, res, 3, 3, 1);
	double gx,gy,gz;
	gx = m_leftPanoCam.t[0] + res[0];
	gy = m_leftPanoCam.t[1] + res[1];
	gz = m_leftPanoCam.t[2] + res[2];

	//transform to the right camera
	double gp[3];
	gp[0] = gx - m_rightPanoCam.t[0];
	gp[1] = gy - m_rightPanoCam.t[1];
	gp[2] = gz - m_rightPanoCam.t[2];
	mult(m_rightPanoCam.R, gp, res, 3, 3, 1);

	double radius = (double)(m_rightPanoCam.cols) / (2*PI);
	GrdToSphere_center(res[0], res[1], res[2], radius, cx, cy);

	rightPt.p[0] = cx + m_rightPanoCam.cols*0.5;
	rightPt.p[1] = m_rightPanoCam.rows*0.5 - cy;

	return 0;
}
