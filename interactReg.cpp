
#include "interactReg.hpp"

#include "register.hpp"
#include "relativepose.hpp"
#include "sift.hpp"
#include "panorama.hpp"

#include "CommonFuncs.h"
#include "matrix.h"


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

		//calculate the normal of epipolar plane
		double pn[3];
		mult(m_EM, leftPt, pn, 3, 3, 1);

		for(int k=0; k<m_pRightPlaneImages.size(); k++)
		{
			int ht = m_pRightPlaneImages[k]->height;
			int wd = m_pRightPlaneImages[k]->width;
			
			double pR[9];
			memcpy(pR, m_rightPlaneCams[k].R, sizeof(double)*9);
			invers_matrix(pR, 3);

			double focalLen = m_rightPlaneCams[k].focus;
			
			//calculate the function of epipolar line
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


		}
	}
	else if(nImageIndex==1)
	{

	}

	return 0;
}