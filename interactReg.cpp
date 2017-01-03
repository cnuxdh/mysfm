
#include "interactReg.hpp"

#include "register.hpp"
#include "relativepose.hpp"
#include "sift.hpp"
#include "panorama.hpp"


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
	pFeatDetect->Detect(pLeftFile, lFeats, 640);
	ImgFeature rFeats;
	pFeatDetect->Detect(pRightFile,rFeats, 640);
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

		lpts[i] = lFeats.GetCenteredPt(li);
		rpts[i] = rFeats.GetCenteredPt(ri);
	}

	//relative pose estimation
	CRelativePoseBase* pRP = new CEstimatePose5PointPano();
	CameraPara leftCam, rightCam;
	leftCam.rows = m_pLeft->height;
	leftCam.cols = m_pLeft->width;
	rightCam.rows = m_pRight->height;
	rightCam.cols = m_pRight->width;
	pRP->EstimatePose(lpts, rpts, leftCam, rightCam );  

	//generate perspective images
	PanoToPlanes(m_pLeft, 60, 90, 90, 1, leftCam.R, leftCam.t, m_pLeftPlaneImages, m_leftPlaneCams);
	if(1)
	{
		char file[256];
		for(int i=0; i<m_leftPlaneCams.size(); i++)
		{
			sprintf(file, "c:\\temp\\persperctive_%d.jpg", i);
			cvSaveImage(file, m_pLeftPlaneImages[i]);
		}
	}

	return 0;
}

int CIPanoReg::PtReg(Point2DDouble srcPt, Point2DDouble& dstPt)
{
	//split the spherical image to several perspective images




	return 0;
}