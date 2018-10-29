#ifndef PANO_REG_USING_POS
#define PANO_REG_USING_POS

#include "defines.hpp"


//for interactive panorama image registration
class DLL_EXPORT CIPanoRegUsingPos
{
public:
	CIPanoRegUsingPos();
	~CIPanoRegUsingPos();

	//int Init(char* pLeftFile, char* pRightFile);
	//int Init(IplImage* pLeft, IplImage* pRight, CameraPara leftCam, CameraPara rightCam);

	int Init(CameraPara leftCam, CameraPara rightCam);
	
	//find the initial ground point according to the input point based on the POS 
	int GroundPointReg(Point2DDouble leftPt, Point2DDouble& rightPt);
	
	//int GetEpipolarLinePts(Point2DDouble srcPt, int nImageIndex, vector<Point2DDouble>& epts);
	//int MatchingPoint(Point2DDouble leftPt, Point2DDouble& rightPt);
	//int PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex);
	
private:

	//for panorama images
	//IplImage*  m_pLeft;
	//IplImage*  m_pRight;
	CameraPara m_leftPanoCam, m_rightPanoCam;
	double     m_EM[9]; //essential matrix
	double     m_R[9];
	double     m_T[3];
};

#endif
