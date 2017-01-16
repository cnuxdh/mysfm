
#ifndef INTERACTIVE_REG_HPP
#define INTERACTIVE_REG_HPP

#include "defines.hpp"

#ifdef OPENCV_1X 
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#else
//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
//#include "opencv2/calib3d/calib3d.hpp"
#endif


//corelib 
#include "Triangles.h"

//#include "pos.hpp"



/* base class for image interactive registration
*/
class DLL_EXPORT CIRegBase
{
public:
	CIRegBase(){}
	virtual ~CIRegBase(){}

	virtual int Init(char* pLeftFile, char* pRightFile){return 0;}
	virtual int Init(IplImage* pLeft, IplImage* pRight, CameraPara leftCam, CameraPara rightCam){return 0;}
	virtual int PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex){return 0;}
	virtual int GetEpipolarLinePts(Point2DDouble srcPt, int nImageIndex, vector<Point2DDouble>& epts){return 0;}
};


//for interactive panorama image registration
class DLL_EXPORT CIPanoReg: public CIRegBase
{
public:
	CIPanoReg();
	~CIPanoReg();

	int Init(char* pLeftFile, char* pRightFile);
	int Init(IplImage* pLeft, IplImage* pRight, CameraPara leftCam, CameraPara rightCam){return 0;}
	int PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex);
	int GetEpipolarLinePts(Point2DDouble srcPt, int nImageIndex, vector<Point2DDouble>& epts){return 0;}

private:

	//for panorama images
	IplImage*  m_pLeft;
	IplImage*  m_pRight;
	CameraPara m_leftPanoCam, m_rightPanoCam;
	double     m_EM[9]; //essential matrix
	double     m_R[9];
	double     m_T[3];

	//for perspective images
	vector<IplImage*>  m_pLeftPlaneImages;
	vector<CameraPara> m_leftPlaneCams;
	vector<IplImage*>  m_pRightPlaneImages;
	vector<CameraPara> m_rightPlaneCams;
};


//for interactive panorama image registration
class DLL_EXPORT CIPanoRegTri: public CIRegBase
{
public:
	CIPanoRegTri();
	~CIPanoRegTri();

	int Init(char* pLeftFile, char* pRightFile);
	int Init(IplImage* pLeft, IplImage* pRight, CameraPara leftCam, CameraPara rightCam){return 0;}
	int PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex);
	int GetEpipolarLinePts(Point2DDouble srcPt, int nImageIndex, vector<Point2DDouble>& epts){return 0;}

private:

	//for panorama images
	IplImage*  m_pLeft;
	IplImage*  m_pRight;
	CameraPara m_leftPanoCam, m_rightPanoCam;
	double     m_EM[9]; //essential matrix
	double     m_R[9];
	double     m_T[3];

	//feature points
	vector<Point2DDouble> m_lpts,m_rpts;

	//for perspective images
	vector<IplImage*>  m_pLeftPlaneImages;
	vector<CameraPara> m_leftPlaneCams;
	vector<IplImage*>  m_pRightPlaneImages;
	vector<CameraPara> m_rightPlaneCams;

	//for triangles
	CTINClass* m_pTin;
};



//for interactive panorama image registration
class DLL_EXPORT CIPanoRegDirect: public CIRegBase
{
public:
	CIPanoRegDirect();
	~CIPanoRegDirect();

	int Init(char* pLeftFile, char* pRightFile);
	int Init(IplImage* pLeft, IplImage* pRight, CameraPara leftCam, CameraPara rightCam);
	int PtReg(Point2DDouble srcPt, Point2DDouble& dstPt, int nImageIndex);
	int GetEpipolarLinePts(Point2DDouble srcPt, int nImageIndex, vector<Point2DDouble>& epts);

private:

	//for panorama images
	IplImage*  m_pLeft;
	IplImage*  m_pRight;
	CameraPara m_leftPanoCam, m_rightPanoCam;
	double     m_EM[9]; //essential matrix
	double     m_R[9];
	double     m_T[3];
};


#endif
