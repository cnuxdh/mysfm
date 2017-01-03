
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




/* base class for image interactive registration
*/
class DLL_EXPORT CIRegBase
{
public:
	CIRegBase(){}
	virtual ~CIRegBase(){}

	virtual int Init(char* pLeftFile, char* pRightFile){return 0;}
	virtual int PtReg(IplImage* pLeft, IplImage* pRight, Point2DDouble srcPt, Point2DDouble& dstPt){return 0;}
};


//for interactive panorama image registration
class DLL_EXPORT CIPanoReg: public CIRegBase
{
public:
	CIPanoReg();
	~CIPanoReg();

	int Init(char* pLeftFile, char* pRightFile);
	int PtReg(Point2DDouble srcPt, Point2DDouble& dstPt);

private:
	IplImage* m_pLeft;
	IplImage* m_pRight;
	vector<IplImage*>  m_pLeftPlaneImages;
	vector<CameraPara> m_leftPlaneCams;
	vector<IplImage*>  m_pRightPlaneImages;
	vector<CameraPara> m_rightPlaneCams;
};


#endif
