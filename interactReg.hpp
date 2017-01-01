
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
class CIRegBase
{
public:
	CIRegBase();
	~CIRegBase();

	virtual int PtReg(IplImage* pLeft, IplImage* pRight, Point2DDouble srcPt, Point2DDouble& dstPt);
};


//for interactive panorama image registration
class CIPanoReg: public CIRegBase
{
public:
	CIPanoReg();
	~CIPanoReg();

	int PtReg(IplImage* pLeft, IplImage* pRight, Point2DDouble srcPt, Point2DDouble& dstPt);

};



#endif
