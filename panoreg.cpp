
#include "panoreg.hpp"

#include "panorama.hpp"

//corelib
#include "Corelib/Matrix.h"
#include "CommonFuncs.h"

CIPanoRegUsingPos::CIPanoRegUsingPos()
{

}
CIPanoRegUsingPos::~CIPanoRegUsingPos()
{

}

int CIPanoRegUsingPos::Init(CameraPara leftCam, CameraPara rightCam)
{
	m_leftPanoCam  = leftCam;
	m_rightPanoCam = rightCam;

	//calculate the essential matrix
	CalculateEssentialMatrix(m_leftPanoCam.R, m_leftPanoCam.T, m_rightPanoCam.R, m_rightPanoCam.T, m_leftPanoCam.bIsExplicit, 
		m_R, m_T, m_EM);

	return 0;
}

/* inputs:
		leftpt,rightpt: top-left image coordinate
*/
int CIPanoRegUsingPos::GroundPointReg(Point2DDouble leftPt, Point2DDouble& rightPt)
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
	gx = m_leftPanoCam.T[0] + res[0];
	gy = m_leftPanoCam.T[1] + res[1];
	gz = m_leftPanoCam.T[2] + res[2];

	//transform to the right camera
	double gp[3];
	gp[0] = gx - m_rightPanoCam.T[0];
	gp[1] = gy - m_rightPanoCam.T[1];
	gp[2] = gz - m_rightPanoCam.T[2];
	mult(m_rightPanoCam.R, gp, res, 3, 3, 1);

	double radius = (double)(m_rightPanoCam.cols) / (2*PI);
	GrdToSphere_center(res[0], res[1], res[2], radius, cx, cy);

	rightPt.p[0] = cx + m_rightPanoCam.cols*0.5;
	rightPt.p[1] = m_rightPanoCam.rows*0.5 - cy;

	return 0;
}
