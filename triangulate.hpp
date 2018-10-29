
#ifndef TRIANGULATE_HPP
#define TRIANGULATE_HPP


#include "defines.hpp"


DLL_EXPORT Point3DDouble TriangulatePt(Point2DDouble p, Point2DDouble q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error);


//////////////////////////////////////////////////////////////////////////
//recover 3D coordinates of point
class DLL_EXPORT CTriangulateBase
{
public:
	CTriangulateBase(){}
	virtual ~CTriangulateBase(){}

	//for track sequence with two projections
	virtual void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
		CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps){} 

	virtual void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
		CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, vector<double>& errorArray){} 

	//for one track with multiple projections
	virtual void Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, 
		Point3DDouble& gps,bool explicit_camera_centers,double& ferror){} 

};

class DLL_EXPORT CTriangulateCV: public CTriangulateBase
{
public:
	CTriangulateCV();
	~CTriangulateCV();
	void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts,
		CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps);
	void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
		CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, vector<double>& errorArray);
	void Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, 
		Point3DDouble& gps,bool explicit_camera_centers,double& ferror);
};

class DLL_EXPORT CTriangulatePano: public CTriangulateBase
{
public:
	CTriangulatePano();
	~CTriangulatePano();
	void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts,
		CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps);
	void Triangulate(vector<Point2DDouble> lPts, vector<Point2DDouble> rPts, 
		CameraPara cam1, CameraPara cam2, vector<Point3DDouble>& gps, vector<double>& errorArray);
	void Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, 
		Point3DDouble& gps,bool explicit_camera_centers,double& ferror);
};

//////////////////////////////////////////////////////////////////////////




#endif