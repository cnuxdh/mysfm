

#ifndef	TRIANGLES_H
#define TRIANGLES_H

#include "export.hpp"


#include<vector>
using namespace std;


typedef struct STRUCT_TRIANGLE
{
	int  p1,p2,p3;
}stTriangle;

typedef struct STRUCT_VERTEXT
{
	double x,y;
}stVertex;


class DLL_EXPORT CTrianglesBase
{
public:
	CTrianglesBase(){}
	virtual ~CTrianglesBase(){}

	virtual int Generate( vector<stVertex>& pts, vector<stTriangle>& triangles){ return 0; }
};

//algorithm from CMU: http://www.cs.cmu.edu/~quake/triangle.html
class DLL_EXPORT  CTrianglesCMU: public CTrianglesBase
{
public:
	CTrianglesCMU();
	~CTrianglesCMU();

	int Generate( vector<stVertex>& pts, vector<stTriangle>& triangles);
};

//algorithm using TinLib: http://jca3.freeshell.org/tinlib/tinlib.html
class DLL_EXPORT  CTinLib: public CTrianglesBase
{
public:
	CTinLib();
	~CTinLib();

	int Generate( vector<stVertex>& pts, vector<stTriangle>& triangles);
};



#endif