
#ifndef BUNDLER_IO
#define BUNDLER_IO

#include "export.hpp"
#include "defines.hpp"



//corelib
#include "commondata.h"


//coredll
#include "defs.h"


#include <vector>
using namespace std;

typedef struct STRU_TRACK
{
	double x,y,z;
	unsigned char r,g,b;
	vector<POINT2> imgpt;
	vector<int> imgid;
	vector<int> ptid;
}stTrack;



DLL_EXPORT int ReadBundlerOutFile(char* filename, vector<stPOS>& camParas, vector<stTrack>& tracks );




//////////////////////////////////////////////////////////////////////////
//output 3D points into model file
class DLL_EXPORT CModelFileBase
{
public:
	CModelFileBase(){}
	virtual ~CModelFileBase(){}
	virtual int Save(char* modelFile, vector<Point3DDouble> pts){return 0;}
	virtual int Save(char* modelFile, vector<Point3DDouble> pts, vector<Point3DDouble> colors){return 0;}
};

class DLL_EXPORT CPlyModel: public CModelFileBase
{
public:
	CPlyModel();
	~CPlyModel();
	int Save(char* modelFile, vector<Point3DDouble> pts);
	int Save(char* modelFile, vector<Point3DDouble> pts, vector<Point3DDouble> colors);
};





#endif