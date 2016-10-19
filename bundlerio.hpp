
#ifndef BUNDLER_IO
#define BUNDLER_IO

#include "export.hpp"
#include "defines.hpp"

//corelib
//#include "commondata.h"

//coredll
#include "defs.h"

#include <vector>
using namespace std;




DLL_EXPORT int ReadBundlerOutFile(char* filename, vector<stPOS>& camParas, vector<stTrack>& tracks );



DLL_EXPORT int ReadPMVSPly(char* filename, vector<stTrack>& tracks);
DLL_EXPORT int WritePMVSPly(char* filename, vector<stTrack>& tracks);


DLL_EXPORT int ReadPMVSPly(char* filename, stTrack** tracks, int* nTrack);
DLL_EXPORT int WritePMVSPly(char* filename, stTrack* tracks, int nTrack);


DLL_EXPORT int WritePMVSPly(char* filename, 
	double* px, double* py,  double* pz, int nPt);



DLL_EXPORT int WritePMVSPly(char* filename, const vector<Point3DDouble>& gpts);




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