
#ifndef BUNDLER_IO
#define BUNDLER_IO

#include "export.hpp"

//corelib
#include "corelib/commondata.h"

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





#endif