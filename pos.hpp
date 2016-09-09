

#ifndef POS_H
#define POS_H

#include "export.hpp"

#include <vector>
using namespace std;

typedef struct stPOSInfo
{
	double lat,lon;
	double height;
	double rx,ry,rz;
	int    id;
}POSInfo;


class DLL_EXPORT CReadPosBase
{
public:
	CReadPosBase(){}
	virtual ~CReadPosBase(){}
	
	virtual int ReadPOSData(char* filename){return 0;}
};


//read the pos data of UAV
class DLL_EXPORT CReadUAVPose: public CReadPosBase
{
public:
	CReadUAVPose();
	~CReadUAVPose();
  
	int ReadPOSData(char* filename);

private:
	vector<POSInfo> m_posArray;
};


#endif