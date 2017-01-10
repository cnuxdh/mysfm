

#ifndef POS_H
#define POS_H

#include "export.hpp"

#include <vector>
using namespace std;

typedef struct stPOSInfo
{
	double lat,lon;
	double height;
	double roll,pitch,yaw; //Eular angle
	int    id;
	double gx,gy;          //UTM x,y
	char   filetitle[256]; //file title
}POSInfo;


class DLL_EXPORT CReadPosBase
{
public:
	CReadPosBase(){}
	virtual ~CReadPosBase(){}
	
	virtual int ReadPOSData(char* filename){return 0;}
	virtual int GetPOS(int index, POSInfo& pos){return 0;}
};


//read the pos data of UAV
class DLL_EXPORT CReadUAVPose: public CReadPosBase
{
public:
	CReadUAVPose();
	~CReadUAVPose();
  
	int ReadPOSData(char* filename);
	int GetPOS(int index, POSInfo& pos);

private:
	vector<POSInfo> m_posArray;
};


//read the pos data installed in the car
class DLL_EXPORT CReadCarPos: public CReadPosBase
{
public:
	CReadCarPos();
	~CReadCarPos();

	int ReadPOSData(char* filename);
	int GetPOS(int index, POSInfo& pos);

private:
	vector<POSInfo> m_posArray;
};

#endif