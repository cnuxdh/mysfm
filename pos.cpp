
#include "pos.hpp"


//corelib
#include "CommonFuncs.h"


//##############################################################
CReadUAVPose::CReadUAVPose()
{

}

CReadUAVPose::~CReadUAVPose()
{

}

int CReadUAVPose::ReadPOSData(char* filename)
{

	return 0;
}
int CReadUAVPose::GetPOS(int index, POSInfo& pos)
{

	return 0;
}
int CReadUAVPose::GetPOSNum()
{
	return m_posArray.size();
}


//##################################################################
CReadCarPos::CReadCarPos()
{

}
CReadCarPos::~CReadCarPos()
{

}

int CReadCarPos::GetPOSNum()
{
	return m_posArray.size();
}

int CReadCarPos::ReadPOSData(char* filepath)
{
	char sline[512];

	int rows = GetFileRows(filepath);
	
	m_posArray.resize(rows-1);

	char   filename[256];
	double lat=0,lon=0,hei=0;
	double gx=0,gy=0;
	double roll, pitch, yaw; 
	double dtime;
	char   day[256];
	char   stime[256];

	FILE* fp = fopen(filepath, "r");
	fgets(sline, 512, fp);
	for(int i=0; i<rows-1; i++)
	{
		fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s ", 
			filename, &lat, &lon, &gx, &gy, &hei, &roll, &pitch, &yaw, &dtime, day, stime );
	
		roll  = roll  / PI * 180;
		pitch = pitch / PI * 180;
		yaw   = yaw   / PI * 180;

		//fscanf(fp, "%s  %lf %lf %lf  %lf %lf %lf  %lf ", 
		//	filename, &gx, &gy, &hei, &yaw, &pitch, &roll, &dtime);
		
		m_posArray[i].lat = lat;
		m_posArray[i].lon = lon;
		m_posArray[i].gx  = gx;
		m_posArray[i].gy  = gy;
		m_posArray[i].height = hei;
		m_posArray[i].roll   = roll;
		m_posArray[i].pitch  = pitch;
		m_posArray[i].yaw    = yaw;

		strcpy(m_posArray[i].filetitle, filename);
	}
	fclose(fp);

	return 0;
}
int CReadCarPos::GetPOS(int index, POSInfo& pos)
{
	pos = m_posArray[index];
	
	return 0;
}


////////////////////////////////////////////////////////////////
CReadCarPosOrbit::CReadCarPosOrbit()
{

}
CReadCarPosOrbit::~CReadCarPosOrbit()
{

}

int CReadCarPosOrbit::GetPOSNum()
{
	return m_posArray.size();
}

int CReadCarPosOrbit::ReadPOSData(char* filepath)
{
	char sline[512];

	int rows = GetFileRows(filepath);

	m_posArray.resize(rows - 1);

	char   filename[256];
	double lat = 0, lon = 0, hei = 0;
	double gx = 0, gy = 0;
	double roll, pitch, yaw;
	double dtime;
	char   day[256];
	char   stime[256];

	FILE* fp = fopen(filepath, "r");
	fgets(sline, 512, fp);
	for (int i = 0; i<rows - 1; i++)
	{
		//fscanf(fp, "%s %lf %lf %lf %lf %lf %lf %lf %lf %lf %s %s ", 
		//	filename, &lat, &lon, &gx, &gy, &hei, &roll, &pitch, &yaw, &dtime, day, stime );

		fscanf(fp, "%s  %lf %lf %lf  %lf %lf %lf  %lf ",
			filename, &gx, &gy, &hei, &yaw, &pitch, &roll, &dtime);

		m_posArray[i].lat = lat;
		m_posArray[i].lon = lon;
		m_posArray[i].gx = gx;
		m_posArray[i].gy = gy;
		m_posArray[i].height = hei;
		m_posArray[i].roll = roll;
		m_posArray[i].pitch = pitch;
		m_posArray[i].yaw = yaw;

		strcpy(m_posArray[i].filetitle, filename);
	}
	fclose(fp);

	return 0;
}
int CReadCarPosOrbit::GetPOS(int index, POSInfo& pos)
{
	pos = m_posArray[index];

	return 0;
}
