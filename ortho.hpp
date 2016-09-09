

#ifndef ORTHO_IMAGE_H
#define ORTHO_IMAGE_H

#include "export.hpp"
#include "defines.hpp"


//corelib
#include "commondata.h"


class DLL_EXPORT COrthoBase
{
public:
	COrthoBase(){}
	virtual ~COrthoBase(){}

	virtual int GenerateOrthoImage(CameraPara camInfo, char* inputFile, char* outFile){ return 0;}

};

class DLL_EXPORT COrthoUAV: public COrthoBase
{
public:
	COrthoUAV();
	~COrthoUAV();

	int GenerateOrthoImage(CameraPara camInfo, char* inputFile, char* outFile);

private:
		
};


void  GenerateOrtho(double* R, double* T, 
					double f, double x0, double y0,
					double l, double r, double t, double b,
					double resolution,
					COLOR* image, int ht, int wd,						 
					unsigned char* ir, 
					unsigned char* ig, 
					unsigned char* ib, 
					int oht, int owd);

#endif
