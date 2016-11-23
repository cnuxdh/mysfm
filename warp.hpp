

#ifndef WARP_H
#define WARP_H

#include "export.hpp"




DLL_EXPORT int GenerateHomographyMatrix(double omiga, double kapa, double phi, 
	double tx, double ty, double tz,
	double nx, double ny, double nz, 
	double distance, double focal,
	double* H );


DLL_EXPORT int HomographyWarp(unsigned char* pSrc, int sht, int swd, 
	unsigned char* pDst, int& dht, int& dwd, double* H);


#endif