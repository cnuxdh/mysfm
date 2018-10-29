

#ifndef  ABSORI_HPP
#define  ABSORI_HPP

#include <vector>
using namespace std;

#include "commondata.h"

//coredll
#include "defs.h"

//cvlib
#include "defines.hpp"
#include "badata.hpp"


DLL_EXPORT double AbsOriOrthogonal(stAbsPOS& absPosParams, vector<stPOS>& camParas, vector<stTrack>& tracks);
DLL_EXPORT double AbsOriOrthogonal(stAbsPOS& absPosParams,vector<CameraPara>& camParas, vector<TrackInfo>& tracks);
DLL_EXPORT double AbsOriOrthogonal(stAbsPOS& absPosParams, vector<CameraPara>& camParas);


DLL_EXPORT int    GenerateProjectAxis(POINT3D za, POINT3D& xa, POINT3D& ya);
DLL_EXPORT int    RotationAlign(vector<POINT3D> srcPts, vector<POINT3D> dstPts, double* rotationMatrix);




#endif