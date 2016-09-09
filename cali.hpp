

#ifndef CV_CALI_H
#define CV_CALI_H

#include "sfm.h"

#include "defines.hpp"
#include "export.hpp"


void InitializeCameraParams(CameraPara camPara, camera_params_t &camera);
void InitializeCameraParams(camera_params_t &camera);

DLL_EXPORT void CreateChessBoard();



#endif