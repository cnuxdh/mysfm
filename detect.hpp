
#ifndef CV_DETECT_HPP
#define CV_DETECT_HPP

#include"defines.hpp"

#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc.hpp"

#include <vector>
using namespace std;

using namespace cv;



DLL_EXPORT int  DetectJYZLoss(Mat Img, vector<Rect>& lossRects);




#endif

