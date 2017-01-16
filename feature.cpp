
#include "math.h"

#include "feature.hpp"


//get the image dimension after scale
void GetResizeDimension(int srcHt, int srcWd, int& dstHt, int& dstWd)
{

	if( srcHt<DETECT_IMAGE_HT )
	{
		dstHt = srcHt;
		dstWd = srcWd;
	}
	else
	{
		double ratio = (double)(DETECT_IMAGE_HT) / (double)(srcHt);
		dstHt = ceil(srcHt*ratio);
		dstWd = ceil(srcWd*ratio);
	}
}

