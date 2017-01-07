
#include "badata.hpp"


ImgFeature::ImgFeature()
{
	featPts.clear();
}

ImgFeature::~ImgFeature()
{

}

Point2DDouble ImgFeature::GetCenteredPt(int ptIndex)
{
	Point2DDouble cp;

	cp.p[0] = featPts[ptIndex].cx;
	cp.p[1] = featPts[ptIndex].cy;

	return cp;
}
Point2DDouble ImgFeature::GetTopLeftPt(int ptIndex)
{
	Point2DDouble cp;

	cp.p[0] = featPts[ptIndex].x;
	cp.p[1] = featPts[ptIndex].y;

	return cp;
}



TrackInfo::TrackInfo()
{
	views.clear();
}

TrackInfo::~TrackInfo()
{

}