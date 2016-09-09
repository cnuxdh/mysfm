
#include "assert.h"

#include "ImageData.hpp"

#include "feature.hpp"
#include "sift.hpp"
//#include "surf.hpp"


#ifdef OPENCV_1X 
#include "cv.h"
#include "highgui.h"
#include "cxcore.h"
#else
//opencv
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/calib3d/calib3d.hpp"
#endif


CImageFeature::CImageFeature()
{
	m_imageFeature.featPts.clear();
	m_visible_points.clear();
	m_visible_keys.clear();

	m_initFocus = 0;
	m_isInitFocus = 0;
}
CImageFeature::~CImageFeature()
{

}

int CImageFeature::Load(char* filename)
{
	strcpy(m_strFileName, filename);	
	
	IplImage* pImage = cvLoadImage(filename);
	int ht = pImage->height;
	int wd = pImage->width;
    cvReleaseImage(&pImage);

	//get the dimension after resize
	int dstHt,dstWd;
	GetResizeDimension(ht, wd, dstHt, dstWd);

	m_ht = dstHt;
	m_wd = dstWd;

	return 1;
}


int CImageFeature::DetectPtFeature(int featureType)
{
	CFeatureBase* pFeatDetect = NULL;

	switch (featureType)
	{
	case 0:
		pFeatDetect = new CSIFT();
		break;
	case 1:
		pFeatDetect = new CSIFTFloat();
		break;
	case 2:
		pFeatDetect = NULL; //new CSURF();
		break;
	}

	pFeatDetect->Detect(m_strFileName, m_ht, m_wd, m_imageFeature);

	delete pFeatDetect;

	return 1;
}

int CImageFeature::GetHt()
{
	return m_ht;
}

int CImageFeature::GetWd()
{
	return m_wd;
}

int CImageFeature::GetFeatNumber()
{
	return m_imageFeature.featPts.size();
}

FeatPoint CImageFeature::GetKeyPoint(int index)
{
	FeatPoint pt;

	pt.cx = m_imageFeature.featPts[index].cx;
	pt.cy = m_imageFeature.featPts[index].cy;
	pt.x  = m_imageFeature.featPts[index].x;
	pt.y  = m_imageFeature.featPts[index].y;
	pt.track = m_imageFeature.featPts[index].trackIdx;
	pt.extra = m_imageFeature.featPts[index].extra;
	
	return pt;
}


/*
int CImageFeature::AddTrackIdx(int idx)
{
	//m_trackVec.push_back(idx);
	return 1;
}
*/

vector<int> CImageFeature::GetTrackSeq()
{
	return m_visible_points;
}
vector<int> CImageFeature::GetPtSeq()
{
	return m_visible_keys;
}

//establish the connection between the track and image point
int CImageFeature::SetTrackIdx(int ptIndex, int idx)
{
	assert( ptIndex<m_imageFeature.featPts.size() );
	//assert( m_imageFeature.featPts[ptIndex].trackIdx < 0 );

	//if this image point has been connected to one track, ignore next steps, added by Donghai, Xie, 2014.1.10
	if( m_imageFeature.featPts[ptIndex].trackIdx >= 0 )
		return 0;

	m_imageFeature.featPts[ptIndex].trackIdx = idx;	
	m_visible_points.push_back(idx);
	m_visible_keys.push_back(ptIndex);
	
	return 1;
}
void CImageFeature::SetFeatExtra(int index, int extra)
{
	assert( index<m_imageFeature.featPts.size() );
	m_imageFeature.featPts[index].extra = extra;	
}

int CImageFeature::GetPointTrackIdx(int ptIndex)
{
	assert( ptIndex<m_imageFeature.featPts.size() );

	return m_imageFeature.featPts[ptIndex].trackIdx;
}

ImgFeature& CImageFeature::GetImageFeature()
{
	return m_imageFeature;
}

int    CImageFeature::IsHasInitFocus()
{
	return m_isInitFocus;
}
void   CImageFeature::SetInitFocus(double focus)
{
	m_initFocus = focus;
	m_isInitFocus = 1;
}

double CImageFeature::GetInitFocus()
{
	return m_initFocus;
}