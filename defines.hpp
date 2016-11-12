#ifndef CV_DEFINES_H
#define CV_DEFINES_H

#include "commondata.h"


// interface for feature detection
#include <vector>
using namespace std;

#define PI 3.1415926



typedef struct stMatrixByte
{
	int ht,wd;
	unsigned char* pbuffer;
}MatrixByte;


typedef struct stiPoint
{
	int x,y;
}iPoint;

//point structure
typedef struct stFeatPoint
{
	double x, y;   //left-top
	double cx,cy;  //normalized to the center
	int extra;     /* 4 bytes of extra storage */
	int track;  /* Track index this point corresponds to */
}FeatPoint;


/* 2D vector of doubles */
typedef struct stPoint2DDouble
{
	double p[2];
	double x,y;
} Point2DDouble;

/* 3D vector of doubles */
typedef struct stPoint3DDouble
{
	double p[3];
	int    extra;
} Point3DDouble;


// feature structure for point description and matching
typedef struct stPtFeature
{
	int    id;
	float  ori;
	int    key_octave,key_intvl;  // 特征点所在的组和层 
	float  sub_intvl;             // 精确定位后特征点所在层的改正数
	float  scl;                   // 关键点的尺度
	float  scl_octave;            //
	float  row,col;               // image coordinate, for backup
	float  x,y;                   // coordinates of point  ( origin is topleft )
	float  cx,cy;                 // normalized coordinate ( origin is the center of image, y axis is up)
	vector<float> feat;           // feature vector
	int    trackIdx;              // corresponding track index
	int    extra;                 // 
}PtFeature;



/*
//feature pts in one image
typedef struct stImgFeature
{
	int id;
	int ht,wd;
	vector<stPtFeature> featPts;
}ImgFeature;
*/

//match pair
typedef struct stMatchPair
{
	int l,r; //the image point index of left image and right image
}MatchPairIndex;


typedef std::pair<int,int> ImagePair;



/*
//simple point class for SBA
class Keypoint 
{
public:    
Keypoint(){ m_x = 0.0; m_y = 0.0; m_extra = -1; m_track = -1; }
Keypoint(float x, float y):m_x(x), m_y(y)
{ m_r = 0; m_g = 0; m_b = 0; }
virtual ~Keypoint() {}   
virtual unsigned char *GetDesc(){return NULL;}
float m_x, m_y;                   // Subpixel location of keypoint.
unsigned char m_r, m_g, m_b;      // Color of this key 
int m_extra;                      // 4 bytes of extra storage 
int m_track;                      // Track index this point corresponds to 
};
*/


//structure for saving the matches of two images
typedef struct stPairMatchRes
{
	int    lId,rId;                 // image id
	double inlierRatio;             // the inlier ratio 
	vector<MatchPairIndex> matchs;  // matching point indexs 
}PairMatchRes;


//define one track
typedef std::pair<int,int>    ImageKey;         //the projection  of a track on one image, <imageId, featureID>
typedef std::vector<ImageKey> ImageKeyVector;   //all projections of a track 
//typedef struct stTrackInfo
//{
//	int id;
//	int extra;              //save extra information (such as the index of image feature point)
//	int valid;              //1:valid , 0:invalid 
//	double derror;    
//	Point3DDouble  grd;
//	ImageKeyVector views;
//}TrackInfo;


//interior parameters




//structure for Camera
typedef struct stCameraPara
{
	double focus;
	double k1,k2;
	double u0,v0;
	double K[9];
	double R[9];
	double t[3];
	double xs,ys,zs; //absolute position of camera
	double lon,lat;  //the lon/lat of position of camera
	double ax,ay,az; //rotation angle around x,y,z (Eular degree)
	int    rows,cols;    //the size of image in pixels
}CameraPara;


#define MIN_INLIERS_EST_PROJECTION 6 /* 7 */ /* 30 */ /* This constant needs
// adjustment
//#define INIT_REPROJECTION_ERROR 16.0 /* 6.0 */ /* 8.0 */


//////////////////////////////////////////////////////////////////////////
//fro Bundle Adjustment

typedef struct STRU_TRACK
{
	double x,y,z;        //position
	double nx,ny,nz;     //normal 
	unsigned char r,g,b; //color
	vector<POINT2> imgpt; 
	vector<int> imgid;
	vector<int> ptid;
}stTrack;


#define HOMOGRAPHY_THRESHOLD 8
#define HOMOGRAPHY_ROUNDS    256

typedef enum 
{
	MotionNone,
	MotionRigid,
	MotionHomography,
	MotionFundametalMatrix,
}MotionModel;



#endif







