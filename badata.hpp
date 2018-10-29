/* Header file for data structure definition for bundle adjustment
   In the bundle adjustment, there are several important structures:
   feature point, track point,...

   feature point: saving and managing the data structure for feature points
		detected from the images
   track point:   saving and managing the data structure for track points 
   
*/


#ifndef BADATA_HPP
#define BADATA_HPP

#include "defines.hpp"
#include "export.hpp"


//opencv
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <opencv2/imgproc/imgproc.hpp>
using namespace cv;


 class DLL_EXPORT ImgFeature
{
public:
	ImgFeature();
	~ImgFeature();

	//get the sum of feature points in this image
	int GetFeatPtSum(){return featPts.size();}

	//get the track index connecting with the feature point
	int GetPtTrackIndex(int ptIndex){return featPts[ptIndex].extra;}
	void SetTrackIndex(int ptIndex, int nTrackIndex){ featPts[ptIndex].extra = nTrackIndex;}

	//get the centered x,y of feaure point
	Point2DDouble GetCenteredPt(int ptIndex);
	Point2DDouble GetTopLeftPt(int ptIndex);
	
	//set the feature points data
	void SetFeaturePoints(vector<stPtFeature>& srcFeatPts){ featPts=srcFeatPts; }

	//
	void SetImageDimension(int height, int width){ ht=height; wd=width; }
	
	int  AddFeatPt( stPtFeature fp ){ featPts.push_back(fp); return 0; }

	void  Clear()
	{    
		featPts.clear(); 
		mDescriptors.empty();
	}

	stPtFeature GetFeatPt(int i){return featPts[i];}

	//save the feature into file, added by xdh, 2018.5.7
	int SaveImgFeaturePts(string filepath);
	int ReadImgFeaturePts(string filepath);


public:
	int id;
	int ht,wd;
	vector<stPtFeature> featPts;  //
	//vector<KeyPoint> mKeys;       //for orb 
	Mat mDescriptors;               //for orb

	
};


//a track includes all projections to the images, and the 3d point,...
class DLL_EXPORT TrackInfo
{
public:
	TrackInfo();
	~TrackInfo();
	
	int GetRemapTrackId(){ return extra; } //the id to the track in other track sequence
	int GetValidValue(){ return valid;}

	void SetGround(Point3DDouble pt){ grd = pt; }
	Point3DDouble GetGround(){return grd;}

	int GetImageKeySum(){return views.size();}
	
	//get and add the image-key of index
	ImageKey GetImageKey(int index){ return views[index]; } 
	void AddImageKey(ImageKey ik){ views.push_back(ik); }
	
	double GetProjectionError(){return derror;}

	void   Clear(){views.clear();}

	int    GetViews(){ return views.size(); }

	void   SetAsCtrlPt(){ ctrl = 1; }
	int    IsCtrl(){ return ctrl; }

public:
	int id;
	int extra;              //save extra information (such as the index of  track)
	int valid;              //1:valid , 0:invalid 
	int ctrl;               //1:control point, 0: free points
	double derror;    
	Point3DDouble  grd;
	ImageKeyVector views;
};


int GetGoodTracks(vector<TrackInfo> srcTracks, vector<TrackInfo>& goodTracks);


#endif



