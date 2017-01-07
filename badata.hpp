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
	
public:
	int id;
	int ht,wd;
	vector<stPtFeature> featPts;
};


class DLL_EXPORT TrackInfo
{
public:
	TrackInfo();
	~TrackInfo();
	
	int GetRemapTrackId(){ return extra; }
	int GetValidValue(){ return valid;}

	void SetGround(Point3DDouble pt){ grd = pt; }
	Point3DDouble GetGround(){return grd;}

	int GetImageKeySum(){return views.size();}
	
	//get and add the image-key of index
	ImageKey GetImageKey(int index){ return views[index]; } 
	void AddImageKey(ImageKey ik){ views.push_back(ik); }
	
	double GetProjectionError(){return derror;}

	void   Clear(){views.clear();}

public:
	int id;
	int extra;              //save extra information (such as the index of image feature point)
	int valid;              //1:valid , 0:invalid 
	double derror;    
	Point3DDouble  grd;
	ImageKeyVector views;
};


#endif



