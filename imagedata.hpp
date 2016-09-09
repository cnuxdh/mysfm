

#ifndef IMAGEDATA_H
#define IMAGEDATA_H

#include "export.hpp"
#include "defines.hpp"



#define SIFT_FEATURE       0
#define SIFT_FLOAT_FEATURE 1
#define SURF_FEATURE       2

class DLL_EXPORT CImageDataBase
{
public:
	CImageDataBase(){}
	virtual ~CImageDataBase(){}

	virtual int Load(char* filename){return 0;}
	virtual int DetectPtFeature(int featureType){return 0;}  


	virtual int GetHt(){return 0;}
	virtual int GetWd(){return 0;}

	virtual int GetFeatNumber(){return 0;}   //the number of feature points in the image
	virtual FeatPoint GetKeyPoint(int index) //return the point information 
	{
		FeatPoint pt;
		return pt;
	}
	virtual void SetFeatExtra(int index, int extra){}
	//virtual int AddTrackIdx(int idx){return 0;}
	virtual vector<int> GetTrackSeq()
	{
		vector<int> vi;
		return vi;
	}

	virtual int SetTrackIdx(int ptIndex, int idx){return 0;}
	virtual int GetPointTrackIdx(int ptIndex){return 0;}
	virtual ImgFeature& GetImageFeature()
	{
		ImgFeature features;
		return features;
	}

	virtual vector<int> GetPtSeq()
	{
		vector<int> vi;
		return vi;
	}

	//for initial focus
	virtual int    IsHasInitFocus(){return 0;}
	virtual void   SetInitFocus(double focus){}
	virtual double GetInitFocus(){return 0;}
};


class DLL_EXPORT CImageFeature: public CImageDataBase
{
public:
	CImageFeature();
	~CImageFeature();

	int Load(char* filename);
	int DetectPtFeature(int featureType);    

	int GetHt();
	int GetWd();

	//for feature
	int GetFeatNumber();   
	FeatPoint GetKeyPoint(int index);
	//int  GetExtra();
	void SetFeatExtra(int index, int extra);
	ImgFeature& GetImageFeature();

	//for track
	//int AddTrackIdx(int idx);  //save the track index into the vector
	vector<int> GetTrackSeq(); //get all track index
	int SetTrackIdx(int ptIndex, int idx); //set the track of single point
	int GetPointTrackIdx(int ptIndex);

	vector<int> GetPtSeq();
	
	int    IsHasInitFocus();
	void   SetInitFocus(double focus);
	double GetInitFocus();
	
private:
	ImgFeature   m_imageFeature;
	char         m_strFileName[256];
	vector<int>  m_visible_points;  //the track id to each feature point
    vector<int>  m_visible_keys;    //the point id to each feature point

	//resized height and width
	int          m_ht;
	int          m_wd;
	

	double       m_initFocus;
	int          m_isInitFocus;
};


#endif