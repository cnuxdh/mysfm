

#ifndef  CV_REGISTER_H
#define  CV_REGISTER_H

#include "export.hpp"
#include "feature.hpp"
#include "defines.hpp"
#include "imagedata.hpp"


#define HOMOGRAPHY_THRESHOLD 4
#define HOMOGRAPHY_ROUNDS    256

enum MotionModel
{
	MotionRigid,
	MotionHomography,
};

//base class for feature matching
class DLL_EXPORT CMatchBase
{
public:
	CMatchBase(){}
	virtual ~CMatchBase(){}

	virtual int Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats ){return 0;}
	virtual int Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes ){return 0;}
	virtual int Match(char* file1, char* file2){return 0;};
	virtual int Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch ){return 0;}
	virtual int Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch, double* m ){return 0;}
};

class DLL_EXPORT CRationMatch: public CMatchBase
{
public:
	CRationMatch();
	~CRationMatch();

	//matching of image pair
	int Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats );
	int Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes);
	int Match(char* file1, char* file2);
};



class DLL_EXPORT CKNNMatch: public CMatchBase
{
public:
	CKNNMatch();
	~CKNNMatch();

	//matching of image pair
	int Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats );
	int Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes);
	int Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch );
};


class DLL_EXPORT CPanoMatch: public CMatchBase
{
public:
	CPanoMatch();
	~CPanoMatch();

	//matching of image pair
	int Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats );
	int Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes);
	int Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch );
};




class DLL_EXPORT CSiftMatch: public CMatchBase
{
public:
	CSiftMatch();
	~CSiftMatch();

	//matching of image pair
	int Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats );
	int Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes);
	int Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch );
	int Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch, double* m);
};



//////////////////////////////////////////////////////////////////////////
//base class for generating tracks
class DLL_EXPORT CGenerateTracksBase
{
public:
	CGenerateTracksBase(){}
	virtual ~CGenerateTracksBase(){}
	virtual  int GenerateTracks( vector<PairMatchRes> pairMatchs, vector<TrackInfo>& tracks ){return 0;}
	virtual  int GenerateTracks( vector<PairMatchRes> pairMatchs, vector<CImageDataBase*> imageData, vector<TrackInfo>& tracks ){return 0;}
};

class DLL_EXPORT CMyGenerateTrack: public CGenerateTracksBase
{
public:
	CMyGenerateTrack();
	~CMyGenerateTrack();
	int GenerateTracks( vector<PairMatchRes> pairMatchs, vector<TrackInfo>& tracks );
	int GenerateTracks( vector<PairMatchRes> pairMatchs, vector<CImageDataBase*> imageData, vector<TrackInfo>& tracks );
};




/////////////////  functions for retrieving information from tracks ///////////////////
//read tracks from the file
DLL_EXPORT int ReadTracks(char* filename, vector<TrackInfo>& tracks);

//output the tracks into the file
DLL_EXPORT int PrintTracks(vector<TrackInfo> tracks, char* filename);

/* input the image pair ( imgId1, imgId2), get the track including this pair
*/
DLL_EXPORT void GetMatch(int imgId1, int imgId2, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq);

/*search the projections in special image corresponding to 3D points
input:
    imgId:    image index for search
    tracks:   raw tracks including all images
    trackSeq: current tracks for bundle adjustment
*/
DLL_EXPORT int FindProjections(int imgId,  vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq);


//find the new image to be added 
DLL_EXPORT int FindNewImage(int numCam, vector<int> cameraIDOrder, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq);


//update the current tracks, i.e. add the new image projection into the current tracks
DLL_EXPORT int UpdateTracks(int imgId, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq);
DLL_EXPORT int UpdateTracks(int imgId, int newId, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq, 
							vector<ImgFeature>& imgFeatures, vector<Point3DDouble>& pt3, vector<Point2DDouble>& pt2);


//add new tracks
DLL_EXPORT int AddNewTracks(int imgId, vector<int> cameraIDOrder, vector<CameraPara> cameras,
				 vector<ImgFeature>& imgFeatures, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq);


//prune suspicious tracks
void PruneBadTracks(vector<TrackInfo>& trackSeq);

//set track id for projected feature points



/* Estimate a transform between two sets of keypoints */
vector<int> EstimateTransform(const vector<stPtFeature> &k1, 
							  const vector<stPtFeature> &k2, 
							  const vector<MatchPairIndex> &matches, 
							  MotionModel mm,
							  int nRANSAC, //rounds
							  double RANSACthresh, //threshold
							  double *Mout);





#endif
