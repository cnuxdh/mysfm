#ifndef CV_BA_H
#define CV_BA_H

#include "export.hpp"
#include "defines.hpp"
#include "imagedata.hpp"

#include "sfm.h"

#define  PROJECTION_ESTIMATION_THRESHOLD 4.0
#define  MIN_PROJ_ERROR_THRESHOLD  8.0
#define  MAX_PROJ_ERROR_THRESHOLD  16.0
#define  MIN_MAX_MATCHES  16
#define  MIN_POINTS 20
#define  NUM_STDDEV 2.0 // 3.0 // 6.0
#define  OUTPUT_VERBOSE_STATS
#define  NUM_ERROR_BINS 10
#define  RAY_ANGLE_THRESHOLD 2.0
#define  CONSTRAIN_FOCAL_WEIGHT 0.0001




int  SetConstraints(camera_params_t* cam, bool bIsEstimateFocus, bool bIsEstimateDistort);

int  BundleAdjustAddAllNewPoints(int num_points, int num_cameras,
								int *added_order,
								camera_params_t *cameras,
								v3_t *points, v3_t *colors,
								double reference_baseline,
								vector<CImageDataBase*> imageData,
								vector<TrackInfo>& trackSeq,
								vector<ImageKeyVector> &pt_views,
								double max_reprojection_error = 16.0,
								int min_views = 2);

vector<ImagePair> FindCamerasWithNMatches(int match_threshold, 
										  int num_cameras, 
										  int *added_order,
										  vector<CImageDataBase*> imageData,
										  const vector<TrackInfo>& trackSeq,
										  const vector<ImageKeyVector> &pt_views);

int  FindCameraWithMostMatches(int num_cameras, int *added_order,
							  int &parent_idx, int &max_matches,
							  vector<CImageDataBase*> imageData,
							  const vector<TrackInfo>& trackSeq,
							  const vector<ImageKeyVector> &pt_views);

double RefinePoints(int num_points, v3_t *points, v2_t *projs,
					int *pt_idxs, camera_params_t *cameras,
					int *added_order,
					vector<CImageDataBase*> imageData,
					const vector<ImageKeyVector> &pt_views,
					camera_params_t *camera_out);

vector<int> RefineCameraParameters(int num_points, 
								   v3_t *points, v2_t *projs, 
								   int *pt_idxs, camera_params_t *camera,
								   double *error_out, 
								   bool adjust_focal,
								   bool remove_outliers,
								   bool optimize_for_fisheye,
								   bool estimate_distortion,
								   double min_proj_error_threshold,
								   double max_proj_error_threshold);

vector<int> RefineCameraAndPoints(int num_points,
								  v3_t *points, v2_t *projs,
								  int *pt_idxs, 
								  camera_params_t *cameras,
								  int *added_order,
								  vector<CImageDataBase*> imageData,
								  const vector<ImageKeyVector> &pt_views,
								  camera_params_t *camera_out,
								  bool remove_outliers);


bool FindAndVerifyCamera(int num_points, v3_t *points_solve, v2_t *projs_solve,
						 int *idxs_solve,
						 double *K, double *R, double *t, 
						 double proj_estimation_threshold,
						 double proj_estimation_threshold_weak,
						 std::vector<int> &inliers,
						 std::vector<int> &inliers_weak,
						 std::vector<int> &outliers);

camera_params_t BundleInitializeImage(const vector<TrackInfo>& trackSeq, 
									  int image_idx, int camera_idx,
									  int num_cameras, int num_points,
									  int *added_order, v3_t *points,
									  camera_params_t *parent,
									  camera_params_t *cameras, 
									  vector<CImageDataBase*> imageData,
									  vector<ImageKeyVector> &pt_views,
									  bool *success_out = NULL,
									  bool refine_cameras_and_points = false);


void InitBundleAdjust(vector<CImageDataBase*> imageData, vector<TrackInfo>& trackSeq);

/*interface function to invoke sba, revised by Donghai Xie, 2014.1.6
*/
double runSFMApi(int num_pts, int num_cameras, int start_camera,
			  bool fix_points, vector<CImageDataBase*> imageData,
			  camera_params_t *init_camera_params,
			  v3_t *init_pts, int *added_order, v3_t *colors,
			  std::vector<ImageKeyVector> &pt_views, double eps2 = 1.0e-12,
			  double *S = NULL, double *U = NULL, double *V = NULL,
			  double *W = NULL, bool remove_outliers = true);


//SFM based on track structure
int DoSFM( vector<TrackInfo> trackSeq, vector<ImgFeature> imageFeatures, 
		  vector<int> cameraIDOrder,
		  vector<CameraPara> &cameras);

//SFM based on the direct connection between 3d and 2d points
int DoSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
		   vector<CImageDataBase*> imageData, vector<int> cameraIDOrder,
		   vector<CameraPara>& cameras);

int TrackSeqImageOrder(vector<TrackInfo>& trackSeq, vector<int> cameraIDOrder);


/*
get the image matching result based on the connection between tracks and image points
input parameters:
	i,j: the image index of left and right image
*/
vector<MatchPairIndex> GetMatchList(int imgId1, int imgId2, vector<CImageDataBase*> imageData);


/* get the projections of one image among the current adjusting tracks  
*/
vector<int>  GetPointList(int imageId, vector<CImageDataBase*> imageData, vector<TrackInfo>& tracks);


//initialize the camera parameter based on current 3D structure
void BundleInitCamera();


class DLL_EXPORT CBABase
{
public:
	CBABase(){}
	virtual ~CBABase(){}
	

	//whole images optimization
	//virtual int RunAll(  ){return 0;} 
	
	//for one optimization 
	virtual int RunSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
		        vector<ImgFeature> imageFeatures,  vector<int> cameraIDOrder,
				vector<CameraPara>& cameras) {return 0;}

	//whole images optimization
	virtual int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<ImgFeature> imageFeatures, 
							 vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir){return 0;}
    //bundle adjustment: new interface
	virtual int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
		                     vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir){return 0;}

};

class DLL_EXPORT CSBA: public CBABase
{
public:
	CSBA();
	~CSBA();

	int RunSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
		vector<ImgFeature> imageFeatures,  vector<int> cameraIDOrder,
		vector<CameraPara>& cameras);

	//
	int BundleAdjust( int numCameras, vector<CameraPara>& cameras,vector<ImgFeature> imageFeatures, 
		              vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir);


	//new interface, including more input parameters
	int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
				      vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir);

private:
	//bool m_estimate_focal_length;
	//bool m_estimate_distortion;
	//bool m_constrain_focal;
	//bool m_use_point_constraints;
};




#endif