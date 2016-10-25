

#ifndef CERES_BA
#define CERES_BA


#include "export.hpp"
#include "defines.hpp"
#include "imagedata.hpp"

#include "sfm.h"
//
#include"CalcAngle.h"

#include "ba.hpp"


#ifdef CERES_LIB
//ceres lib
#include "ceres/ceres.h"
#include "ceres/rotation.h"
#endif

#include <vector>
using namespace std;



/////////////////////////////////////////////////////////////////
//for ceres
// Templated pinhole camera model for used with Ceres.  The camera is
// parameterized using 9 parameters: 3 for rotation, 3 for translation, 1 for
// focal length and 2 for radial distortion. The principal point is not modeled
// (i.e. it is assumed be located at the image center).

#ifdef CERES_LIB
struct SFMReprojectionError 
{
	SFMReprojectionError(double observed_x, double observed_y)
		: observed_x(observed_x), observed_y(observed_y) {}


	template <typename T>
	bool operator()(const T* const cameraIn,  //3 interior camera parameters
					const T* const cameraOut, //6 outer camera parameters
					const T* const point,     //3 parameters for ground point
					T* residuals) const       //2 output residual parameters
	{
		T focal = cameraIn[0];
		T k1    = cameraIn[1];
		T k2    = cameraIn[2];

		T omiga  = cameraOut[0];
		T phi    = cameraOut[1];
		T kapa   = cameraOut[2];
		T t[3];
		t[0] = cameraOut[3];
		t[1] = cameraOut[4];
		t[2] = cameraOut[5];

		//int ht, wd;
		T R[9];
		GenerateRMatrixDirect(omiga, phi, kapa, R);

		T ix1,iy1;
		GrdToImgWithDistort(point[0], point[1], point[2], &ix1, &iy1, R, t, focal, T(0), T(0), k1, k2);
		
		residuals[0] = ix1 - T(observed_x);
		residuals[1] = iy1 - T(observed_y);
		
		//printf("rx: %lf , ry: %lf \n", residuals[0], residuals[1]);

		return true;
	}
	

	/*
	template <typename T>
	bool operator()(const T* const cameraParas, //9 parameters: 3 inner, 6 outer
		const T* const point,     //3 parameters for ground point
		T* residuals) const       //2 output residual parameters
	{
		T focal = cameraParas[0];
		T k1    = cameraParas[1];
		T k2    = cameraParas[2];

		T omiga  = cameraParas[3];
		T phi    = cameraParas[4];
		T kapa   = cameraParas[5];
		T t[3];
		t[0] = cameraParas[6];
		t[1] = cameraParas[7];
		t[2] = cameraParas[8];

		//printf("focal length: %lf \n", focal);

		//int ht, wd;
		T R[9];
		GenerateRMatrixDirect(omiga, phi, kapa, R);

		T ix1,iy1;
		GrdToImgWithDistort(point[0], point[1], point[2], &ix1, &iy1, R, t, focal, T(0), T(0), k1, k2);

		//double dx = ix1;
		//double dy = iy1;
		//cout<<ix1<<" "<<iy1<<"\n";
		//printf("ix: %ld  iy:%lf \n", ix1, iy1);

		residuals[0] = ix1 - T(observed_x);
		residuals[1] = iy1 - T(observed_y);


		//printf("rx: %lf , ry: %lf \n", residuals[0], residuals[1]);

		return true;
	}*/

	// Factory to hide the construction of the CostFunction object from
	// the client code.
	static ceres::CostFunction* Create(const double observed_x,
		const double observed_y) {
			return (new ceres::AutoDiffCostFunction<SFMReprojectionError, 2, 3, 6, 3>(
				new SFMReprojectionError(observed_x, observed_y)));
	}

	double observed_x;
	double observed_y;
};

int CeresBA( vector<TrackInfo> trackSeq, vector<ImgFeature> imageFeatures, 
			/*vector<int> cameraIDOrder,*/	vector<CameraPara> &cameras);


class DLL_EXPORT CCeresBA: public CBABase
{
public:
	CCeresBA();
	~CCeresBA();

	int RunSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
		vector<ImgFeature> imageFeatures,  vector<int> cameraIDOrder,
		vector<CameraPara>& cameras);

	//
	int BundleAdjust( int numCameras, vector<CameraPara>& cameras,vector<ImgFeature> imageFeatures, 
		vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir);


	//new interface, including more input parameters
	int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
		vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir);


	int BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<ImgFeature>& imageFeatures, 
		vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks);

private:
};

#endif


#endif
