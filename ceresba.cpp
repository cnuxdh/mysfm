
#include "ceresba.hpp"

#include "float.h"
#include "math.h"

#include "ba.hpp"
#include "cali.hpp"
#include "relativepose.hpp"
#include "distortion.hpp"
#include "bundlerio.hpp"


//matrix
#include "matrix/matrix.h"



///////////////////////////////////////////////////////////////////////////////////////
//select the new image to insert into the bundle adjustment
int SelectNewImage(vector<int> cameraVisited,
	vector<TrackInfo>& tracks, 
	vector<TrackInfo>& trackSeqNew )
{
	int newCamera = -1;

	int numCamera = cameraVisited.size();

	vector<int> projFreq;
	projFreq.resize(numCamera, 0);

	for(int i=0; i<trackSeqNew.size(); i++)
	{
		int nTrackId = trackSeqNew[i].extra;

		for(int j=0; j<tracks[nTrackId].views.size(); j++)
		{
			int nImageId = tracks[nTrackId].views[j].first;
			int nPtId = tracks[nTrackId].views[j].second;

			//ignore the cameras already added
			if( cameraVisited[nImageId] > 0 )
				continue;

			projFreq[nImageId] ++;
		}
	}

	int nMaxFreq = 1;
	for(int i=0; i<numCamera; i++)
	{
		if(nMaxFreq < projFreq[i])
		{
			nMaxFreq = projFreq[i];
			newCamera = i;
		}
	}
	return newCamera;
}

//update the current track sequence
int UpdateBATracks( int newCameraIndex, vector<int> cameraVisited,
	vector<ImgFeature>&  imageFeatures, 
	vector<TrackInfo>& tracks, 
	vector<TrackInfo>& trackSeqNew)
{
	int nFeatPtNum = imageFeatures[newCameraIndex].featPts.size();
	for(int i=0; i<nFeatPtNum; i++)
	{
		int nTrackIndex = imageFeatures[newCameraIndex].featPts[i].extra;

		if(nTrackIndex<0)
			continue;

		int nNewTrackIndex = tracks[nTrackIndex].extra;

		ImageKey ik;
		ik.first  = newCameraIndex;
		ik.second = i; 

		if(nNewTrackIndex>=0) //corresponding to the new Track sequence
		{
			trackSeqNew[nNewTrackIndex].views.push_back(ik);
		}
		else //adding a new track into the new Track sequence
		{
			TrackInfo newTrack;
			newTrack.extra = nTrackIndex;


			newTrack.views.push_back(ik);

			for(int k=0; k<tracks[nTrackIndex].views.size(); k++)
			{
				int imageId = tracks[nTrackIndex].views[k].first;
				int nPtId = tracks[nTrackIndex].views[k].second;
				if( cameraVisited[imageId]>0)
				{
					ik.first = imageId;
					ik.second = nPtId;
					newTrack.views.push_back(ik);
				}
			}

			if(newTrack.views.size()>1)
			{
				tracks[nTrackIndex].extra = trackSeqNew.size(); //added by xdh, 2016.10.31
				trackSeqNew.push_back(newTrack);
			}
		}
	}		

	return 0;
}

//calculate the ground point of track seq
int CaculateTrackSeqGrd(vector<ImgFeature>&  imageFeatures, 
	vector<TrackInfo>&   tracks, 
	vector<CameraPara>&  cameras,
	bool explicit_camera_centers)
{
	CTriangulateBase* pTriangulate = new CTriangulateCV();
	
	//vector<TrackInfo> newTracks;

	for(int i=0; i<tracks.size(); i++)
	{
		bool estimate_distortion = true;

		int num_views = (int) tracks[i].views.size();

		vector<Point2DDouble> pts;
		vector<CameraPara> cams;
		Point3DDouble gps; 
		double ferror;

		//collect all projection points of one track
		for (int j = 0; j < num_views; j++) 
		{
			int camera_idx = tracks[i].views[j].first;
			int key_idx    = tracks[i].views[j].second;

			double cx = imageFeatures[camera_idx].featPts[key_idx].cx;
			double cy = imageFeatures[camera_idx].featPts[key_idx].cy;

			Point2DDouble pt;
			pt.x = cx;
			pt.y = cy;
			pt.p[0] = cx;
			pt.p[1] = cy;
			pts.push_back(pt);

			cams.push_back( cameras[camera_idx] );
		}

		pTriangulate->Triangulate(pts, cams, gps, true, ferror);

		tracks[i].derror = ferror;
		tracks[i].grd = gps;
	}

	return 0;
}

int CalculateNewCamParas(int nCameraId, 
	vector<ImgFeature>&  imageFeatures,
	vector<TrackInfo>& trackSeqNew, 
	vector<TrackInfo>& tracks,
	CameraPara& cam )
{
	vector<Point3DDouble> pt3;
	vector<Point2DDouble> pt2;

	//collect the track corresponding to the new camera
	for(int i=0; i<trackSeqNew.size(); i++)
	{
		int nTrackId = trackSeqNew[i].extra;              //find the original track

		//ignore the track with large projection error
		if( trackSeqNew[i].derror > 16 )
			continue;

		for(int j=0; j<tracks[nTrackId].views.size(); j++)
		{
			int nImageId = tracks[nTrackId].views[j].first;
			int nPtId = tracks[nTrackId].views[j].second;

			if(nImageId == nCameraId)
			{
				pt3.push_back( trackSeqNew[i].grd );

				Point2DDouble p2;
				p2.p[0] = imageFeatures[nImageId].featPts[nPtId].cx;
				p2.p[1] = imageFeatures[nImageId].featPts[nPtId].cy;
				pt2.push_back(p2);
			}
		}
	}

	//ransac DLT calculation
	CPoseEstimationBase* pPE = new CDLTPose();
	int r = pPE->EstimatePose(pt3, pt2, cam);
	delete pPE;

	//if successful, do bundle adjustment
	if(r==0)
	{
		printf("\n optimization for DLT results: \n");

		//remove the wrong projections
		vector<Point3DDouble> inlierPt3;
		vector<Point2DDouble> inlierPt2;
		for(int i=0; i<pt3.size(); i++)
		{
			Point2DDouble projPt;
			GrdToImg(pt3[i], projPt, cam);

			double dx = projPt.p[0] - pt2[i].p[0];
			double dy = projPt.p[1] - pt2[i].p[1];
			double projError = sqrt(dx*dx + dy*dy);

			if(projError<16)
			{
				inlierPt3.push_back(pt3[i]);
				inlierPt2.push_back(pt2[i]);
			}
		}		

		CeresBA(inlierPt3, inlierPt2, cam);
		printf("\n");
	}
	else
	{
		return -1;
	}

	return 0;
}



#ifdef CERES_LIB

int CeresBA( vector<TrackInfo> trackSeq, vector<ImgFeature> imageFeatures, 
			 /*vector<int> cameraIDOrder,*/	vector<CameraPara> &cameras)
{
	//int  num_pts = trackSeq.size();
	int  num_cameras = cameras.size(); 

	//collect camera parameters
	double* pInteriorParams = new double[3];          //focal length, k1, k2
	pInteriorParams[0] = cameras[0].focus;
	pInteriorParams[1] = 0;
	pInteriorParams[2] = 0;	
	
	double* pOuterParams = new double[num_cameras*6]; //omiga, phi, kapa, t0,t1,t2
	for(int i=0; i<num_cameras; i++)
	{
		pOuterParams[i*6]   = cameras[i].ax;
		pOuterParams[i*6+1] = cameras[i].ay;
		pOuterParams[i*6+2] = cameras[i].az;
		pOuterParams[i*6+3] = cameras[i].t[0];
		pOuterParams[i*6+4] = cameras[i].t[1];
		pOuterParams[i*6+5] = cameras[i].t[2];
	}
	

	/*
	double* pOuterParams = new double[num_cameras*9]; //omiga, phi, kapa, t0,t1,t2
	for(int i=0; i<num_cameras; i++)
	{
		pOuterParams[i*9]   = cameras[i].focus;
		pOuterParams[i*9+1] = 0;
		pOuterParams[i*9+2] = 0;
		pOuterParams[i*9+3] = cameras[i].ax;
		pOuterParams[i*9+4] = cameras[i].ay;
		pOuterParams[i*9+5] = cameras[i].az;
		pOuterParams[i*9+6] = cameras[i].t[0];
		pOuterParams[i*9+7] = cameras[i].t[1];
		pOuterParams[i*9+8] = cameras[i].t[2];
	}*/

	//find the reasonable tracks 
	vector<int> goodTrackIndex;
	for(int i=0; i<trackSeq.size(); i++)
	{
		if(trackSeq[i].derror<4)
			goodTrackIndex.push_back(i);
	}
	int num_pts = goodTrackIndex.size();

	//collect ground points
	double* grdPt = new double[num_pts*3]; //(double*)malloc( _pts*3*sizeof(double) );
	for(int i=0; i<num_pts; i++)
	{
		int index = goodTrackIndex[i];
		grdPt[i*3]   = trackSeq[index].grd.p[0];
		grdPt[i*3+1] = trackSeq[index].grd.p[1];
		grdPt[i*3+2] = trackSeq[index].grd.p[2];
	}
	
	//generate observation points
	int nProjection = 0;
	for(int i=0; i<goodTrackIndex.size(); i++)
	{
		int index = goodTrackIndex[i];
		int nview = trackSeq[index].views.size();
		nProjection += nview;
	}

	vector<int> vecCamIndex;    // camera index for each projection point
	vector<int> vecTrackIndex;  // track point index for each projection point
	vecCamIndex.resize(nProjection);
	vecTrackIndex.resize(nProjection);

	double* projections = new double[nProjection*2]; //(double*)malloc(nProjection*2*sizeof(double));
	int ip = 0;
	//for(int i=0; i<trackSeq.size(); i++)
	for(int i=0; i<goodTrackIndex.size(); i++)
	{   
		int index = goodTrackIndex[i];
		int nview = trackSeq[index].views.size();

		for(int j=0; j<nview; j++)
		{
			int cameraID = trackSeq[index].views[j].first;
			int ptIndex  = trackSeq[index].views[j].second;

			projections[ip*2]   = imageFeatures[cameraID].featPts[ptIndex].cx; 
			projections[ip*2+1] = imageFeatures[cameraID].featPts[ptIndex].cy; 			
			
			vecTrackIndex[ip]  = i; 
			vecCamIndex[ip]    = cameraID;

			ip++;
		}
	}

	printf("3D Points: %d   projection number: %d \n", num_pts, nProjection);

	//google::InitGoogleLogging(NULL);

	//nProjection = nProjection - 1;

	ceres::Problem problem;
	//invoke ceres functions, each time adding one projection
	for(int i=0; i<nProjection; i++)
	{
		ceres::CostFunction* cost_function =
			SFMReprojectionError::Create(projections[2 * i + 0], projections[2 * i + 1]);

		int cameraId = vecCamIndex[i];
		int trackId  = vecTrackIndex[i];
		
		//printf("camera id: %d  track id: %d \n", cameraId, trackId);

		problem.AddResidualBlock(cost_function,
			NULL /* squared loss */,
			pInteriorParams,			    //inner camera parameters: 3
			pOuterParams + cameraId*6,      //outer camera parameters: 6
			grdPt + trackId*3);			    //parameters for point :   3
	}

	//save the optimized results
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	std::cout<<"optimized results: "<<"\n";

	printf("Inner parameters... \n");
	for(int i=0; i<3; i++)
		std::cout<<pInteriorParams[i]<<"  ";
	std::cout<<endl;
	
	printf("Ounter parameters ...\n");
	for(int i=0; i<num_cameras; i++)
	{
		for(int j=0; j<6; j++)
		{
			std::cout<<pOuterParams[i*6+j]<<"   ";
		}
		
		cameras[i].focus = pInteriorParams[0];
		cameras[i].k1    = pInteriorParams[1];
		cameras[i].k2    = pInteriorParams[2];
		
		cameras[i].ax   = pOuterParams[i*6] ;
		cameras[i].ay   = pOuterParams[i*6+1];
		cameras[i].az   = pOuterParams[i*6+2];
		cameras[i].t[0] = pOuterParams[i*6+3];
		cameras[i].t[1] = pOuterParams[i*6+4];
		cameras[i].t[2] = pOuterParams[i*6+5];

		//update the rotation matrix
		double R[9];
		GenerateRMatrixDirect(cameras[i].ax, cameras[i].ay, cameras[i].az, R);
		memcpy(cameras[i].R, R, sizeof(double)*9);

		std::cout<<endl;
	}

	for(int i=0; i<num_pts; i++)
	{   
		int index = goodTrackIndex[i];

		trackSeq[index].grd.p[0] = grdPt[i*3];
		trackSeq[index].grd.p[1] = grdPt[i*3+1];
		trackSeq[index].grd.p[2] = grdPt[i*3+2];
	}

	return 0;
}


int CeresBA( vector<Point3DDouble>& grdPts, vector<Point2DDouble>& imgPts, CameraPara& camera )
{

	double  pInteriorParams[3];
	double  pOuterParams[6];

	int nProjection = imgPts.size();

	pInteriorParams[0] = camera.focus;
	pInteriorParams[1] = camera.k1;
	pInteriorParams[2] = camera.k2;

	pOuterParams[0] = camera.ax;
	pOuterParams[1] = camera.ay;
	pOuterParams[2] = camera.az;
	pOuterParams[3] = camera.t[0];
	pOuterParams[4] = camera.t[1];
	pOuterParams[5] = camera.t[2];
		
	ceres::Problem problem;
	//invoke ceres functions, each time adding one projection
	for(int i=0; i<nProjection; i++)
	{
		ceres::CostFunction* cost_function =
			RefineCameraError::Create( imgPts[i].p[0] , imgPts[i].p[1], 
			grdPts[i].p[0], grdPts[i].p[1], grdPts[i].p[2]);

		problem.AddResidualBlock(cost_function,
			NULL /* squared loss */,
			pInteriorParams,   //inner camera parameters: 3
			pOuterParams       //outer camera parameters: 6
			);			   
	}

	//save the optimized results
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	std::cout << summary.FullReport() << "\n";

	std::cout<<"optimized results: "<<"\n";
	printf("Inner parameters... \n");
	for(int i=0; i<3; i++)
		std::cout<<pInteriorParams[i]<<"  ";
	std::cout<<endl;

	printf("Ounter parameters ...\n");

	for(int i=0; i<6; i++)
	{
		std::cout<<pOuterParams[i]<<"   ";
	}
	
	//save the optimized results
	camera.focus = pInteriorParams[0];
	camera.k1    = pInteriorParams[1];
	camera.k2    = pInteriorParams[2];
	camera.ax   = pOuterParams[0] ;
	camera.ay   = pOuterParams[1];
	camera.az   = pOuterParams[2];
	camera.t[0] = pOuterParams[3];
	camera.t[1] = pOuterParams[4];
	camera.t[2] = pOuterParams[5];

	//update the rotation matrix
	double R[9];
	GenerateRMatrixDirect(camera.ax, camera.ay, camera.az, R);
	memcpy(camera.R, R, sizeof(double)*9);

	std::cout<<endl;

	return 0;
}


#endif


CCeresBA::CCeresBA()
{
}
CCeresBA::~CCeresBA()
{
}
int CCeresBA::RunSFM( vector<Point3DDouble> pt3, vector<ImageKeyVector> ptViews, 
	vector<ImgFeature> imageFeatures,  vector<int> cameraIDOrder,
	vector<CameraPara>& cameras)
{
	return 0;
}
//
int CCeresBA::BundleAdjust( int numCameras, vector<CameraPara>& cameras,vector<ImgFeature> imageFeatures, 
	vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir)
{
	return 0;
}
//new interface, including more input parameters
int CCeresBA::BundleAdjust(int numCameras, vector<CameraPara>& cameras, vector<CImageDataBase*> imageData, 
	vector<PairMatchRes> pairMatchs, vector<TrackInfo> tracks, char* outDir)
{
	return 0;
}
int CCeresBA::BundleAdjust(int numCameras, 
	vector<CameraPara>&   cameras,     
	vector<ImgFeature>&   imageFeatures, //feature points
	vector<PairMatchRes>& pairMatchs,    //matching results
	vector<TrackInfo>& tracks)           //track points ( connecting points )
{
	//1. select the initial pair
	double maxInlier = 0;
	int    index = 0;
	for(int i=0; i<pairMatchs.size(); i++)
	{
		if( maxInlier<pairMatchs[i].matchs.size() )
		{
			maxInlier = pairMatchs[i].matchs.size();
			index = i;
		}
	}

	vector<int> cameraVisited;
	cameraVisited.resize(numCameras, 0);

	int leftImageId  = pairMatchs[index].lId;
	int rightImageId = pairMatchs[index].rId;
	vector<TrackInfo> trackSeq;  //current tracks for bundle adjustment
	GetMatch(leftImageId, rightImageId, tracks, trackSeq);    


	cameraVisited[leftImageId]  = 1;
	cameraVisited[rightImageId] = 1;

	//2. relative pose estimation
	//2.1 pose estimation
	vector<Point2DDouble> lpts,rpts;	
	int nMatch = trackSeq.size();
	for(int i=0; i<nMatch; i++)
	{
		int imgId1 = trackSeq[i].views[0].first;
		int ptId1  = trackSeq[i].views[0].second;
		int imgId2 = trackSeq[i].views[1].first;
		int ptId2  = trackSeq[i].views[1].second;		
		Point2DDouble pl,pr;

		if( imgId1==leftImageId && imgId2==rightImageId )
		{
			pl.p[0] = imageFeatures[imgId1].featPts[ptId1].cx;
			pl.p[1] = imageFeatures[imgId1].featPts[ptId1].cy;
			lpts.push_back(pl);
			pr.p[0] = imageFeatures[imgId2].featPts[ptId2].cx;
			pr.p[1] = imageFeatures[imgId2].featPts[ptId2].cy;
			rpts.push_back(pr);
		}
		else
		{
			pl.p[0] = imageFeatures[imgId2].featPts[ptId2].cx;
			pl.p[1] = imageFeatures[imgId2].featPts[ptId2].cy;
			lpts.push_back(pl);
			pr.p[0] = imageFeatures[imgId1].featPts[ptId1].cx;
			pr.p[1] = imageFeatures[imgId1].featPts[ptId1].cy;
			rpts.push_back(pr);		
		}
	}

	CRelativePoseBase* pRP = new CEstimatePose5Point();
	pRP->EstimatePose(lpts, rpts, cameras[leftImageId], cameras[rightImageId] );   

	//2.2 triangulation
	CTriangulateBase* pTri = new CTriangulateCV();
	vector<Point3DDouble> gpts;
	vector<double> errorarray;
	pTri->Triangulate(lpts, rpts, cameras[leftImageId], cameras[rightImageId], gpts, errorarray);


	//2.3 bundle adjust for pair
	vector<int> cameraIDOrder;
	cameraIDOrder.push_back(leftImageId);
	cameraIDOrder.push_back(rightImageId);
	for(int i=0; i<gpts.size(); i++)
	{
		trackSeq[i].grd = gpts[i];
		trackSeq[i].derror = errorarray[i];
	}
	WritePMVSPly("c:\\temp\\pair.ply", gpts);

#ifdef CERES_LIB
	CeresBA(trackSeq, imageFeatures, cameras);
#endif	
	
	//update the track coordinates
	CaculateTrackSeqGrd(imageFeatures, trackSeq, cameras, true);


	//3. adding new images
	while(1)
	{
		int newCameraIndex = SelectNewImage(cameraVisited, tracks, trackSeq);

		if(newCameraIndex<0)
			break;

		//calculate the camera parameters of new selected
		CalculateNewCamParas(newCameraIndex, imageFeatures, trackSeq, tracks, cameras[newCameraIndex]);
		 
		//update tracks according to the new image 
		UpdateBATracks(newCameraIndex, cameraVisited, imageFeatures, tracks, trackSeq);

		//update the track point 3D coordinate
		CaculateTrackSeqGrd(imageFeatures, trackSeq, cameras, true);
		
		cameraVisited[newCameraIndex] = 1;

#ifdef CERES_LIB
		printf("BA for all cameras... \n");
		CeresBA(trackSeq, imageFeatures, cameras);
#endif	
	}

	//update the track coordinates
	CaculateTrackSeqGrd(imageFeatures, trackSeq, cameras, true);
	
	//save the ba results: camera position, track points
	vector<Point3DDouble> goodGrds;
	printf("output the optimized track points.... \n");
	for(int i=0; i<trackSeq.size(); i++)
	{
		if( trackSeq[i].derror < 4 )
		{
			goodGrds.push_back( trackSeq[i].grd );
			for(int j=0; j<trackSeq[i].views.size(); j++)
			{
				printf("%d %d ", trackSeq[i].views[j].first, trackSeq[i].views[j].second);
			}
			printf("\n");
		}
	}
	for(int i=0; i<cameras.size(); i++)
	{
		Point3DDouble cp;
		cp.p[0] = cameras[i].t[0];
		cp.p[1] = cameras[i].t[1];
		cp.p[2] = cameras[i].t[2];
		goodGrds.push_back( cp );
	}

	WritePMVSPly("c:\\temp\\ba.ply", goodGrds);


	return 0;
}