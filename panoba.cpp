

/*
   



*/


#include"defines.hpp"
#include"panoba.hpp"
#include"relativepose.hpp"
#include"bundlerio.hpp"
#include "triangulate.hpp"

#include "panorama.hpp"
#include "ceresba.hpp"




CPanoBAWithPos::CPanoBAWithPos()
{

}

CPanoBAWithPos::~CPanoBAWithPos()
{

}

int CPanoBAWithPos::SetCenter(double ax, double ay, double az)
{
	mAx = ax;
	mAy = ay;
	mAz = az;

	return 0;
}

int CPanoBAWithPos::BundleAdjust(int numCameras, vector<CameraPara>& cameras,
	vector<ImgFeature>& imageFeatures,	
	vector<TrackInfo>& tracks)
{
	if (cameras.size() < 2)
		return -1;

	//calculate the threshold
	int projThreshold = 8;
	if (cameras[0].camtype == PerspectiveCam)
	{
		projThreshold = cameras[0].cols * 0.003; //
		projThreshold = min(16, max(8, projThreshold));
	}
	else if (cameras[0].camtype == PanoramCam)
	{
		double radius = (double)(cameras[0].cols) / (2 * PI);
		projThreshold = radius*(10.0 / 180.0*PI); //angle: degree
	}

	//because we use "trackSeq" to save the current tracks, so wo need to generate connections
	//between the original "tracks" and "trackSeq"
	for (int i = 0; i < tracks.size(); i++)
	{
		tracks[i].extra = i;
	}

	//calculate the 3D coordinates
	CaculateTrackSeqGrd(imageFeatures, tracks, cameras, true);
	/*printf("Track grd: \n");
	for (int i = 0; i<tracks.size(); i++)
	{
		printf("%lf %lf %lf  %lf \n", 
			tracks[i].grd.p[0], tracks[i].grd.p[1], tracks[i].grd.p[2],
			tracks[i].GetProjectionError());

		double height = tracks[i].GetGround().p[2] + mAz;

		if ( tracks[i].GetProjectionError() > projThreshold || height < 0)
		{
			tracks[i].Clear();
		}
	}*/


	//
	printf("\n BA for all cameras... \n");
	//optimization for all cameras
	vector<int> cameraIDOrder;
	cameraIDOrder.resize(cameras.size());
	for (int i = 0; i < cameras.size(); i++)
		cameraIDOrder[i] = i;

	SaveTracksToPly("c:\\temp\\tracks-raw_pos.ply", tracks, cameraIDOrder, cameras);


	//RefineAllParametersWithPos(tracks, imageFeatures, cameraIDOrder, cameras, mAz);
	RefineAllParameters(tracks, imageFeatures, cameraIDOrder, tracks, cameras);

	SaveTracksToPly("c:\\temp\\final-ba-withpos.ply", tracks, cameraIDOrder, cameras);

	//output the camera parameters
	printf("\n ****************** Final results: ********************** \n");
	printf("Camera Number: %d \n", cameraIDOrder.size());
	for (int i = 0; i<cameraIDOrder.size(); i++)
	{
		int camId = cameraIDOrder[i];
		printf("camera: %d ", camId);
		printf(" focal: %lf  ", cameras[camId].focalLen);
		printf(" angle: %lf %lf %lf  ", cameras[camId].ax, cameras[camId].ay, cameras[camId].az);
		printf(" position: %lf %lf %lf \n", cameras[camId].T[0], cameras[camId].T[1], cameras[camId].T[2]);
	}

	return 0;
}



/////////////////// free ba ///////////////////////////
CPanoBA::CPanoBA()
{

}

CPanoBA::~CPanoBA()
{

}

int CPanoBA::BundleAdjust(int numCameras, vector<CameraPara>& cameras, 
						  vector<ImgFeature>& imageFeatures, 
						  vector<PairMatchRes>& pairMatchs, 
						  vector<TrackInfo>& tracks)
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
		
		/*
		if( maxInlier<pairMatchs[i].inlierRatio )
		{
			maxInlier = pairMatchs[i].inlierRatio;
			index = i;
		}*/
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

	printf("initialized pair: %d - %d : \n", leftImageId, rightImageId);
	CRelativePoseBase* pRP = new CEstimatePose5PointPano();
	pRP->EstimatePose(lpts, rpts, cameras[leftImageId], cameras[rightImageId] );   


	//2.2 triangulation
	CTriangulateBase* pTri = new CTriangulateCV(); //new CTriangulatePano();
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
	SaveTracksToPly("c:\\temp\\pair-raw.ply", trackSeq, cameraIDOrder, cameras);
	
	
	RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);
  	SaveTracksToPly("c:\\temp\\pair-raw-inliers.ply", trackSeq, cameraIDOrder, cameras);
	RunBA(trackSeq, imageFeatures, cameraIDOrder, cameras);
	RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);
	//int num_pruned = RemoveBadPoints(trackSeq, imageFeatures, cameras);
	

	//RefineAllParameters(trackSeq, imageFeatures, cameraIDOrder, tracks, cameras);

	SaveTracksToPly("c:\\temp\\pair.ply", trackSeq, cameraIDOrder, cameras);

	//3. adding new images
	while(1)
	{
		int newCameraIndex = SelectNewImage(cameraVisited, tracks, trackSeq);
		if(newCameraIndex<0)
			break;

		//calculate the camera parameters of new selected
		printf("\n\n************** get initial parameters for camera: %d ******************** \n", newCameraIndex);
		int initFocalLen = cameras[cameraIDOrder[0]].focalLen;
		int res = CalculateNewCamParas(newCameraIndex, imageFeatures, trackSeq, tracks, 
			initFocalLen, cameras[newCameraIndex]);

		//DLT pose estimation failed
		if( res < 0 )
		{
			cameraVisited[newCameraIndex] = 1;
			printf("\n discard camera: %d \n", newCameraIndex);
			continue;
		}

		//adding the image points of new image into the current tracks 
		printf("\n\n\n ******************* adding image %d ******************** \n", newCameraIndex);
		UpdateBATracks(newCameraIndex, cameraVisited, cameraIDOrder, imageFeatures, tracks, trackSeq, cameras);


		printf("\n BA for all cameras... \n");
		//optimization for all cameras
		RefineAllParameters(trackSeq, imageFeatures, cameraIDOrder, tracks, cameras);

		//save temperary file
		char file[256];
		sprintf(file, "c:\\temp\\bundle_%d.ply", cameraIDOrder.size());
		SaveTracksToPly(file, trackSeq, cameraIDOrder, cameras);
	}

	SaveTracksToPly("c:\\temp\\final-ba.ply", trackSeq, cameraIDOrder, cameras);

	//output the camera parameters
	printf("\n ****************** Final results: ********************** \n");
	printf("Camera Number: %d \n", cameraIDOrder.size());
	for(int i=0; i<cameraIDOrder.size(); i++)
	{
		int camId = cameraIDOrder[i];
		printf("camera: %d ", camId);
		printf(" focal: %lf  ", cameras[camId].focalLen);
		printf(" angle: %lf %lf %lf  ", cameras[camId].ax, cameras[camId].ay, cameras[camId].az);
		printf(" position: %lf %lf %lf \n", cameras[camId].T[0], cameras[camId].T[1], cameras[camId].T[2]);
	}
	
	return 0;
}


////////////////////////////////////////////////////////////////////////////////////
CPanoGlobalBA::CPanoGlobalBA()
{
	mnProjNum = 0;
	mnGrdNum = 0;
	mnCamNum = 0;
	mpProjPt   = NULL;
	mpGrdPt    = NULL;
	mpCamPoseParas = NULL;
}

CPanoGlobalBA::~CPanoGlobalBA()
{
	delete mpProjPt;
	delete mpGrdPt;
	delete mpCamPoseParas;
}

int CPanoGlobalBA::Init(vector<CameraPara>& cameras,
	vector<ImgFeature>& imageFeatures,
	vector<TrackInfo>&  tracks)
{

	int nBACameras = cameras.size(); //the cameras attending the BA

	if (nBACameras<2)
		return -1;

	//intrinsic camera parameters
	mCamCaliParas[0] = cameras[0].focalLen;
	mCamCaliParas[1] = 0;
	mCamCaliParas[2] = 0;

	//exterior camera parameters
	mnCamNum = nBACameras;
	mpCamPoseParas = new double[nBACameras * 6];
	for (int i = 0; i<nBACameras; i++)
	{
		double aa[3];
		rot2aa(cameras[i].R, aa);
		mpCamPoseParas[i * 6] = aa[0];
		mpCamPoseParas[i * 6 + 1] = aa[1];
		mpCamPoseParas[i * 6 + 2] = aa[2];

		//from R(X-T) to RX + T
		//revised by xiedonghai, 2016.12.30
		double it[3];
		mult(cameras[i].R, cameras[i].T, it, 3, 3, 1);

		mpCamPoseParas[i * 6 + 3] = -it[0];
		mpCamPoseParas[i * 6 + 4] = -it[1];
		mpCamPoseParas[i * 6 + 5] = -it[2];
	}

	//ground point parameters
	mnGrdNum = tracks.size();
	mpGrdPt = new double[mnGrdNum * 3];
	for (int i = 0; i<mnGrdNum; i++)
	{
		Point3DDouble gp = tracks[i].GetGround();
		mpGrdPt[i * 3]     = gp.p[0];
		mpGrdPt[i * 3 + 1] = gp.p[1];
		mpGrdPt[i * 3 + 2] = gp.p[2];
	}

	//calculate the projections
	int nProjection = 0;
	int validGrdPt = 0;
	for (int i = 0; i<tracks.size(); i++)
	{
		int nview = tracks[i].GetImageKeySum(); //trackSeq[i].views.size();
		nProjection += nview;
		if (nview>0)
			validGrdPt++;
	}

	//fill projections and generate connections
	mnProjNum = nProjection;
	mVecCamIndex.resize(nProjection);
	mVecTrackIndex.resize(nProjection);
	//double* projections = new double[nProjection * 2]; 
	mpProjPt = new double[nProjection * 2];
	int ip = 0;
	for (int i = 0; i<tracks.size(); i++)
	{
		int nview = tracks[i].views.size();

		for (int j = 0; j<nview; j++)
		{
			ImageKey ik = tracks[i].GetImageKey(j);

			int cameraID = ik.first;
			int ptIndex = ik.second;

			Point2DDouble ipt = imageFeatures[cameraID].GetCenteredPt(ptIndex);
			mpProjPt[ip * 2] = ipt.p[0];
			mpProjPt[ip * 2 + 1] = ipt.p[1];

			//connections between projection and track & camera
			mVecTrackIndex[ip] = i;
			mVecCamIndex[ip]   = cameraID;

			ip++;
		}
	}

	mCameras = cameras;

	return 0;
}


int CPanoGlobalBA::AddFreePtBlock()
{

	for (int i = 0; i<mnProjNum; i++)
	{
		int cameraId = mVecCamIndex[i];

		if (mCameras[cameraId].camtype == PerspectiveCam)
		{
			ceres::CostFunction* cost_function =
				SFMReprojectionError::Create(mpProjPt[2 * i + 0], mpProjPt[2 * i + 1]);

			int trackId = mVecTrackIndex[i];

			mProblem.AddResidualBlock(cost_function,
				NULL /* squared loss */,
				mCamCaliParas,			    //inner camera parameters: 3
				mpCamPoseParas + cameraId * 6,      //outer camera parameters: 6
				mpGrdPt + trackId * 3);			    //parameters for point :   3
		}
		else if (mCameras[cameraId].camtype == PanoramCam)
		{
			double radius = (double)(mCameras[cameraId].cols) / (2 * PI);
			int trackId = mVecTrackIndex[i];

			ceres::CostFunction* cost_function =
				SFMPanoReprojectionError::Create(mpProjPt[2 * i + 0],
				mpProjPt[2 * i + 1],
				radius);

			//ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
			ceres::LossFunction* loss_function = new ceres::CauchyLoss(0.5);
			mProblem.AddResidualBlock(cost_function,
				loss_function,
				mpCamPoseParas + cameraId * 6,      //outer camera parameters: 6
				mpGrdPt + trackId * 3);			    //parameters for point : 3
		}
	}
	
	return 0;
}

int CPanoGlobalBA::AddCtrlPtBlock(vector<stCtrlPt>& ctrlPts)
{
	for (int i = 0; i < ctrlPts.size(); i++)
	{
		for (int k = 0; k < ctrlPts[i].projPts.size(); k++)
		{
			int cameraId = ctrlPts[i].camIds[k];

			double radius = (double)(mCameras[cameraId].cols) / (2 * PI);

			double px = ctrlPts[i].projPts[k].p[0];
			double py = ctrlPts[i].projPts[k].p[1];

			double gx, gy, gz;
			gx = ctrlPts[i].grd.p[0];
			gy = ctrlPts[i].grd.p[1];
			gz = ctrlPts[i].grd.p[2];

			ceres::CostFunction* cost_function =
				SFMPanoCtrlPtReprojectionError::Create(px, py, radius, gx, gy, gz);

			//ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
			//ceres::LossFunction* loss_function = new ceres::CauchyLoss(0.5);
			mProblem.AddResidualBlock(cost_function,
				NULL,
				mpCamPoseParas + cameraId * 6     //outer camera parameters: 6
				);
		}
	}

	return 0;
}

int CPanoGlobalBA::Run()
{
	//save the optimized results
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	//options.linear_solver_type = ceres::SPARSE_SCHUR;
	//options.parameter_tolerance = 1e-9;
	options.gradient_tolerance = 1e-16;
	options.function_tolerance = 1e-16;
	options.minimizer_progress_to_stdout = true;

	cout << "NumParameterBlocks: " << mProblem.NumParameterBlocks() << "\n";
	cout << "NumResidualBlocks: " << mProblem.NumResidualBlocks() << "\n";
	cout << "NumParameters: " << mProblem.NumParameters() << "\n";

	ceres::Solver::Summary summary;
	ceres::Solve(options, &mProblem, &summary);

	std::cout << summary.FullReport() << "\n";
	std::cout << "final cost: " << summary.final_cost << "\n";
	std::cout << "average error: " << (summary.final_cost / mnProjNum) << "\n";
	std::cout << "successful steps: " << summary.num_successful_steps << "\n";
	std::cout << "unsuccessful steps: " << summary.num_unsuccessful_steps << "\n";

	std::cout << "all cameras optimized results: " << "\n";

	return 0;
}

int CPanoGlobalBA::Output(vector<CameraPara>& cameras, vector<TrackInfo>&  tracks)
{
	for (int i = 0; i<cameras.size(); i++)
	{
		int camId = i;

		cameras[camId].focalLen = mCamCaliParas[0];
		cameras[camId].k1 = mCamCaliParas[1];
		cameras[camId].k2 = mCamCaliParas[2];

		//from axis-angle to rotation matrix
		double aa[3];
		aa[0] = mpCamPoseParas[i * 6];
		aa[1] = mpCamPoseParas[i * 6 + 1];
		aa[2] = mpCamPoseParas[i * 6 + 2];

		double Rt[9];
		aa2rot(aa, Rt);
		double R[9];
		transpose(Rt, R, 3, 3);
		memcpy(cameras[camId].R, R, sizeof(double)* 9);
		double ea[3];
		rot2eular(R, ea);

		cameras[camId].ax = ea[0];
		cameras[camId].ay = ea[1];
		cameras[camId].az = ea[2];

		//cout<<"Angle: "<< cameras[camId].ax <<" "<< 
		//cameras[camId].ay<<" " << cameras[camId].az << "\n"; 

		//from RX+T to R(X-T)
		double iR[9];
		memcpy(iR, R, sizeof(double)* 9);
		invers_matrix(iR, 3);
		double it[3];
		it[0] = mpCamPoseParas[i * 6 + 3];
		it[1] = mpCamPoseParas[i * 6 + 4];
		it[2] = mpCamPoseParas[i * 6 + 5];
		double et[3];
		mult(iR, it, et, 3, 3, 1);

		cameras[camId].T[0] = -et[0];
		cameras[camId].T[1] = -et[1];
		cameras[camId].T[2] = -et[2];

	}
	std::cout << endl;

	for (int i = 0; i<mnGrdNum; i++)
	{
		Point3DDouble gp;
		gp.p[0] = mpGrdPt[i * 3];
		gp.p[1] = mpGrdPt[i * 3 + 1];
		gp.p[2] = mpGrdPt[i * 3 + 2];

		tracks[i].SetGround(gp);
	}


	return 0;
}