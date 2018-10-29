
#include "ceresba.hpp"


//#include "float.h"
#include "math.h"

#include "ba.hpp"
#include "cali.hpp"
#include "relativepose.hpp"
#include "distortion.hpp"
#include "bundlerio.hpp"
#include "rotation.hpp"
#include "badata.hpp"
#include "panorama.hpp"
#include "triangulate.hpp"

//baselib
#include "baselib.h"

//corelib matrix
#include "Matrix.h"

// find the track index in the current track sequence:  
// image feature --> all tracks --> current track seq
vector<int> GetTracksInImage(ImgFeature& imageFeat, vector<TrackInfo>& allTracks)
{
	vector<int>  tracksSeenByCurrentImage;
	for (int i = 0; i < imageFeat.GetFeatPtSum(); i++)
	{

		int trackid = imageFeat.GetPtTrackIndex(i);
		if (trackid < 0)
			continue;

		int currentTrackId = allTracks[trackid].GetRemapTrackId();
		if (currentTrackId < 0)
			continue;

		tracksSeenByCurrentImage.push_back(currentTrackId);
	}

	return tracksSeenByCurrentImage;
}



//find the tracks seen by the "camid"
vector<TrackInfo> FindLocalTracks(int camid, vector<TrackInfo>& curTracks,
	vector<TrackInfo>& allTracks, vector<ImgFeature>& imageFeatures)
{
	vector<TrackInfo> localTracks;


	//find the tracks seen by the "camid"
	vector<int>  tracksSeenByCurrentImage = GetTracksInImage(imageFeatures[camid], allTracks);


	//find the images seeing tracks
	vector<int> cameraid;
	for (int i = 0; i < tracksSeenByCurrentImage.size(); i++)
	{
		int tid = tracksSeenByCurrentImage[i];

		//search the images seeing the tracks
		for (int k = 0; k < curTracks[tid].GetImageKeySum(); k++)
		{
			ImageKey ik = curTracks[tid].GetImageKey(k);

			int imageid = ik.first;
			int featid = ik.second;

			cameraid.push_back(imageid);
		}
	}
	//delete the repeated number
	sort(cameraid.begin(), cameraid.end());
	vector<int>::iterator pos;
	pos = unique(cameraid.begin(), cameraid.end());
	cameraid.erase(pos, cameraid.end());


	//collect the tracks seen by the neighbor images
	vector<int> neighborTrackIds;
	for (int i = 0; i < cameraid.size(); i++)
	{
		int cid = cameraid[i];

		vector<int> imageTracks = GetTracksInImage(imageFeatures[cid], allTracks);

		for (int k = 0; k < imageTracks.size(); k++)
			neighborTrackIds.push_back(imageTracks[k]);
	}
	//delete the repeated number
	sort(neighborTrackIds.begin(), neighborTrackIds.end());
	pos = unique(neighborTrackIds.begin(), neighborTrackIds.end());
	neighborTrackIds.erase(pos, neighborTrackIds.end());


	//copy the tracks
	for (int i = 0; i < neighborTrackIds.size(); i++)
	{
		int idInCurrentTrackSeq = neighborTrackIds[i];
		TrackInfo ti = curTracks[idInCurrentTrackSeq];
		ti.extra = idInCurrentTrackSeq;
		localTracks.push_back(ti);
	}

	return localTracks;
}

//collect the camera index in the tracks
vector<int> GetCameraInTracks(vector<TrackInfo>& tracks)
{
	vector<int> cameraid;

	for (int i = 0; i < tracks.size(); i++)
	{
		//determine if the track is seen by "camid"
		for (int k = 0; k < tracks[i].GetImageKeySum(); k++)
		{
			ImageKey ik = tracks[i].GetImageKey(k);

			int imageid = ik.first;
			int featid = ik.second;

			cameraid.push_back(imageid);
		}
	}

	//delete the repeated number
	sort(cameraid.begin(), cameraid.end());
	vector<int>::iterator pos;
	pos = unique(cameraid.begin(), cameraid.end());
	cameraid.erase(pos, cameraid.end());

	return cameraid;
}

//////////////////////////////////////// new apis/////////////////////////////////////////////
/*
   trackseq:       the track points for BA
   imageFeatures:  the image points 
   cameraIDOrder:  cameras to be adjusted
   cameras:        all camera parameters
*/
int RunBA( vector<TrackInfo>& trackSeq, vector<ImgFeature>& imageFeatures, 
		   vector<int> cameraIDOrder, vector<CameraPara> &cameras)
{
	int nBACameras = cameraIDOrder.size(); //the cameras attending the BA

	if(nBACameras<2)
		return -1;

	int nAllCameras = cameras.size();
	vector<int> cameraMap;
	cameraMap.resize(nAllCameras, 0);      //save the index of reordered cameras 
	
	//intrinsic camera parameters
	double* pInteriorParams = new double[3];          //focal length, k1, k2
	pInteriorParams[0] = cameras[cameraIDOrder[0]].focalLen;
	pInteriorParams[1] = 0;
	pInteriorParams[2] = 0;	

	//extrinsic camera parameters
	double* pOuterParams = new double[nBACameras*6]; //axis-angle vector, t0,t1,t2
	for(int i=0; i<nBACameras; i++)
	{		
		int ci = cameraIDOrder[i]; //camera index
		cameraMap[ci] = i;

		double aa[3];
		rot2aa(cameras[ci].R, aa);
		pOuterParams[i*6]   = aa[0];
		pOuterParams[i*6+1] = aa[1];
		pOuterParams[i*6+2] = aa[2];
		
		//from R(X-T) to RX + T
		//revised by xiedonghai, 2016.12.30
		double it[3];
		mult(cameras[ci].R, cameras[ci].T, it, 3, 3, 1);

		pOuterParams[i*6+3] = -it[0];
		pOuterParams[i*6+4] = -it[1];
		pOuterParams[i*6+5] = -it[2];		
	}
	
	//ground point parameters
	int num_pts = trackSeq.size();
	double* grdPt = new double[num_pts*3];
	for(int i=0; i<num_pts; i++)
	{   
		Point3DDouble gp = trackSeq[i].GetGround();
		grdPt[i*3]   = gp.p[0];
		grdPt[i*3+1] = gp.p[1];
		grdPt[i*3+2] = gp.p[2];
	}

	//generate observation points
	int nProjection = 0;
	int validGrdPt = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{
		int nview = trackSeq[i].GetImageKeySum(); //trackSeq[i].views.size();
		nProjection += nview;
		if(nview>0)
			validGrdPt++;
	}
	vector<int> vecCamIndex;    // camera index for each projection point
	vector<int> vecTrackIndex;  // track point index for each projection point
	vecCamIndex.resize(nProjection);
	vecTrackIndex.resize(nProjection);
	double* projections = new double[nProjection*2]; //(double*)malloc(nProjection*2*sizeof(double));
	int ip = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{   
		int nview = trackSeq[i].views.size();

		for(int j=0; j<nview; j++)
		{
			ImageKey ik = trackSeq[i].GetImageKey(j);

			int cameraID = ik.first;
			int ptIndex  = ik.second;

			Point2DDouble ipt = imageFeatures[cameraID].GetCenteredPt(ptIndex);
			projections[ip*2]   = ipt.p[0]; 
			projections[ip*2+1] = ipt.p[1]; 			

			vecTrackIndex[ip]  = i; 

			int mapId = cameraMap[cameraID];
			vecCamIndex[ip]    = mapId;

			ip++;
		}
	}
	
	printf("3D Points: %d  valid points: %d  projection number: %d \n", num_pts, validGrdPt, nProjection);
	//google::InitGoogleLogging(NULL);
	ceres::Problem problem;
	//invoke ceres functions, each time adding one projection
	for(int i=0; i<nProjection; i++)
	{
		int cameraId = vecCamIndex[i];
		
		if( cameras[cameraId].camtype == PerspectiveCam )
		{
			ceres::CostFunction* cost_function =
				SFMReprojectionError::Create(projections[2 * i + 0], projections[2 * i + 1]);

			int trackId  = vecTrackIndex[i];
			
			ceres::LossFunction* loss_function = new ceres::CauchyLoss(0.5);
			
			problem.AddResidualBlock(cost_function,
				//NULL /* squared loss */,
				loss_function,
				pInteriorParams,			    //inner camera parameters: 3
				pOuterParams + cameraId*6,      //outer camera parameters: 6
				grdPt + trackId*3);			    //parameters for point :   3
		}
		else if(cameras[cameraId].camtype == PanoramCam)
		{
			double radius = (double)(cameras[cameraId].cols) / (2*PI); 
			
			//if (cameraId == 0 || cameraId == nBACameras - 1)
			//{
			//	
			//	double a0 = pOuterParams[cameraId * 6];
			//	double a1 = pOuterParams[cameraId * 6+1];
			//	double a2 = pOuterParams[cameraId * 6+2];
			//	double tx = pOuterParams[cameraId * 6+3];
			//	double ty = pOuterParams[cameraId * 6+4];
			//	double tz = pOuterParams[cameraId * 6+5];

			//	ceres::CostFunction* cost_function =
			//		RefinePanoErrorFixedCamera::Create(projections[2 * i + 0],
			//		projections[2 * i + 1],
			//		a0,a1,a2,tx,ty,tz,
			//		radius);
			//	int trackId = vecTrackIndex[i];
			//	problem.AddResidualBlock(cost_function,
			//		NULL /* squared loss */,
			//		grdPt + trackId * 3);			    //parameters for point : 3
			//}
			//else




			int trackId = vecTrackIndex[i];

			if (trackSeq[trackId].IsCtrl())
			{
				
				ceres::CostFunction* cost_function =
					SFMPanoCtrlPtReprojectionError::Create(projections[2 * i + 0],
					projections[2 * i + 1],
					radius,
					grdPt[trackId * 3],
					grdPt[trackId * 3 + 1],
					grdPt[trackId * 3 + 2]);

				//ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
				//ceres::LossFunction* loss_function = new ceres::CauchyLoss(0.5);
				problem.AddResidualBlock(cost_function,
					NULL,
					//loss_function,
					pOuterParams + cameraId * 6     //outer camera parameters: 6
					);
			}
			else
			{
				ceres::CostFunction* cost_function =
					SFMPanoReprojectionError::Create(projections[2 * i + 0],
					projections[2 * i + 1],
					radius);

				//ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
				ceres::LossFunction* loss_function = new ceres::CauchyLoss(0.5);
				problem.AddResidualBlock(cost_function,
					loss_function,
					pOuterParams + cameraId * 6,      //outer camera parameters: 6
					grdPt + trackId * 3);			    //parameters for point : 3
			}

			//ceres::CostFunction* cost_function =
			//	SFMPanoErrorFixedAngle::Create(projections[2 * i + 0],
			//	projections[2 * i + 1],
			//	pOuterParams[cameraId*6],
			//	pOuterParams[cameraId*6+1],
			//	pOuterParams[cameraId*6+2],
			//	radius);

			//int trackId = vecTrackIndex[i];
			//problem.AddResidualBlock(cost_function,
			//	NULL, // squared loss
			//	pOuterParams + cameraId * 6 + 3,      //exterior position camera parameters: 3
			//	grdPt + trackId * 3);			      //parameters for point :   3
			//
				
		}
	}

	//save the optimized results
	ceres::Solver::Options options;
	//options.linear_solver_type = ceres::DENSE_SCHUR;
	options.linear_solver_type = ceres::SPARSE_SCHUR;
	//options.parameter_tolerance = 1e-9;
	//options.gradient_tolerance = 1e-16;
	//options.function_tolerance = 1e-16;
	options.minimizer_progress_to_stdout = true;

	cout<<"NumParameterBlocks: "<<problem.NumParameterBlocks()<<"\n";
	cout<<"NumResidualBlocks: "<<problem.NumResidualBlocks()<<"\n";
	cout<<"NumParameters: "<<problem.NumParameters()<<"\n";

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	//std::cout << summary.FullReport() << "\n";
	std::cout<<"final cost: "<<summary.final_cost<<"\n";
	std::cout<<"average error: "<<(summary.final_cost/nProjection)<<"\n";
	std::cout<<"successful steps: "<<summary.num_successful_steps<<"\n";
	std::cout<<"unsuccessful steps: "<<summary.num_unsuccessful_steps<<"\n";
	

	//if the error is large, delete the new added data
	

	std::cout<<"all cameras optimized results: "<<"\n";

	printf("Intrinsic parameters... \n");
	for(int i=0; i<3; i++)
		std::cout<<pInteriorParams[i]<<"  ";
	std::cout<<endl;
	
	
	printf("Extrinsic parameters ...\n");
	for(int i=0; i<nBACameras; i++)
	{
		int camId = cameraIDOrder[i];
		int mapId = cameraMap[camId];

		cameras[camId].focalLen = pInteriorParams[0];
		cameras[camId].k1    = pInteriorParams[1];
		cameras[camId].k2    = pInteriorParams[2];
		
		//from axis-angle to rotation matrix
		double aa[3];
		aa[0] = pOuterParams[mapId*6];
		aa[1] = pOuterParams[mapId*6+1];
		aa[2] = pOuterParams[mapId*6+2];

		double Rt[9];
		aa2rot(aa, Rt);
		double R[9];
		transpose(Rt, R, 3, 3);
		memcpy(cameras[camId].R, R, sizeof(double)*9);
		double ea[3];
		rot2eular(R, ea);

		cameras[camId].ax = ea[0];
		cameras[camId].ay = ea[1];
		cameras[camId].az = ea[2];
		
		//cout<<"Angle: "<< cameras[camId].ax <<" "<< 
		//cameras[camId].ay<<" " << cameras[camId].az << "\n"; 

		//from RX+T to R(X-T)
		double iR[9];
		memcpy(iR, R, sizeof(double)*9);
		invers_matrix(iR, 3);
		double it[3];
		it[0] = pOuterParams[mapId*6+3];
		it[1] = pOuterParams[mapId*6+4];
		it[2] = pOuterParams[mapId*6+5];
		double et[3];
		mult(iR, it, et, 3, 3, 1);

		cameras[camId].T[0] = -et[0];
		cameras[camId].T[1] = -et[1];
		cameras[camId].T[2] = -et[2];

		cameras[camId].xs = -et[0];
		cameras[camId].ys = -et[1];
		cameras[camId].zs = -et[2];

		//cout<<"Translation: "<<cameras[camId].t[0]<<" "<<
			//cameras[camId].t[1]<<" "<<cameras[camId].t[2]<<"\n";

		//std::cout<<endl;
	}
	std::cout<<endl;

	for(int i=0; i<num_pts; i++)
	{   
		Point3DDouble gp;
		gp.p[0] = grdPt[i*3];
		gp.p[1] = grdPt[i*3+1];
		gp.p[2] = grdPt[i*3+2];

		trackSeq[i].SetGround(gp);		
	}

	return 0;
}


/* remove the track whose projection error is large
*/
int RemoveOutlierPtsWithPos(vector<TrackInfo>& tracks,
	vector<ImgFeature>& imageFeatures, vector<int> cameraIDOrder,
	vector<CameraPara>& cameras, double az)
{
	int numCam = cameraIDOrder.size();

	vector<int> outliers; //save the index of outliers for ground point
	vector<int> inliers;
	for (int i = 0; i<numCam; i++)
	{
		int camId = cameraIDOrder[i];
		int nFeatPt = imageFeatures[camId].GetFeatPtSum();

		double radius = (double)(cameras[camId].cols) / (2 * PI);

		vector<double> dists;
		for (int k = 0; k<nFeatPt; k++)
		{
			//the projection point corresponds to no track
			int nTrackId = imageFeatures[camId].GetPtTrackIndex(k);
			if (nTrackId<0)
				continue;

			//the projection point has not been added to the current track seq
			int nCurrentTrack = tracks[nTrackId].GetRemapTrackId();
			if (nCurrentTrack<0)
				continue;

			//the track has been discarded
			if (tracks[nCurrentTrack].GetImageKeySum() == 0)
				continue;

			Point3DDouble gp = tracks[nCurrentTrack].GetGround();
			Point2DDouble ip = imageFeatures[camId].GetCenteredPt(k);

			double dist = 0;
			//if(cameras[camId].camtype == PerspectiveCam)
			{
				Point2DDouble projPt;
				GrdToImg(gp, projPt, cameras[camId]);
				double dx = projPt.p[0] - ip.p[0];
				double dy = projPt.p[1] - ip.p[1];
				dist = sqrt(dx*dx + dy*dy);
			}
			/*else if(cameras[camId].camtype == PanoramCam)
			{
			double gx,gy,gz;
			SphereTo3D_center(ip.p[0], ip.p[1], radius, gx, gy, gz);
			double n1 = sqrt(gp.p[0]*gp.p[0] + gp.p[1]*gp.p[1] + gp.p[2]*gp.p[2]) + 0.000001;
			double n2 = sqrt(gx*gx + gy*gy + gz*gz) + 0.000001;
			double dotv = gp.p[0]*gx + gp.p[1]*gy + gp.p[2]*gz;
			double angle = acos( dotv / (n1*n2) );
			dist = angle*radius;
			}
			*/

			dists.push_back(dist);
		}

		//Estimate the median of the distances and calculate the threshold
		double *pDist = new double[dists.size()];
		for (int i = 0; i<dists.size(); i++)
			pDist[i] = dists[i];
		double med = dll_kth_element_copy(dists.size(),
			dll_iround(0.85 * dists.size()),
			pDist);
		dll_median_copy(dists.size(), pDist);
		double thresh = 4.8 * NUM_STDDEV * med; //1.2 * NUM_STDDEV * med;  // k * stddev

		if (cameras[camId].camtype == PerspectiveCam)
		{
			thresh = CLAMP(thresh, MIN_PROJ_ERROR_THRESHOLD, MAX_PROJ_ERROR_THRESHOLD);
		}
		else if (cameras[camId].camtype == PanoramCam)
		{
			//double radius = (double)(cameras[camId].cols) / (2*PI);
			double minT = radius*0.4 / 180.0*3.14;
			double maxT = radius*5.0 / 180.0*3.14;
			printf("minT: %lf  maxT: %lf  thresh: %lf \n", minT, maxT, thresh);
			thresh = CLAMP(thresh, minT, maxT);
		}

		delete[] pDist;

		//collect the outliers
		int nIndex = 0;
		for (int k = 0; k<nFeatPt; k++)
		{
			int nTrackId = imageFeatures[camId].GetPtTrackIndex(k);
			if (nTrackId<0)
				continue;

			int nCurrentTrack = tracks[nTrackId].GetRemapTrackId();
			if (nCurrentTrack<0)
				continue;

			//the track has been discarded
			if (tracks[nCurrentTrack].GetImageKeySum() == 0)
				continue;

			if (dists[nIndex]>thresh)
			{
				bool bFound = false;
				for (int ki = 0; ki<outliers.size(); ki++)
				{
					if (nCurrentTrack == outliers[ki])
						bFound = true;
				}
				if (!bFound)
				{
					outliers.push_back(nCurrentTrack);

					//cut the connection between track and projection point
					imageFeatures[camId].SetTrackIndex(k, -1);
				}
			}
			nIndex++;
		}
	}

	//remove outliers
	//vector<int> outlierMask;
	//outlierMask.resize( trackSeq.size(), 0);

	for (int i = 0; i<outliers.size(); i++)
	{
		//do not remove the control points
		if (tracks[outliers[i]].IsCtrl())
			continue;

		//outlierMask[ outliers[i] ] = 1;
		tracks[outliers[i]].Clear();
	}

	//remove outliers according height
	for (int i = 0; i < tracks.size(); i++)
	{
		//do not remove the control points
		if (tracks[i].IsCtrl())
			continue;

		if (tracks[i].GetImageKeySum() == 0)
			continue;

		double gz = tracks[i].GetGround().p[2];

		//if (gz > 60 || (gz+az)<0 )
		if ((gz+az)<0)
			tracks[i].Clear();
	}

	/*
	vector<TrackInfo> goodTracks;
	for(int i=0; i<trackSeq.size(); i++)
	{
	if(outlierMask[i]>0)
	continue;
	goodTracks.push_back(trackSeq[i]);
	}
	trackSeq = goodTracks;
	*/

	//remove the tracks whose intersection angle is two small or large


	printf("remove %d outliers \n", outliers.size());

	return outliers.size();
}


/* remove the track whose projection error is large
*/
int RemoveOutlierPts( vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq, 
	vector<ImgFeature>& imageFeatures, vector<int> cameraIDOrder, 
	vector<CameraPara>& cameras)
{
	int numCam = cameraIDOrder.size();

	
	vector<int> outliers; //save the index of outliers for ground point
	vector<int> inliers; 
	for(int i=0; i<numCam; i++)
	{
		int camId   = cameraIDOrder[i];
		int nFeatPt = imageFeatures[camId].GetFeatPtSum();

		double radius = (double)(cameras[camId].cols) / (2*PI);

		vector<double> dists;
		for(int k=0; k<nFeatPt; k++)
		{
			//the projection point corresponds to no track
			int nTrackId = imageFeatures[camId].GetPtTrackIndex(k);
			if(nTrackId<0)
				continue;
			
			//the projection point has not been added to the current track seq
			int nCurrentTrack = tracks[nTrackId].GetRemapTrackId();
			if(nCurrentTrack<0)
				continue;

			//the track has been discarded
			if( trackSeq[nCurrentTrack].GetImageKeySum()==0 )
				continue;

			Point3DDouble gp = trackSeq[nCurrentTrack].GetGround();
			Point2DDouble ip = imageFeatures[camId].GetCenteredPt(k);

			double dist = 0;
			//if(cameras[camId].camtype == PerspectiveCam)
			{		
				Point2DDouble projPt;
				GrdToImg(gp, projPt, cameras[camId]);
				double dx = projPt.p[0] - ip.p[0];
				double dy = projPt.p[1] - ip.p[1];
				dist = sqrt(dx*dx + dy*dy);
			}
			//else if(cameras[camId].camtype == PanoramCam)
			//{
			//	double gx,gy,gz;
			//	SphereTo3D_center(ip.p[0], ip.p[1], radius, gx, gy, gz);
			//	double n1 = sqrt(gp.p[0]*gp.p[0] + gp.p[1]*gp.p[1] + gp.p[2]*gp.p[2]) + 0.000001;
			//	double n2 = sqrt(gx*gx + gy*gy + gz*gz) + 0.000001;
			//	double dotv = gp.p[0]*gx + gp.p[1]*gy + gp.p[2]*gz;
			//	double angle = acos( dotv / (n1*n2) );
			//	dist = angle*radius;
			//}
			
			
			dists.push_back(dist);
		}
		
		//Estimate the median of the distances and calculate the threshold
		double *pDist = new double[dists.size()];
		for(int i=0; i<dists.size(); i++)
			pDist[i] = dists[i];
		double med = dll_kth_element_copy(dists.size(), 
			dll_iround(0.85 * dists.size()),
			pDist);
		dll_median_copy(dists.size(), pDist);
		double thresh = 1.2 * NUM_STDDEV * med;  // k * stddev

		if(cameras[camId].camtype==PerspectiveCam)
		{
			thresh = CLAMP(thresh, MIN_PROJ_ERROR_THRESHOLD, MAX_PROJ_ERROR_THRESHOLD); 
		}
		else if(cameras[camId].camtype==PanoramCam)
		{ 
			//double radius = (double)(cameras[camId].cols) / (2*PI);
			double minT = radius*0.1/180.0*3.14;
			double maxT = radius*5.0/180.0*3.14;
			printf("minT: %lf  maxT: %lf  thresh: %lf \n", minT, maxT, thresh);
			thresh = CLAMP(thresh, minT, maxT);
		}

		delete[] pDist;

		//collect the outliers
		int nIndex = 0;
		for(int k=0; k<nFeatPt; k++)
		{
			int nTrackId = imageFeatures[camId].GetPtTrackIndex(k);
			if(nTrackId<0)
				continue;
		
			int nCurrentTrack = tracks[nTrackId].GetRemapTrackId();
			if(nCurrentTrack<0)
				continue;

			//the track has been discarded
			if( trackSeq[nCurrentTrack].GetImageKeySum()==0 )
				continue;

			if(dists[nIndex]>thresh)
			{					
				bool bFound = false;
				for(int ki=0; ki<outliers.size(); ki++)
				{
					if( nCurrentTrack == outliers[ki] )
						bFound = true;
				}
				if(!bFound)
				{
					outliers.push_back(nCurrentTrack);

					//cut the connection between track and projection point
					imageFeatures[camId].SetTrackIndex(k, -1); 
				}
			}
			nIndex ++;
		}
	}

	//remove outliers
	//vector<int> outlierMask;
	//outlierMask.resize( trackSeq.size(), 0);
	
	for(int i=0; i<outliers.size(); i++)
	{
		//outlierMask[ outliers[i] ] = 1;
		trackSeq[ outliers[i] ].Clear();
	}

	////remove outliers according height
	for (int i = 0; i < trackSeq.size(); i++)
	{
		if (trackSeq[i].GetImageKeySum() == 0)
			continue;

		double gz = trackSeq[i].GetGround().p[2];
		
		if (cameras[0].camtype == PerspectiveCam){
			if (gz>0) trackSeq[i].Clear(); //in perspective, all points should be blow to the cameras
		}			
	}

	/*
	vector<TrackInfo> goodTracks;
	for(int i=0; i<trackSeq.size(); i++)
	{
		if(outlierMask[i]>0)
			continue;
		goodTracks.push_back(trackSeq[i]);
	}
	trackSeq = goodTracks;
	*/

	//remove the tracks whose intersection angle is two small or large


	//printf("remove %d outliers \n", outliers.size());

	//return outliers.size();

	return 0;
}

/*
	remove the tracks whose angles between projection points are small
*/
int RemoveBadPoints(vector<TrackInfo>& trackSeq, 
	vector<ImgFeature>& imageFeatures, vector<CameraPara>& cameras )
{
	//CTriangulateBase* pTriangulate = new CTriangulateCV();
	
	int num_pruned = 0;

	for(int i=0; i<trackSeq.size(); i++)
	{
		bool estimate_distortion = true;

		int num_views = (int) trackSeq[i].GetImageKeySum();

		Point3DDouble gps = trackSeq[i].GetGround();

		//calculate the angle between the track point and camera center
		double maxangle = 0;
		for (int j = 0; j < num_views; j++) 
		{
			ImageKey ik1 = trackSeq[i].GetImageKey(j);
			int camera_idx = ik1.first;

			Point3DDouble camCenter;
			camCenter.p[0] = cameras[camera_idx].T[0];
			camCenter.p[1] = cameras[camera_idx].T[1];
			camCenter.p[2] = cameras[camera_idx].T[2];

			Point3DDouble v1;
			v1.p[0] = gps.p[0] - camCenter.p[0];
			v1.p[1] = gps.p[1] - camCenter.p[1];
			v1.p[2] = gps.p[2] - camCenter.p[2];

			double norm = dll_matrix_norm(3, 1, v1.p);
			v1.p[0] /= norm;
			v1.p[1] /= norm;
			v1.p[2] /= norm;

			for (int k = j+1; k < num_views; k++) 
			{
				ImageKey ik2 = trackSeq[i].GetImageKey(k);
				int camera_idx = ik2.first;

				Point3DDouble camCenter;
				camCenter.p[0] = cameras[camera_idx].T[0];
				camCenter.p[1] = cameras[camera_idx].T[1];
				camCenter.p[2] = cameras[camera_idx].T[2];

				Point3DDouble v2;
				v2.p[0] = gps.p[0] - camCenter.p[0];
				v2.p[1] = gps.p[1] - camCenter.p[1];
				v2.p[2] = gps.p[2] - camCenter.p[2];

				double norm = dll_matrix_norm(3, 1, v2.p);
				v2.p[0] /= norm;
				v2.p[1] /= norm;
				v2.p[2] /= norm;

				double dotvalue = v1.p[0]*v2.p[0] + v1.p[1]*v2.p[1] + v1.p[2]*v2.p[2]; 
				double angle = acos(dotvalue);
				if(angle>maxangle)
					maxangle = angle;
			}
		}

		//when the max angle is less than the threshold, then do not use it in the optimization
		maxangle = maxangle / PI * 180;
		if(maxangle<5)
		{
			//cut the connections between the track and projection points
			int num_views = (int) trackSeq[i].GetImageKeySum();
			for(int ki=0; ki<num_views; ki++)
			{
				ImageKey ik = trackSeq[i].GetImageKey(ki);
				int camId = ik.first;
				int ptId  = ik.second;
				imageFeatures[camId].SetTrackIndex(ptId, -1);
			}

			trackSeq[i].Clear();
			num_pruned ++;
		}
	}

	return num_pruned;
}



//////////////////////////////////////////////////////////////////////////////////////////////////




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

int CalculateCtrlPt(vector<stCtrlPt>& ctrlpts,
	vector<CameraPara>&  cameras,
	bool explicit_camera_centers){

	CTriangulateBase* pTriangulate = new CTriangulateCV();

	for (int i = 0; i<ctrlpts.size(); i++)
	{
		bool estimate_distortion = true;

		int num_views = ctrlpts[i].projPts.size(); //(int) tracks[i].views.size();

		if (num_views<2)
			continue;

		vector<Point2DDouble> pts;
		vector<CameraPara> cams;
		Point3DDouble gps;
		double ferror;

		//collect all projection points of one track
		for (int j = 0; j < num_views; j++)
		{
			int camera_idx = ctrlpts[i].camIds[j];
			
			double cx = ctrlpts[i].projPts[j].p[0];
			double cy = ctrlpts[i].projPts[j].p[1];

			Point2DDouble pt;
			pt.x = cx;
			pt.y = cy;
			pt.p[0] = cx;
			pt.p[1] = cy;
			pts.push_back(pt);

			cams.push_back(cameras[camera_idx]);
		}

		pTriangulate->Triangulate(pts, cams, gps, true, ferror);

		//ferror is the projection error of all points average in one image
		ctrlpts[i].derror = ferror;
		ctrlpts[i].grd = gps;
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
		if (tracks[i].IsCtrl())
			continue;

		bool estimate_distortion = true;

		int num_views = tracks[i].GetImageKeySum(); //(int) tracks[i].views.size();

		if(num_views<2)
			continue;

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

		//ferror is the projection error of all points average in one image
		tracks[i].derror = ferror;
		tracks[i].grd    = gps;
		
		/*
		if( gps.p[2]>0 )
		{
			printf("%lf \n", gps.p[2]);
		}*/

		/*
		//calculate the angle between the track point and camera center
		double maxangle = 0;
		for (int j = 0; j < num_views; j++) 
		{
			int camera_idx = tracks[i].views[j].first;

			Point3DDouble camCenter;
			camCenter.p[0] = cameras[camera_idx].t[0];
			camCenter.p[1] = cameras[camera_idx].t[1];
			camCenter.p[2] = cameras[camera_idx].t[2];

			Point3DDouble v1;
			v1.p[0] = gps.p[0] - camCenter.p[0];
			v1.p[1] = gps.p[1] - camCenter.p[1];
			v1.p[2] = gps.p[2] - camCenter.p[2];

			double norm = dll_matrix_norm(3, 1, v1.p);
			v1.p[0] /= norm;
			v1.p[1] /= norm;
			v1.p[2] /= norm;

			for (int k = j+1; k < num_views; k++) 
			{
				int camera_idx = tracks[i].views[k].first;

				Point3DDouble camCenter;
				camCenter.p[0] = cameras[camera_idx].t[0];
				camCenter.p[1] = cameras[camera_idx].t[1];
				camCenter.p[2] = cameras[camera_idx].t[2];

				Point3DDouble v2;
				v2.p[0] = gps.p[0] - camCenter.p[0];
				v2.p[1] = gps.p[1] - camCenter.p[1];
				v2.p[2] = gps.p[2] - camCenter.p[2];

				double norm = dll_matrix_norm(3, 1, v2.p);
				v2.p[0] /= norm;
				v2.p[1] /= norm;
				v2.p[2] /= norm;

				double dotvalue = v1.p[0]*v2.p[0] + v1.p[1]*v2.p[1] + v1.p[2]*v2.p[2]; 
				double angle = acos(dotvalue);
				if(angle>maxangle)
					maxangle = angle;
			}
		}

		//when the max angle is less than the threshold, then do not use it in the optimization
		maxangle = maxangle / PI * 180;
		if(maxangle<5)
		{
			tracks[i].valid = 0;
			//printf("invald point: %lf %lf %lf   error:%lf \n", 
			//	tracks[i].grd.p[0], tracks[i].grd.p[1], tracks[i].grd.p[2],
			//	tracks[i].derror);
		}
		*/
	}

	return 0;
}

//update the current track sequence
int UpdateBATracks( int newCameraIndex, 
					vector<int>&  cameraVisited,
					vector<int>&  cameraIDOrder,
					vector<ImgFeature>&  imageFeatures, 
					vector<TrackInfo>& tracks, 
					vector<TrackInfo>& trackSeqNew,
					vector<CameraPara>& cameras)
{

	int projThreshold = 8; 

	if( cameras[newCameraIndex].camtype == PerspectiveCam )
	{
		projThreshold = cameras[newCameraIndex].cols * 0.003 ; //
		projThreshold = min(16, max(8, projThreshold));
	}
	else if ( cameras[newCameraIndex].camtype == PanoramCam )
	{
		double radius = (double)(cameras[newCameraIndex].cols) / (2*PI);
		projThreshold = radius*(5.0/180.0*PI); //angle: degree
	}

	int nOldTracks = trackSeqNew.size();

	int nFeatPtNum = imageFeatures[newCameraIndex].featPts.size();

	for(int i=0; i<nFeatPtNum; i++)
	{
		int nTrackIndex  = imageFeatures[newCameraIndex].GetPtTrackIndex(i);

		if(nTrackIndex<0)
			continue;

		int nNewTrackIndex = tracks[nTrackIndex].extra;

		ImageKey ik;
		ik.first  = newCameraIndex;
		ik.second = i; 

		if(nNewTrackIndex>=0) //corresponding to the new Track sequence
		{
			int nviews = trackSeqNew[nNewTrackIndex].GetImageKeySum();
			
			//only valid track can be added 
			if(nviews>0)
			{
				//calculate the projection error
				Point2DDouble projPt;
				GrdToImg( trackSeqNew[nNewTrackIndex].GetGround(), projPt, cameras[newCameraIndex] );
				
				Point2DDouble ip = imageFeatures[newCameraIndex].GetCenteredPt(i);
				double dx = (projPt.p[0]-ip.p[0]);
				double dy = (projPt.p[1]-ip.p[1]);
				double error = sqrt( dx*dx+dy*dy );
				
				if(error<projThreshold)
					trackSeqNew[nNewTrackIndex].AddImageKey(ik);
			}
		}
		else //adding a new track into the new Track sequence
		{
			TrackInfo newTrack;
			newTrack.extra = nTrackIndex;
			newTrack.valid = 1;
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


	//evaluate the newly added tracks, and delete those with large projection errors
	vector<TrackInfo> newTracks;
	for(int i=nOldTracks; i<trackSeqNew.size(); i++)
	{
		newTracks.push_back( trackSeqNew[i] );
	}
	CaculateTrackSeqGrd(imageFeatures, newTracks, cameras, true);
	for(int i=nOldTracks; i<trackSeqNew.size(); i++)
	{
		if( newTracks[i-nOldTracks].GetProjectionError() > projThreshold )
		{
			trackSeqNew[i].Clear();
		}
		else
		{
			trackSeqNew[i].SetGround( newTracks[i-nOldTracks].GetGround() );
		}
	}

	cameraIDOrder.push_back(newCameraIndex);
	cameraVisited[newCameraIndex] = 1;

	return 0;
}

//#ifdef CERES_LIB
int CeresBA( vector<TrackInfo> trackSeq, vector<ImgFeature> imageFeatures, 
			 /*vector<int> cameraIDOrder,*/	vector<CameraPara> &cameras)
{
	int  num_cameras = cameras.size(); 

	//insert intrinsic camera parameters
	double* pInteriorParams = new double[3];          //focal length, k1, k2
	pInteriorParams[0] = cameras[0].focalLen;
	pInteriorParams[1] = 0;
	pInteriorParams[2] = 0;	
	
	//insert extrinsic camera parameters
	double* pOuterParams = new double[num_cameras*6]; //axis-angle vector, t0,t1,t2
	for(int i=0; i<num_cameras; i++)
	{		
		double aa[3];
		rot2aa(cameras[i].R, aa);
		pOuterParams[i*6]   = aa[0];
		pOuterParams[i*6+1] = aa[1];
		pOuterParams[i*6+2] = aa[2];
		
		/*
		pOuterParams[i*6]   = cameras[i].ax;
		pOuterParams[i*6+1] = cameras[i].ay;
		pOuterParams[i*6+2] = cameras[i].az;
		*/

		double iR[9];
		memcpy(iR, cameras[i].R, sizeof(double)*9);
		invers_matrix(iR, 3);
		double it[3];
		mult(iR, cameras[i].T, it, 3, 3, 1);

		pOuterParams[i*6+3] = -it[0];
		pOuterParams[i*6+4] = -it[1];
		pOuterParams[i*6+5] = -it[2];
		
		//pOuterParams[i*6+3] = cameras[i].t[0];
		//pOuterParams[i*6+4] = cameras[i].t[1];
		//pOuterParams[i*6+5] = cameras[i].t[2];
	}
	
	//find the reasonable tracks 
	vector<int> goodTrackIndex;
	for(int i=0; i<trackSeq.size(); i++)
	{
		//the track is invalid according to the 
		if(trackSeq[i].valid<1)
			continue;

		if(trackSeq[i].derror<4)
			goodTrackIndex.push_back(i);
	}

	int num_pts = goodTrackIndex.size();

	//insert ground point parameters
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

	ceres::Problem problem;
	//invoke ceres functions, each time adding one projection
	for(int i=0; i<nProjection; i++)
	{
		ceres::CostFunction* cost_function =
			SFMReprojectionError::Create(projections[2 * i + 0], projections[2 * i + 1]);

		//printf("%lf %lf \n", projections[2 * i + 0], projections[2 * i + 1] );

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

	cout<<"NumParameterBlocks: "<<problem.NumParameterBlocks()<<"\n";
	cout<<"NumResidualBlocks: "<<problem.NumResidualBlocks()<<"\n";
	cout<<"NumParameters: "<<problem.NumParameters()<<"\n";

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	//std::cout << summary.FullReport() << "\n";

	std::cout<<"\n"<<" all cameras optimized results: "<<"\n";

	printf("Intrisic parameters... \n");
	for(int i=0; i<3; i++)
		std::cout<<pInteriorParams[i]<<"  ";
	std::cout<<endl;
	
	printf("Extrisic parameters ...\n");
	for(int i=0; i<num_cameras; i++)
	{
		/*
		for(int j=0; j<6; j++)
		{
			std::cout<<pOuterParams[i*6+j]<<"   ";
		}
		*/
		
		cameras[i].focalLen = pInteriorParams[0];
		cameras[i].k1    = pInteriorParams[1];
		cameras[i].k2    = pInteriorParams[2];
		
		//from axis-angle to rotation matrix
		double aa[3];
		aa[0] = pOuterParams[i*6];
		aa[1] = pOuterParams[i*6+1];
		aa[2] = pOuterParams[i*6+2];
		double Rt[9];
		aa2rot(aa, Rt);
		double R[9];
		transpose(Rt, R, 3, 3);
		memcpy(cameras[i].R, R, sizeof(double)*9);
		double ea[3];
		rot2eular(R, ea);
		cameras[i].ax = ea[0];
		cameras[i].ay = ea[1];
		cameras[i].az = ea[2];
		
		cout<<"Angle: "<< cameras[i].ax <<" "<< cameras[i].ay<<" " << cameras[i].az << "\n"; 

		/*
		double ax,ay,az;
		ax = pOuterParams[i*6];
		ay = pOuterParams[i*6+1];
		az = pOuterParams[i*6+2];
		GenerateRMatrixDirect(ax,ay,az,cameras[i].R);
		double ea[3];
		rot2eular(cameras[i].R, ea);
		cameras[i].ax = ea[0];
		cameras[i].ay = ea[1];
		cameras[i].ay = ea[2];
		*/

		//from RX+T to RX-T
		double iR[9];
		memcpy(iR, R, sizeof(double)*9);
		invers_matrix(iR, 3);
		double it[3];
		it[0] = pOuterParams[i*6+3];
		it[1] = pOuterParams[i*6+4];
		it[2] = pOuterParams[i*6+5];
		double et[3];
		mult(iR, it, et, 3, 3, 1);
		cameras[i].T[0] = -et[0];
		cameras[i].T[1] = -et[1];
		cameras[i].T[2] = -et[2];

		cout<<"Translation: "<<cameras[i].T[0]<<" "<<cameras[i].T[1]<<" "<<cameras[i].T[2]<<"\n";

		std::cout<<endl;
	}

	for(int i=0; i<num_pts; i++)
	{   
		int index = goodTrackIndex[i];

		trackSeq[index].grd.p[0] = grdPt[i*3];
		trackSeq[index].grd.p[1] = grdPt[i*3+1];
		trackSeq[index].grd.p[2] = grdPt[i*3+2];
	}

	return num_pts;
}


int CeresBA( vector<Point3DDouble>& grdPts, vector<Point2DDouble>& imgPts, CameraPara& camera, 
	bool bIsFixedFocalLen )
{
	double  pInteriorParams[3];
	double  pOuterParams[6];

	int nProjection = imgPts.size();

	pInteriorParams[0] = camera.focalLen;
	pInteriorParams[1] = camera.k1;
	pInteriorParams[2] = camera.k2;

	//calculate the axis-angle vector
	double aa[3];
	rot2aa(camera.R, aa);
	pOuterParams[0] = aa[0];
	pOuterParams[1] = aa[1];
	pOuterParams[2] = aa[2];

	/*
	double iR[9];
	memcpy(iR, camera.R, sizeof(double)*9);
	invers_matrix(iR, 3);
	double it[3];
	mult(iR, camera.t, it, 3, 3, 1);
	*/

	//from R(X-T) to RX+T, revised by xiedonghai, 2016.12.30
	double it[3];
	mult(camera.R, camera.T, it, 3, 3, 1);

	pOuterParams[3] = -it[0];
	pOuterParams[4] = -it[1];
	pOuterParams[5] = -it[2];

	printf("\n 3D Points: %d \n", grdPts.size());

	ceres::Problem problem;
	//invoke ceres functions, each time adding one projection
	for(int i=0; i<nProjection; i++)
	{
		if( camera.camtype == PerspectiveCam )
		{
			if(bIsFixedFocalLen)
			{
				ceres::CostFunction* cost_function =
					RefineCameraFixedFocalLen::Create( imgPts[i].p[0] , imgPts[i].p[1], 
					camera.focalLen, grdPts[i].p[0], grdPts[i].p[1], grdPts[i].p[2]);

				problem.AddResidualBlock(cost_function,
					NULL /* squared loss */,
					pOuterParams       //camera extrinsic parameters: 6
					);
			}
			else
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
		}
		else if( camera.camtype == PanoramCam )
		{
			double radius = (double)(camera.cols) / (2*PI);

			ceres::CostFunction* cost_function =
				RefinePanoCameraError::Create( imgPts[i].p[0] , imgPts[i].p[1], 
				grdPts[i].p[0], grdPts[i].p[1], grdPts[i].p[2], radius);

			problem.AddResidualBlock(cost_function,
				NULL /* squared loss */,
				pOuterParams       //camera extrinsic parameters: 6
				);
		}
	}

	//save the optimized results
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_SCHUR;
	options.minimizer_progress_to_stdout = true;

	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);
	//std::cout << summary.FullReport() << "\n";

	std::cout<<"optimized results: "<<"\n";
	printf("Inner parameters... \n");
	for(int i=0; i<3; i++)
		std::cout<<pInteriorParams[i]<<"  ";
	std::cout<<endl;

	printf("Extrinsic parameters ...\n");

	for(int i=0; i<6; i++)
	{
		std::cout<<pOuterParams[i]<<"   ";
	}
	std::cout<<"\n";
	
	//save the optimized results
	camera.focalLen = pInteriorParams[0];
	camera.k1    = pInteriorParams[1];
	camera.k2    = pInteriorParams[2];
	
	//aa[3];
	aa[0] = pOuterParams[0];
	aa[1] = pOuterParams[1];
	aa[2] = pOuterParams[2];
	double Rt[9];
	aa2rot(aa, Rt);
	double R[9];
	transpose(Rt, R, 3, 3);
	memcpy(camera.R, R, sizeof(double)*9);
	//calculate the Eular angle
	double ea[3];
	rot2eular(R, ea);
	camera.ax = ea[0];
	camera.ay = ea[1];
	camera.az = ea[2];
	

	//from RX+T to R(X-T)
	double iR[9];
	memcpy(iR, R, sizeof(double)*9);
	invers_matrix(iR, 3);
	//double it[3];
	it[0] = pOuterParams[3];
	it[1] = pOuterParams[4];
	it[2] = pOuterParams[5];
	double et[3];
	mult(iR, it, et, 3, 3, 1);
	camera.T[0] = -et[0];
	camera.T[1] = -et[1];
	camera.T[2] = -et[2];

	return 0;
}

//#endif

////////////////////////////////////////////////////////////////////////////

/* calculate the projection error of each image point included in current track sequence,
   and judge if the point is a wrong point.if the point is a bad point, then erase it
   from the corresponding track. Do this one image after image, the threshold is from 
   the statistics of projection errors distributed in one image.
*/
int FindProjectionOutliersByImagePt(vector<TrackInfo>& trackSeq, vector<ImgFeature>& imageFeatures, 
	vector<int> cameraIDOrder, vector<TrackInfo>& tracks, vector<CameraPara>& cameras)
{
	for(int i=0; i<cameraIDOrder.size(); i++)
	{
		int imageId  = cameraIDOrder[i];
		int nFeatNum = imageFeatures[imageId].GetFeatPtSum(); //imageFeatures[imageId].featPts.size();

		for(int j=0; j<nFeatNum; j++)
		{
			int trackId    = imageFeatures[imageId].GetPtTrackIndex(j); //imageFeatures[imageId].featPts[j].extra;
			if(trackId<0)
				continue;

			int curTrackId = tracks[trackId].GetRemapTrackId(); 
			if(curTrackId<0)
				continue;

			int nViews = trackSeq[curTrackId].GetImageKeySum();
			for(int k=0; k<nViews; k++)
			{
				ImageKey ik = trackSeq[curTrackId].GetImageKey(k);
				int imageId = ik.first;
				int ptId = ik.second;

				Point3DDouble gp = trackSeq[curTrackId].GetGround();
				Point2DDouble ip = imageFeatures[imageId].GetCenteredPt(ptId);

				Point2DDouble tp;
				GrdToImg(gp, tp, cameras[imageId]);

				double dx = ip.p[0] - tp.p[0];
				double dy = ip.p[1] - ip.p[1];
				double err = sqrt(dx*dx+dy*dy);
			}
		}
	}
	return 0;
}


/*  calculate the total projection error of each track, if the error is large, 
    then discard the track directly.   
	return value: the number of bad tracks
*/
int FindProjectionOutliersByTrack(vector<TrackInfo>& trackSeq, vector<ImgFeature>& imageFeatures, 
	vector<int> cameraIDOrder, vector<TrackInfo>& tracks, vector<CameraPara>& cameras, 
	double threshold)
{	
	vector<TrackInfo> newTracks;
	int nBadTracks = 0;
	for(int i=0; i<trackSeq.size(); i++)
	{
		int nview        = trackSeq[i].GetImageKeySum();
		Point3DDouble gp = trackSeq[i].GetGround();

		double totalError = 0;
		vector<double> errors;
		int nInliers = 0;
		for(int k=0; k<nview; k++)
		{
			ImageKey ik = trackSeq[i].GetImageKey(k);
			int imageId = ik.first;
			int ptId    = ik.second;

			Point2DDouble ip = imageFeatures[imageId].GetCenteredPt(ptId);

			Point2DDouble projPt;
			GrdToImg(gp, projPt, cameras[imageId]);

			double dx = ip.p[0] - projPt.p[0];
			double dy = ip.p[1] - projPt.p[1];
			double err = sqrt(dx*dx+dy*dy);
			totalError += err;
			errors.push_back(err);

			if(err<threshold)
				nInliers++;
		}
		
		if(nInliers>=2)
			newTracks.push_back(trackSeq[i]);
		else
			nBadTracks ++;
	}
	
	trackSeq = newTracks;
	
	return nBadTracks;
}

int RefineAllParametersWithPos(vector<TrackInfo>& tracks, vector<ImgFeature>& imageFeatures,
	vector<int> cameraIDOrder, vector<CameraPara>& cameras, double az)
{

	int nOutliers = 0;
	int nTotalOutliers = 0;
	do
	{
		//1. optimization 
		RunBA(tracks, imageFeatures, cameraIDOrder, cameras);

		//2. find the wrong feature point projections and remove them
		nOutliers = RemoveOutlierPtsWithPos(tracks, imageFeatures, cameraIDOrder, cameras, az);
		nTotalOutliers += nOutliers;

	} while (nOutliers>0);

	//remove the bad points according to the angle
	//int num_pruned = RemoveBadPoints(trackSeq, imageFeatures, cameras);

	return nTotalOutliers;


	return 0;
}

//refine the camera parameters and track points 
int RefineAllParameters(vector<TrackInfo>& trackSeq, vector<ImgFeature>& imageFeatures, 
	vector<int> cameraIDOrder, vector<TrackInfo>& tracks, vector<CameraPara>& cameras)
{
	int nOutliers = 0;
	int nTotalOutliers = 0;
	int nIteration = 0;
	do
	{ 
		//1. optimization 
		RunBA(trackSeq, imageFeatures, cameraIDOrder, cameras);
		
		//2. find the wrong feature point projections and remove them
		nOutliers = RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);
		nTotalOutliers += nOutliers;

        //3. remove the bad images


		nIteration++;

	}while(nOutliers>0 && nIteration<5);
	
	//remove the bad points according to the angle
	//int num_pruned = RemoveBadPoints(trackSeq, imageFeatures, cameras);

	return nTotalOutliers;
}

int FindGoodPoints(vector<Point3DDouble>& grdPts, vector<Point2DDouble>& imgPts, CameraPara& camera,
					double projThreshold)
{
	int nInliers = 0;
	
	//calculate all projection errors and get the threshold
	double* pErrors = new double[grdPts.size()];
	for(int i=0; i<grdPts.size(); i++)
	{
		Point2DDouble projPt;
		GrdToImg(grdPts[i], projPt, camera);

		double dx = imgPts[i].p[0] - projPt.p[0];
		double dy = imgPts[i].p[1] - projPt.p[1];
		double err = sqrt(dx*dx+dy*dy);
		pErrors[i] = err;
	}
	double med = dll_kth_element_copy(grdPts.size(), dll_iround(0.8 * grdPts.size()), pErrors);
	double threshold = 1.2 * NUM_STDDEV * med; 

	if(camera.camtype == PerspectiveCam)
	{
		threshold = CLAMP(threshold, MIN_PROJ_ERROR_THRESHOLD, MAX_PROJ_ERROR_THRESHOLD);  
	}
	else if(camera.camtype == PanoramCam)
	{
		double radius = (double)(camera.cols) / (2*PI);
		double minT = radius*0.01;
		double maxT = radius*0.1;
		threshold = CLAMP(threshold, minT, maxT);
	}
	
	//find the good points
	vector<Point3DDouble> newGrdPts;
	vector<Point2DDouble> newImgPts;
	for(int i=0; i<grdPts.size(); i++)
	{
		/*
		Point2DDouble projPt;
		GrdToImg(grdPts[i], projPt, camera);
		double dx = imgPts[i].p[0] - projPt.p[0];
		double dy = imgPts[i].p[1] - projPt.p[1];
		double err = sqrt(dx*dx+dy*dy);
		*/

		if(pErrors[i]<threshold)
		{
			newGrdPts.push_back(grdPts[i]);
			newImgPts.push_back(imgPts[i]);
			nInliers++;
		}
	}

	delete[] pErrors;

	grdPts = newGrdPts;
	imgPts = newImgPts;

	return nInliers;
}


//refine only the camera parameters iterately after DLT 
//return value: 
//  0-camera is good and can be added into bundle adjustment
//  1-camera is bad 
int RefineCamera( vector<Point3DDouble>& grdPts, vector<Point2DDouble>& imgPts, CameraPara& camera)
{
	//1. optimization with fixed camera focal length
	CeresBA(grdPts, imgPts, camera, true);

	int nCurrentInliers = FindGoodPoints(grdPts, imgPts, camera, 8);
	int nNewInliers = nCurrentInliers;
	while(1)
	{
		//2. optimization for single camera parameters
		CeresBA(grdPts, imgPts, camera, true);
				
		//3. calculate the errors
		nNewInliers = FindGoodPoints(grdPts, imgPts, camera, 8);

		//4. evaluate the optimization
		if(nCurrentInliers==nNewInliers)
			break;

		nCurrentInliers = nNewInliers;
	}
	
	return nNewInliers;
}

//calculate the camera intrinsic and extrinsic parameters of camera
int CalculateNewCamParas(int nCameraId, 
	vector<ImgFeature>&  imageFeatures,
	vector<TrackInfo>& trackSeqNew, 
	vector<TrackInfo>& tracks,
	double initFocalLen,
	CameraPara& cam )
{
	vector<Point3DDouble> pt3;
	vector<Point2DDouble> pt2;

	//collect the track corresponding to the new camera
	for(int i=0; i<trackSeqNew.size(); i++)
	{
		if( trackSeqNew[i].GetImageKeySum() == 0 )
			continue;

		int nTrackId = trackSeqNew[i].GetRemapTrackId();  //find the original track

		int nviews = tracks[nTrackId].GetImageKeySum();
		for(int j=0; j<nviews; j++)
		{
			ImageKey ik = tracks[nTrackId].GetImageKey(j);
			int nImageId = ik.first;
			int nPtId    = ik.second;

			if(nImageId == nCameraId)
			{
				pt3.push_back( trackSeqNew[i].GetGround() );
				Point2DDouble p2 = imageFeatures[nImageId].GetCenteredPt(nPtId);
				pt2.push_back(p2);
			}
		}
	}

	//ransac DLT calculation
	CPoseEstimationBase* pPE = NULL;
	
	if( cam.camtype == PerspectiveCam )
		pPE = new CDLTPose();
	else if(cam.camtype == PanoramCam)
		pPE = new CPanoDLTPose();

	int r = pPE->EstimatePose(pt3, pt2, cam);

	vector<int> inliers = pPE->GetInliers();
	vector<int> inliers_weak = pPE->GetWeakInliers();
	delete pPE;
	
	//if successful, do bundle adjustment
	if(r==0)
	{
		printf("\n optimization for DLT results: \n");

		//remove the wrong projections
		vector<Point3DDouble> inlierPt3;
		vector<Point2DDouble> inlierPt2;

		for(int i=0; i<inliers_weak.size(); i++)
		{
			inlierPt3.push_back( pt3[inliers_weak[i]] );
			inlierPt2.push_back( pt2[inliers_weak[i]] );
		}

		cam.focalLen = initFocalLen;
		int nInliers = RefineCamera(inlierPt3, inlierPt2, cam);
		
		//
		if (nInliers < 8)
		{
			printf("inliers of DLT are too fewer! ,  return -1!  \n");
			return -1;
		}
	}
	else
	{
		return -1;
	}

	return 0;
}

//select the best pair from all match pairs
int PickInitialPair(vector<PairMatchRes>& pairMatchs)
{
	int pid = -1;

	double minHInliers = 10000; 
	double maxMatch = 0;
    double minScore = 1;
	for (int i = 0; i < pairMatchs.size(); i++)
	{
		/*
		if (pairMatchs[i].GetMatchNumber()>128)
		{
			if (minHInliers > pairMatchs[i].inlierRatioHomography)
			{
				minHInliers = pairMatchs[i].inlierRatioHomography;
				pid = i;
			}
		}*/

		/*
		if (   ( pairMatchs[i].GetMatchNumber() > maxMatch ) 
			&& ( pairMatchs[i].inlierRatioHomography<(pairMatchs[i].inlierRatio) ) )
			{
				maxMatch = pairMatchs[i].GetMatchNumber();
				pid = i;
			}*/

		//if (   pairMatchs[i].GetMatchNumber() > maxMatch
		//	&& pairMatchs[i].GetMatchNumber()>100) 
        if(pairMatchs[i].GetMatchNumber()>MIN_MATCH_NUM && pairMatchs[i].GetScore()<minScore)
        {
            minScore = pairMatchs[i].GetScore();
			//maxMatch = pairMatchs[i].GetMatchNumber();
			pid = i;
		}
	}

	return pid;
}


/////////////////////////////////////////////////////////////////////////////
CCeresBA::CCeresBA()
{
	//mpdProgress = NULL;
}

CCeresBA::~CCeresBA()
{

}

void  CCeresBA::SetProgress(double* pdProgress)
{
	//mpdProgress = pdProgress;
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
	vector<TrackInfo>& tracks, //track points ( all connecting points )
    bool& bIsStop)           
{
	if (numCameras < 2)
		return -1;

	int ht = cameras[0].rows;
	int wd = cameras[0].cols;
	printf("height:%d  width:%d \n", ht, wd);
	
	//1. select the initial pair
	int index = PickInitialPair(pairMatchs);
	printf("initial pair index: %d \n", index);
	if (index < 0){
		printf("failed to pick initial pair ! \n");
		return -1;
	}

	vector<int> cameraVisited;
	cameraVisited.resize(numCameras, 0);

	int leftImageId  = pairMatchs[index].lId;
	int rightImageId = pairMatchs[index].rId;
	vector<TrackInfo> trackSeq;  //current tracks for bundle adjustment
	GetMatch(leftImageId, rightImageId, tracks, trackSeq);    


	cameraVisited[leftImageId]  = 1;
	cameraVisited[rightImageId] = 1;
	cameras[leftImageId].bIsAddedIntoNet = true;
	cameras[rightImageId].bIsAddedIntoNet = true;

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
	//vector<Point3DDouble> goodPts;
	for(int i=0; i<gpts.size(); i++)
	{
		trackSeq[i].grd = gpts[i];
		trackSeq[i].derror = errorarray[i];
	}
	
	RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);

#ifdef _DEBUG
	SaveTracksToPly("c:\\temp\\pair-raw.ply", trackSeq, cameraIDOrder, cameras);
#endif

    SaveTracksToPly("c:\\temp\\pair-init.ply", trackSeq, cameraIDOrder, cameras);

	//CeresBA(trackSeq, imageFeatures, cameras);
	RunBA(trackSeq, imageFeatures, cameraIDOrder, cameras);

	RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);

    SaveTracksToPly("c:\\temp\\pair-ba.ply", trackSeq, cameraIDOrder, cameras);


#ifdef _DEBUG
	SaveTracksToPly("c:\\temp\\pair.ply", trackSeq, cameraIDOrder, cameras);
#endif

	//3. adding new images
	vector<int> lostCameras; //cameras lost in DLT 

	//*mpdProgress = 0;
    mdProgress = 0;
	double step = 100.0 / double(cameras.size());
    
	while(1)
	{
        if (bIsStop) {
            return -1;
        }

		//*mpdProgress += step;
        mdProgress += step;

		int newCameraIndex = SelectNewImage(cameraVisited, tracks, trackSeq);
		printf("\n new camera index: %d \n", newCameraIndex);

		if (newCameraIndex < 0)
		{
			printf("can not find new camera, ba is over!");
			break;
		}
		
		//calculate the camera parameters of new selected
		printf("\n\n************** get initial parameters for camera: %d ******************** \n", newCameraIndex);
		int initFocalLen = cameras[cameraIDOrder[0]].focalLen; //use the optimized focal length as initial value
		int res = CalculateNewCamParas(newCameraIndex, imageFeatures, trackSeq, tracks, 
			initFocalLen, cameras[newCameraIndex]);
		
		//DLT pose estimation failed
		if( res < 0 )
		{
			cameraVisited[newCameraIndex] = 1;
			cameras[newCameraIndex].focalLen = 0; //added by xiedonghai, 2018.5.19
			cameras[newCameraIndex].bIsAddedIntoNet = false;
			lostCameras.push_back(newCameraIndex);
			printf("\n discard camera: %d \n", newCameraIndex);
			continue;
		}

		cameras[newCameraIndex].bIsAddedIntoNet = true;

		
		//adding the image points of new image into the current tracks 
		printf("\n\n\n ******************* adding image %d ******************** \n", newCameraIndex);
		UpdateBATracks(newCameraIndex, cameraVisited, cameraIDOrder, imageFeatures,
			tracks, trackSeq, cameras);


		//local BA, added by xiedonghai, 2018.4.23
		//printf("local BA.... \n");
		//vector<TrackInfo> localTracks = FindLocalTracks(newCameraIndex, trackSeq, tracks, imageFeatures);
		//vector<int> localCameraIds = GetCameraInTracks(localTracks);
		//printf("local cameras.....\n");
		//for (int i = 0; i < localCameraIds.size(); i++)
		//	printf("% ", localCameraIds[i]);
		//printf("\n");
		//RunBA(localTracks, imageFeatures, localCameraIds, cameras);
		//for (int i = 0; i < localTracks.size(); i++)
		//{
		//	int ti = localTracks[i].extra;

		//	//revised by xdh, 2018.4.26, only update the 3D coordinates 
		//	trackSeq[ti].SetGround( localTracks[i].GetGround()); 
		//}

		//if (cameraIDOrder.size() % 20 == 0)
		//{
		//	RefineAllParameters(trackSeq, imageFeatures, cameraIDOrder, tracks, cameras);
		//}

		////RemoveOutlierPts(tracks, trackSeq, imageFeatures, localCameraIds, cameras);
		//char tfile[256];
		//sprintf(tfile, "c:\\temp\\final-%d.ply", newCameraIndex);
		//SaveTracksToPly(tfile, trackSeq, cameraIDOrder, cameras);

		//for debug
		//RunBA(trackSeq, imageFeatures, cameraIDOrder, cameras);
		//RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);
		//char tfile[256];
		//sprintf(tfile, "c:\\temp\\final-%d.ply", newCameraIndex);
		//SaveTracksToPly(tfile, trackSeq, cameraIDOrder, cameras);

		//optimization for all cameras
		//printf("\n BA for all cameras... \n");
		RefineAllParameters(trackSeq, imageFeatures, cameraIDOrder, tracks, cameras);

        //char plyfile[256];
        //sprintf(plyfile, "c:\\temp\\%d.ply", newCameraIndex);
        //SaveTracksToPly(plyfile, trackSeq, cameraIDOrder, cameras);
        printf("camera parameters .... \n");
        for (int k = 0; k < cameraIDOrder.size(); k++)
        {
            int id = cameraIDOrder[k];
            printf("image:%d   f:%lf   Ts: %lf %lf %lf \n", 
                id, cameras[id].focalLen, 
                cameras[id].xs, cameras[id].ys, cameras[id].zs);
        }

        mpCallFunc(NULL, mdProgress, "Bundle Adjustment...");
	}
	
	//printf("\n Final BA for all cameras... \n");
	//RefineAllParameters(trackSeq, imageFeatures, cameraIDOrder, tracks, cameras);

	/*
	/////////////   reprocess the images lost in DLT, added by xdh, 2018.4.26 //////////////////
	vector<int> successfulCameras;
	for (int i = 0; i < cameraVisited.size(); i++)
	{
		bool bisfound = false;
		for (int k = 0; k < lostCameras.size(); k++)
		{
			if (cameraVisited[i] == lostCameras[k])
			{
				bisfound = true;
			}
		}
		if (!bisfound)
			successfulCameras.push_back(cameraVisited[i]);
	}
	cameraVisited = successfulCameras;

	for (int i = 0; i < lostCameras.size(); i++)
	{
		int cid = lostCameras[i];
		
		int initFocalLen = cameras[cameraIDOrder[0]].focalLen;
		int res = CalculateNewCamParas(cid, imageFeatures, trackSeq, tracks,
			initFocalLen, cameras[cid]);

		if (res>0)
		{
			UpdateBATracks(cid, cameraVisited, cameraIDOrder, imageFeatures, tracks, trackSeq, cameras);
			cameraIDOrder.push_back(cid);
			cameraVisited.push_back(cid);
		}
	}
	RefineAllParameters(trackSeq, imageFeatures, cameraIDOrder, tracks, cameras);
	////////////////////////////////////////////////////////////////////////////////////////////////////
	*/


	//SaveTracksToPly("c:\\temp\\final-ba.ply", trackSeq, cameraIDOrder, cameras);

	//output the camera parameters
	printf("\n ****************** Final results: ********************** \n");
	printf("Camera Number: %d \n", cameraIDOrder.size());

	for (int i = 0; i < cameras.size(); i++)
	{
		cameras[i].bIsAddedIntoNet = false;
	}
	for(int i=0; i<cameraIDOrder.size(); i++)
	{
		int camId = cameraIDOrder[i];
		cameras[camId].bIsAddedIntoNet = true;
		printf("camera: %d ", camId);
		printf(" focal: %lf  ", cameras[camId].focalLen);
		printf(" angle: %lf %lf %lf  ", cameras[camId].ax, cameras[camId].ay, cameras[camId].az);
		printf(" position: %lf %lf %lf \n", cameras[camId].T[0], cameras[camId].T[1], cameras[camId].T[2]);
	}
	
	//save the good tracks, added by xdh, 2018.1.9
	//tracks.clear();
	//tracks = trackSeq;
	GetGoodTracks(trackSeq, tracks);

    SaveTracksToPly("c:\\temp\\final-ba.ply", trackSeq, cameraIDOrder, cameras);


	return 0;
}