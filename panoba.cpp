

/*
   



*/


#include"defines.hpp"
#include"panoba.hpp"
#include"relativepose.hpp"
#include"bundlerio.hpp"


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
	CTriangulateBase* pTri = new CTriangulatePano();
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
	
	//RemoveOutlierPts(tracks, trackSeq, imageFeatures, cameraIDOrder, cameras);
  
	SaveTracksToPly("c:\\temp\\pair-raw.ply", trackSeq, cameraIDOrder, cameras);


	return 0;
}