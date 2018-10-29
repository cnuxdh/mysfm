

#include"videoMosaic.hpp"


//corelib
#include"commondata.h"
#include"CommonFuncs.h"

//cvlib
#include"badata.hpp"
#include"orbfeat.hpp"
#include"register.hpp"
#include"myui.hpp"


CPlaneMosaicPara::CPlaneMosaicPara()
{
	ma = 1;
	mb = 0;
	mc = 0;
	md = 0;
}

CPlaneMosaicPara::~CPlaneMosaicPara()
{

}

int CPlaneMosaicPara::SetSimilarityParas(double a, double b, double c, double d)
{
	ma = a;
	mb = b;
	mc = c;
	md = d;

	return 0;
}

int CPlaneMosaicPara::GetSimilarityParas(double& a, double& b, double& c, double& d)
{
	a = ma;
	b = mb;
	c = mc;
	d = md;
	return 0;
}


int SimilarityRANSAC(MyPointF* lImagePts, MyPointF* rImagePts, int nPt,
	double& a, double& b, double& c, double& d, vector<int>& vInliers)
{
	int ransac_rounds = 200;
	int maxInliers = 0;
	int DistanceThrshold = 16;
	
	for (int round = 0; round < ransac_rounds; round++)
	{
		int indices[4];
		ChooseRandIndex(nPt, 4, indices);
		
		MyPointF lps[4];
		MyPointF rps[4];
		for (int k = 0; k < 4; k++)
		{
			int randi = indices[k];
			lps[k].x = lImagePts[randi].x;
			lps[k].y = lImagePts[randi].y;
			rps[k].x = rImagePts[randi].x;
			rps[k].y = rImagePts[randi].y;
		}

		double p[4];
		SimilarityTransform(rps, lps, 4, p);

		//error
		double error = 0;
		int inliers = 0;
		vInliers.clear();
		for (int i = 0; i < nPt; i++)
		{
			double tx, ty;
			double dx = lImagePts[i].x;
			double dy = lImagePts[i].y;
			double ox = rImagePts[i].x;
			double oy = rImagePts[i].y;
			tx = SIMILARITY_X(ox, oy, p);
			ty = SIMILARITY_Y(ox, oy, p);

			double len = sqrt((tx - dx)*(tx - dx) + (ty - dy)*(ty - dy));
			if (len < DistanceThrshold)
			{
				vInliers.push_back(i);
				inliers++;
			}
			error += len;
		}
		error /= double(nPt);

		if (maxInliers < inliers)
		{
			maxInliers = inliers;
			a = p[0];
			b = p[1];
			c = p[2];
			d = p[3];
		}
	}

	return maxInliers;
}


//initialize the ground points of pair
int VideoMosaicInitialize(vector<TrackInfo>&  tracks, vector<CVideoKeyFrame>& keyFrames)
{
	for (int i = 0; i < tracks.size(); i++)
	{
		int nview = tracks[i].GetViews();
		if (nview == 0)
			continue;

		double gx = 0;
		double gy = 0;
		for (int k = 0; k < nview; k++)
		{
			ImageKey ik = tracks[i].GetImageKey(k);
			int imageIndex = ik.first;
			int ptIndex = ik.second;

			double cx = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cx;
			double cy = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cy;

			double p[4];
			keyFrames[imageIndex].mCamera.GetSimilarityParas(p[0], p[1], p[2], p[3]);

			gx += SIMILARITY_X(cx, cy, p);
			gy += SIMILARITY_Y(cx, cy, p);
		}
		gx /= double(nview);
		gy /= double(nview);

		//calculate the error
		double error = 0;
		for (int k = 0; k < nview; k++)
		{
			ImageKey ik = tracks[i].GetImageKey(k);
			int imageIndex = ik.first;
			int ptIndex = ik.second;

			double cx = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cx;
			double cy = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cy;

			double p[4];
			keyFrames[imageIndex].mCamera.GetSimilarityParas(p[0], p[1], p[2], p[3]);

			double tx = SIMILARITY_X(cx, cy, p);
			double ty = SIMILARITY_Y(cx, cy, p);

			error += sqrt( (gx - tx)*(gx - tx) + (gy - ty)*(gy - ty) );
		}
		error /= double(nview);

		if (error<16)
		{
			Point3DDouble gp;
			gp.p[0] = gx;
			gp.p[1] = gy;
			gp.p[2] = 0;
			tracks[i].SetGround(gp);
		}
		else
		{
			//outliers
			tracks[i].Clear();
		}
	}	
	return 0;
}


int VisibleTracksInitialize(vector<CVideoKeyFrame>& keyFrames,
	vector<TrackInfo>&  tracks, vector<int> visibleIndex)
{
	for (int i = 0; i < visibleIndex.size(); i++)
	{
		int nTrack = visibleIndex[i];

		int nview = tracks[nTrack].GetViews();
		if (nview == 0)
			continue;

		double gx = 0;
		double gy = 0;
		for (int k = 0; k < nview; k++)
		{
			ImageKey ik = tracks[nTrack].GetImageKey(k);
			int imageIndex = ik.first;
			int ptIndex = ik.second;

			double cx = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cx;
			double cy = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cy;

			double p[4];
			keyFrames[imageIndex].mCamera.GetSimilarityParas(p[0], p[1], p[2], p[3]);

			gx += SIMILARITY_X(cx, cy, p);
			gy += SIMILARITY_Y(cx, cy, p);
		}
		gx /= double(nview);
		gy /= double(nview);

		//calculate the error
		double error = 0;
		for (int k = 0; k < nview; k++)
		{
			ImageKey ik = tracks[nTrack].GetImageKey(k);
			int imageIndex = ik.first;
			int ptIndex = ik.second;

			double cx = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cx;
			double cy = keyFrames[imageIndex].mFeature.GetFeatPt(ptIndex).cy;

			double p[4];
			keyFrames[imageIndex].mCamera.GetSimilarityParas(p[0], p[1], p[2], p[3]);

			double tx = SIMILARITY_X(cx, cy, p);
			double ty = SIMILARITY_Y(cx, cy, p);

			error += sqrt((gx - tx)*(gx - tx) + (gy - ty)*(gy - ty));
		}
		error /= double(nview);

		if (error<16)
		{
			Point3DDouble gp;
			gp.p[0] = gx;
			gp.p[1] = gy;
			gp.p[2] = 0;
			tracks[nTrack].SetGround(gp);
		}
		else
		{
			//outliers
			tracks[nTrack].Clear();
		}
	}

	return 0;
}

//get the tracks seen by the current frame
int GetCurrentFrameVisibleTracks(vector<TrackInfo>& tracks, 
	vector<CVideoKeyFrame>& keyFrames,
	vector<PairMatchRes> pm, 
	int currentFrameIndex,
	vector<int>& visibleTrackIndex,
	vector<int>& currentFramePtIndex)
{
	visibleTrackIndex.clear();
	currentFramePtIndex.clear();

	for (int i = 0; i < pm.size(); i++)
	{   
		int lImageIndex = pm[i].lId;
		int rImageIndex = pm[i].rId;

		//tranverse the match pair
		for (int k = 0; k < pm[i].matchs.size(); k++)
		{
			int l = pm[i].matchs[k].l;
			int r = pm[i].matchs[k].r;
		    
			int nTrackIndex;
			if (lImageIndex == currentFrameIndex)
			{
				nTrackIndex = keyFrames[rImageIndex].mFeature.GetFeatPt(r).extra;
				if (nTrackIndex > 0)
				{
					if (tracks[nTrackIndex].GetViews() == 0)
						continue;

					currentFramePtIndex.push_back(l);
					visibleTrackIndex.push_back(nTrackIndex);
				}
			}
			else
			{
				nTrackIndex = keyFrames[lImageIndex].mFeature.GetFeatPt(l).extra;
				if (nTrackIndex > 0)
				{
					if (tracks[nTrackIndex].GetViews() == 0)
						continue;

					currentFramePtIndex.push_back(r);
					visibleTrackIndex.push_back(nTrackIndex);
				}
			}
		}
	}

	/*
	for (int i = 0; i < tracks.size(); i++)
	{
		int nview = tracks[i].GetViews();

		for (int j = 0; j < nview; j++)
		{
			ImageKey ik = tracks[i].GetImageKey(j);
			int imageIndex = ik.first;
			int ptIndex = ik.second;
			if (imageIndex == currentFrameIndex)
			{
				visibleTrackIndex.push_back(i);
				currentFramePtIndex.push_back(ptIndex);
				break;
			}
		}
	}*/

	
	return 0;
}


//get the tracks seen by the current frame
int GetCurrentVisibleTracksNew(vector<TrackInfo>& tracks,
	ImgFeature& currentFrame, vector<int>& visibleTrackIndex)
{
	visibleTrackIndex.clear();
	
	for (int i = 0; i < currentFrame.GetFeatPtSum(); i++)
	{
		int trackIndex = currentFrame.GetFeatPt(i).extra;

		if (trackIndex < 0)
			continue;

		if (tracks[trackIndex].GetViews() == 0)
			continue;
				
		visibleTrackIndex.push_back(trackIndex);
	}


	/*for (int i = 0; i < tracks.size(); i++)
	{
		int nview = tracks[i].GetViews();

		for (int j = 0; j < nview; j++)
		{
			ImageKey ik = tracks[i].GetImageKey(j);
			
			int imageIndex = ik.first;
			int ptIndex = ik.second;

			if (imageIndex == currentFrameIndex)
			{
				visibleTrackIndex.push_back(i);
				break;
			}
		}
	}*/

	return 0;
}

//////////////////////////////////////////////////////////
// video local bundle adjustment
CVideoFrameBA::CVideoFrameBA()
{

}
CVideoFrameBA::~CVideoFrameBA()
{

}

int CVideoFrameBA::Run()
{
	


	return 0;
}



CVideoMosaicSystem::CVideoMosaicSystem()
{
	mbIsHavingKeyFrame = false;
	mbIsInitialized = false;
}
CVideoMosaicSystem::~CVideoMosaicSystem()
{

}

int CVideoMosaicSystem::MatchingWithClosetFrame()
{
	//detect feature
	CORBFeat orbFeat;
	mCurrentFrame.Clear();
	orbFeat.Detect(mCurrentFrameImage, mCurrentFrame.mFeature, 640);
	printf("feature number: %d \n", mCurrentFrame.mFeature.GetFeatPtSum());

	//registration with the last keyframe
	int nsize = mKeyFrames.size();
	mNewFrameMatch.Clear();
	mNewFrameMatch.lId = nsize - 1; //index of last image
	mNewFrameMatch.rId = nsize;   //index of new image
	CORBPerspectiveMatch orbMatch;
	orbMatch.Match(mKeyFrames[nsize - 1].GetFeat(), mCurrentFrame.mFeature, mNewFrameMatch);
	printf("match -  lIndex:%d  rIndex:%d  num:%d \n",
		mNewFrameMatch.lId, mNewFrameMatch.rId,
		mNewFrameMatch.matchs.size());

	return 0;
}

int CVideoMosaicSystem::SimTransWithLastKeyFrame(
	double& sa, double& sb, double& sc, double& sd, double& inlierRatio)
{
	int nsize = mKeyFrames.size();

	if (nsize < 1)
		return 0;

	int nMatch = mNewFrameMatch.GetMatchNumber();

	//similarity transform
	MyPointF* srcPt = new MyPointF[nMatch];
	MyPointF* dstPt = new MyPointF[nMatch];
	
	for (int i = 0; i < nMatch; i++)
	{
		int lPtIndex = mNewFrameMatch.GetLeftMatchPtIndex(i);
		int rPtIndex = mNewFrameMatch.GetRightMatchPtIndex(i);

		double x1 = mKeyFrames[nsize - 1].mFeature.GetFeatPt(lPtIndex).cx;
		double y1 = mKeyFrames[nsize - 1].mFeature.GetFeatPt(lPtIndex).cy;
		double x2 = mCurrentFrame.mFeature.GetFeatPt(rPtIndex).cx;
		double y2 = mCurrentFrame.mFeature.GetFeatPt(rPtIndex).cy;

		srcPt[i].x = x2;
		srcPt[i].y = y2;
		dstPt[i].x = x1;
		dstPt[i].y = y1;
	}

	//double sa, sb, sc, sd;
	vector<int> vInliers;
	int inliers = SimilarityRANSAC(dstPt, srcPt, nMatch, sa, sb, sc, sd,
		vInliers);
	delete srcPt;
	delete dstPt;

	inlierRatio = double(inliers) / double(nMatch);
	//printf("right ratio: %lf \n", rightRatio);

	return 0;
}


bool CVideoMosaicSystem::Initialize()
{
	int nsize = mKeyFrames.size();

	if (nsize == 0)
	{
		CORBFeat orbFeat;
		//mCurrentFrame.Clear();
		CVideoKeyFrame vkf;
		orbFeat.Detect(mCurrentFrameImage, vkf.mFeature, 640);
		printf("feature number: %d \n", vkf.mFeature.GetFeatPtSum());
		vkf.mCamera.SetSimilarityParas(1, 0, 0, 0);

		if (vkf.mFeature.GetFeatPtSum() > 1024)
		{
			mKeyFrames.push_back(vkf);
		}

		return false;
	}
	else
	{
		MatchingWithClosetFrame();

		//calculate similarity
		double sa, sb, sc, sd, ratio;
		SimTransWithLastKeyFrame(sa, sb, sc, sd, ratio);
		
		if (ratio < 0.7)
			return false;

		mCurrentFrame.mCamera.SetSimilarityParas(sa, sb, sc, sd);
		mKeyFrames.push_back(mCurrentFrame);
		
		//generate tracks
		CFastGenerateTrack generateTracks;
		vector<PairMatchRes> pm;
		pm.push_back(mNewFrameMatch);
		generateTracks.GenerateTracks(mKeyFrames, pm, mTracks);
		
		//initialize the coordinates of tracks
		VideoMosaicInitialize(mTracks, mKeyFrames);
		
		return true;
	}
	return false;
}

bool CVideoMosaicSystem::IsKeyFrame()
{
	int nsize = mKeyFrames.size();
	
	MatchingWithClosetFrame();

	//if (mNewFrameMatch.GetMatchNumber() < 24)
	//	return false;

	//calculate the translation
	vector<double> translateX;
	vector<double> translateY;
	for (int i = 0; i < mNewFrameMatch.GetMatchNumber(); i++)
	{
		int lPtIndex = mNewFrameMatch.GetLeftMatchPtIndex(i);
		int rPtIndex = mNewFrameMatch.GetRightMatchPtIndex(i);

		double x1 = mKeyFrames[nsize - 1].mFeature.GetFeatPt(lPtIndex).x;
		double y1 = mKeyFrames[nsize - 1].mFeature.GetFeatPt(lPtIndex).y;
		double x2 = mCurrentFrame.mFeature.GetFeatPt(rPtIndex).x;
		double y2 = mCurrentFrame.mFeature.GetFeatPt(rPtIndex).y;

		translateX.push_back(fabs(x1 - x2));
		translateY.push_back(fabs(y1 - y2));
	}

	int half = translateX.size() / 2;

	int ht = mCurrentFrame.mFeature.ht;
	int wd = mCurrentFrame.mFeature.wd;
	double threshold = (wd + ht)*0.5*0.005;

	if (translateX[half] < threshold && translateY[half] < threshold)
		return false;

	return true;
}

int CVideoMosaicSystem::CalculateNewCamPara()
{
	//retrieve the tracks seen by current frame befor inserting new image
	int currentIndex = mKeyFrames.size();
	vector<int> visibleTrackIndex;
	vector<int> currentFramePtIndex;

	vector<PairMatchRes> pm;
	pm.push_back(mNewFrameMatch);
	GetCurrentFrameVisibleTracks(mTracks, mKeyFrames, pm,
		currentIndex, visibleTrackIndex, currentFramePtIndex);
	
	//if (visibleTrackIndex.size() < 12)
	//	return 0;

	//CPlaneMosaicPara planeMosaicP;
	double sa=1, sb=0, sc=0, sd=0;
	int nsize = mKeyFrames.size();
	mKeyFrames[nsize - 1].mCamera.GetSimilarityParas(sa, sb, sc, sd);//copy from the last
	
	if (visibleTrackIndex.size() > 18)
	{
		//initialize the parameters of current frame
		MyPointF* groundPts = new MyPointF[visibleTrackIndex.size()];
		MyPointF* framePts = new MyPointF[visibleTrackIndex.size()];
		for (int i = 0; i < visibleTrackIndex.size(); i++)
		{
			int nTrackIndex = visibleTrackIndex[i];
			groundPts[i].x = mTracks[nTrackIndex].GetGround().p[0];
			groundPts[i].y = mTracks[nTrackIndex].GetGround().p[1];

			int nPtIndex = currentFramePtIndex[i];
			framePts[i].x = mCurrentFrame.mFeature.GetFeatPt(nPtIndex).cx;
			framePts[i].y = mCurrentFrame.mFeature.GetFeatPt(nPtIndex).cy;
		}
		vector<int> vInliers;
		int inliers = SimilarityRANSAC(groundPts, framePts, visibleTrackIndex.size(),
			sa, sb, sc, sd, vInliers);
		printf("similarity: %lf %lf %lf %lf \n", sa, sb, sc, sd);
		delete groundPts;
		delete framePts;
		double rightRatio = double(inliers) / double(visibleTrackIndex.size());
		//if (rightRatio < 0.7)
		//	continue;
	}

	mCurrentFrame.mCamera.SetSimilarityParas(sa, sb, sc, sd);

	//addding the current frame into tracks
	mKeyFrames.push_back(mCurrentFrame);
	CFastGenerateTrack generateTracks;
	pm.push_back(mNewFrameMatch);
	generateTracks.GenerateTracks(mKeyFrames, pm, mTracks);

	//calculate the ground point of current visible tracks
	vector<int> currentTrackIndexNew;
	int nCurrentIndex = mKeyFrames.size() - 1;
	GetCurrentVisibleTracksNew(mTracks, mKeyFrames[nCurrentIndex].mFeature, currentTrackIndexNew);
	VisibleTracksInitialize(mKeyFrames, mTracks, currentTrackIndexNew);
	
	mbIsHavingKeyFrame = true;

	return 0;
}

int CVideoMosaicSystem::WarpImage(Mat& mosaic)
{
	if (mbIsHavingKeyFrame == false)
		return 0;

	////////////////  warp image and save  //////////////////////
	int mht = mosaic.rows;
	int mwd = mosaic.cols;

	int minx, miny, maxx, maxy;
	int ht = mCurrentFrameImage.rows;
	int wd = mCurrentFrameImage.cols;

	double sa, sb, sc, sd;
	mCurrentFrame.mCamera.GetSimilarityParas(sa, sb, sc, sd);
	CalculateDstRect(ht, wd, sa, sb, sc, sd, minx, miny, maxx, maxy);
	int dstHt = maxy - miny;
	int dstWd = maxx - minx;
	Mat warpimage(dstHt, dstWd, CV_8UC4, Scalar(0, 0, 0, 0));
	Mat mask(dstHt, dstWd, CV_8UC1, Scalar(0));
	ImageSimilarityWarp(mCurrentFrameImage, warpimage, sa, sb, sc, sd, minx, miny, maxx, maxy, mask);
	Mat roi(mosaic, Rect(minx + mwd*0.5, -miny + mht*0.5, warpimage.cols, warpimage.rows));
	warpimage.copyTo(roi, mask);
	////////////////////////////////////////////////////////////////

	return 0;
}

int CVideoMosaicSystem::LocalBA()
{
	//find the tracks seen by the last frame
	vector<int> currentTrackIndexNew;
	int nCurrentIndex = mKeyFrames.size() - 1;
	GetCurrentVisibleTracksNew(mTracks, mKeyFrames[nCurrentIndex].mFeature, currentTrackIndexNew);
	
	//calculate the initial value
	VisibleTracksInitialize(mKeyFrames, mTracks, currentTrackIndexNew);

	//




	return 0;
}


int CVideoMosaicSystem::Run(Mat& frameImage)
{

	mCurrentFrameImage.release();
	frameImage.copyTo(mCurrentFrameImage);

	//mCurrentFrameImage = frameImage.clone();
	//frameImage.clone();

	mbIsHavingKeyFrame = false;

	//initialize	
	if (!mbIsInitialized)
	{
		mbIsInitialized = Initialize();
		return 0;
	}
		
    //matching with the closest frame image
	if (IsKeyFrame())
	{
		CalculateNewCamPara();
		
		//LocalBA();
	}

	return 0;
}