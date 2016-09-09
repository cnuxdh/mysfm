
#include "windows.h"
//#include "time.h"
#include "mmsystem.h"


#include "register.hpp"
#include "relativepose.hpp"

#include <iostream>
#include <fstream>
#include <string>


//coredll
#include "main.h"

//sift dll
#include "siftmatch.h"
#include "ransac.h"

//matrix
#include "matrix/matrix.h"


//imagelib
#include "homography.h"
#include "horn.h"
#include "matrix.h"
#include "tps.h"
#include "vector.h"


using namespace std;



static int CountInliers(const vector<stPtFeature> &k1, 
						const vector<stPtFeature> &k2, 
						vector<MatchPairIndex> matches,
						double *M, double thresh, vector<int> &inliers)
{
	inliers.clear();
	int count = 0;

	for (unsigned int i = 0; i < matches.size(); i++) 
	{
		/* Determine if the ith feature in f1, when transformed by M,
		* is within RANSACthresh of its match in f2 (if one exists)
		*
		* if so, increment count and append i to inliers */

		double p[3];

		p[0] = k1[matches[i].l].cx;
		p[1] = k1[matches[i].l].cy;
		p[2] = 1.0;

		double q[3];
		matrix_product(3, 3, 3, 1, M, p, q);

		double qx = q[0] / q[2];
		double qy = q[1] / q[2];

		double dx = qx - k2[matches[i].r].cx;
		double dy = qy - k2[matches[i].r].cy;

		double dist = sqrt(dx * dx + dy * dy);

		if (dist <= thresh) 
		{
			count++;
			inliers.push_back(i);
		}
	}

	return count;
}

static int LeastSquaresFit(const vector<stPtFeature> &k1, 
						   const vector<stPtFeature> &k2, 
						   vector<MatchPairIndex> matches, MotionModel mm,
						   const std::vector<int> &inliers, double *M)
{
	v3_t *r_pts = new v3_t[inliers.size()];
	v3_t *l_pts = new v3_t[inliers.size()];
	double *weight = new double[inliers.size()];

	/* Compute residual */
	double error = 0.0;
	for (int i = 0; i < (int) inliers.size(); i++) 
	{
		int idx1 = matches[inliers[i]].l;
		int idx2 = matches[inliers[i]].r;

		double r[3], l[3];
		l[0] = k1[idx1].cx;
		l[1] = k1[idx1].cy;
		l[2] = 1.0;

		r[0] = k2[idx2].cx;
		r[1] = k2[idx2].cy;
		r[2] = 1.0;	

		double rp[3];
		matrix_product(3, 3, 3, 1, M, l, rp);

		rp[0] /= rp[2];
		rp[1] /= rp[2];

		double dx = rp[0] - r[0];
		double dy = rp[1] - r[1];

		error += dx * dx + dy * dy;
	}

	printf("[LeastSquaresFit] Residual error (before) is %0.3e\n", error);    


	for (int i=0; i < (int) inliers.size(); i++) 
	{
		int idx1 = matches[inliers[i]].l;
		int idx2 = matches[inliers[i]].r;

		Vx(l_pts[i]) = k1[idx1].cx;
		Vy(l_pts[i]) = k1[idx1].cy;
		Vz(l_pts[i]) = 1.0;

		Vx(r_pts[i]) = k2[idx2].cx;
		Vy(r_pts[i]) = k2[idx2].cy;
		Vz(r_pts[i]) = 1.0;

		weight[i] = 1.0;
	}

	switch (mm) 
	{
	case MotionRigid: 
		{
			double R[9], T[9], Tout[9], scale;
			align_horn((int) inliers.size(), r_pts, l_pts, R, T, Tout, &scale, weight);
			memcpy(M, Tout, 9 * sizeof(double));
			break;
		}

	case MotionHomography: 
		{
			align_homography((int) inliers.size(), r_pts, l_pts, M, 1);
			break;
		}
	}

	/* Compute residual */
	error = 0.0;
	for (int i = 0; i < (int) inliers.size(); i++) 
	{
		int idx1 = matches[inliers[i]].l;
		int idx2 = matches[inliers[i]].r;

		double r[3], l[3];
		l[0] = k1[idx1].cx;
		l[1] = k1[idx1].cy;
		l[2] = 1.0;

		r[0] = k2[idx2].cx;
		r[1] = k2[idx2].cy;
		r[2] = 1.0;	

		double rp[3];
		matrix_product(3, 3, 3, 1, M, l, rp);

		rp[0] /= rp[2];
		rp[1] /= rp[2];

		double dx = rp[0] - r[0];
		double dy = rp[1] - r[1];

		error += dx * dx + dy * dy;
	}

	printf("[LeastSquaresFit] Residual error (after) is %0.3e\n", error);    

	delete [] r_pts;
	delete [] l_pts;
	delete [] weight;

	return 0;
}

/* Estimate a transform between two sets of keypoints */
vector<int> EstimateTransform(const vector<stPtFeature> &k1, 
							  const vector<stPtFeature> &k2, 
							  const vector<MatchPairIndex> &matches, 
							  MotionModel mm,
							  int nRANSAC, double RANSACthresh, 
							  double *Mout)
{
	int min_matches = -1;
	switch (mm) 
	{
	case MotionRigid:
		min_matches = 3;
		break;
	case MotionHomography:
		min_matches = 4;
		break;
	}

	int *match_idxs = new int[min_matches];

	int num_matches = (int) matches.size();
	int max_inliers = 0;
	double Mbest[9];

	if (num_matches < min_matches) 
	{
		std::vector<int> empty;
		printf("Cannot estimate rigid transform \n");
		return empty;
	}

	v3_t *r_pts = new v3_t[min_matches];
	v3_t *l_pts = new v3_t[min_matches];
	double *weight = new double[min_matches];

	for (int round = 0; round < nRANSAC; round++) 
	{
		for (int i = 0; i < min_matches; i++) 
		{
			bool found;
			int idx;

			do 
			{
				found = true;
				idx = rand() % num_matches;

				for (int j = 0; j < i; j++) 
				{
					if (match_idxs[j] == idx) 
					{
						found = false;
						break;
					}
				}
			} while (!found);

			match_idxs[i] = idx;
		}

		/* Solve for the motion */

		for (int i = 0; i < min_matches; i++) 
		{
			int idx1 = matches[match_idxs[i]].l;
			int idx2 = matches[match_idxs[i]].r;

			Vx(l_pts[i]) = k1[idx1].cx;
			Vy(l_pts[i]) = k1[idx1].cy;
			Vz(l_pts[i]) = 1.0;

			Vx(r_pts[i]) = k2[idx2].cx;
			Vy(r_pts[i]) = k2[idx2].cy;
			Vz(r_pts[i]) = 1.0;

			weight[i] = 1.0;
		}

		double Mcurr[9];

		switch (mm)
		{
		case MotionRigid: 
			{
				double R[9], T[9], Tout[9], scale;
				align_horn(min_matches, r_pts, l_pts, R, T, Tout, &scale, weight);
				memcpy(Mcurr, Tout, 9 * sizeof(double));
				break;
			}

		case MotionHomography: 
			{
				align_homography(min_matches, r_pts, l_pts, Mcurr, 1);
				break;
			}
		}

		std::vector<int> inliers;
		int num_inliers = CountInliers(k1, k2, matches, Mcurr, RANSACthresh, inliers);

		if (num_inliers > max_inliers) 
		{
			max_inliers = num_inliers;
			memcpy(Mbest, Mcurr, 9 * sizeof(double));
		}
	}

	std::vector<int> inliers;
	CountInliers(k1, k2, matches, Mbest, RANSACthresh, inliers);
	memcpy(Mout, Mbest, 9 * sizeof(double));
	LeastSquaresFit(k1, k2, matches, mm, inliers, Mout);

	delete [] match_idxs;
	delete [] r_pts;
	delete [] l_pts;
	delete [] weight;

	return inliers;
}

/*
//convert my data structure of point feature to OpenCV format
void FeatureConvert(vector<PtFeature>& imageFeats, ImageFeatures &features)
{
	int numPt = imageFeats.size();
	int featDim = imageFeats[0].feat.size();
	
	//features.keypoints.resize(numPt);
	features.descriptors.create(numPt, featDim, CV_8U);

	uchar* pBuffer = NULL;
	for(int i=0; i<numPt; i++)
	{
		KeyPoint p;
		p.pt.x = imageFeats[i].x;
		p.pt.y = imageFeats[i].y;
		features.keypoints.push_back(p);

		pBuffer = features.descriptors.ptr<uchar>(i);
		for(int j=0; j<featDim; j++)
		{
			pBuffer[j] =  imageFeats[i].feat[j];
		}
	}
}
*/


//convert my data format to SIFT_float format
void FeatureConvert(PtFeature srcFeat, Key_Point& dstFeat)
{
	dstFeat.index      = srcFeat.id;
	//dstFeat.key_column = srcFeat.x;
	//dstFeat.key_row    = srcFeat.y;
	dstFeat.initl_column = srcFeat.x;
	dstFeat.initl_row    = srcFeat.y;
	
	int nFeatDim = srcFeat.feat.size();

	for(int i=0; i<nFeatDim; i++)
	{
		dstFeat.descriptor[i] = srcFeat.feat[i];
	}
}


CRationMatch::CRationMatch()
{

}

CRationMatch::~CRationMatch()
{

}

int CRationMatch::Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats )
{	
	if( lImageFeats.size()<2 || rImageFeats.size()<2 )
		return 0;

	/*
	//based on OpenCV
    vector<ImageFeatures> features(2);
    
	//format conversion
	FeatureConvert(lImageFeats, features[0]);
	FeatureConvert(rImageFeats, features[1]);	
	
	int64 t = getTickCount();
	
	vector<MatchesInfo> pairwise_matches;
	float match_conf = 0.3;
	BestOf2NearestMatcher matcher(false, match_conf);
	matcher(features, pairwise_matches);
	matcher.collectGarbage();
	
	cout<<"Pairwise matching, time: " << ((getTickCount() - t) / getTickFrequency()) << " sec";

	//printf("%d \n", pairwise_matches.)
	*/

	return 1;
}

int CRationMatch::Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes)
{
	/*

	//based on OpenCV
	vector<ImageFeatures> features(2);

	features[0].img_idx = 0;
	features[1].img_idx = 1;
	features[0].img_size.height = lImage.ht;
	features[0].img_size.width  = lImage.wd;
	features[1].img_size.height = rImage.ht;
	features[1].img_size.width  = rImage.wd;

	//format conversion
	FeatureConvert(lImage.featPts, features[0]);
	FeatureConvert(rImage.featPts, features[1]);	
	
	int64 t = getTickCount();

	vector<MatchesInfo> pairwise_matches;
	float match_conf = 0.3;
	BestOf2NearestMatcher matcher(false, match_conf);
	matcher(features, pairwise_matches);
	matcher.collectGarbage();

	cout<<"Pairwise matching, time: " << ((getTickCount() - t) / getTickFrequency()) << " sec";

	//printf("%d \n", pairwise_matches.)
	*/

	return 1;
}

int CRationMatch::Match(char* file1, char* file2)
{

	/*
	Mat img1 = imread(file1);
	Mat img2 = imread(file2);

    vector<ImageFeatures> features(2);

	Ptr<FeaturesFinder> finder = new SurfFeaturesFinder();

	(*finder)(img1, features[0]);
	features[0].img_idx = 0;

	(*finder)(img2, features[1]);
	features[1].img_idx = 1;

	int64 t = getTickCount();

	vector<MatchesInfo> pairwise_matches;
	float match_conf = 0.3;
	BestOf2NearestMatcher matcher(false, match_conf);
	matcher(features, pairwise_matches);
	matcher.collectGarbage();

	cout<<"Pairwise matching, time: " << ((getTickCount() - t) / getTickFrequency()) << " sec";

	finder->collectGarbage();
	*/

	return 1;
}
//////////////////////////////////////////////////////////////////////////




CPanoMatch::CPanoMatch()
{

}

CPanoMatch::~CPanoMatch()
{

}

int CPanoMatch::Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats )
{
	return 1;
}

int CPanoMatch::Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes)
{
	if(lImage.featPts.size()<8 || rImage.featPts.size()<8)
		return 0;

	int nDim = lImage.featPts[0].feat.size();

	//data format conversion
	int lPtNum = lImage.featPts.size();
	int rPtNum = rImage.featPts.size();

	Key_Point* lPts = (Key_Point*)malloc(lPtNum*sizeof(Key_Point));
	Key_Point* rPts = (Key_Point*)malloc(rPtNum*sizeof(Key_Point));

	int i = 0;
	for(i=0; i<lPtNum; i++)
	{
		FeatureConvert( lImage.featPts[i], lPts[i] );
	}
	for(i=0; i<rPtNum; i++)
	{
		FeatureConvert( rImage.featPts[i], rPts[i] );
	}

	MatchPoint* pMatch;
	int nMatch = 0;
	SiftPairMatch(lPts, lPtNum, rPts, rPtNum, &pMatch, &nMatch, nDim);

	for(i=0; i<nMatch; i++)
	{
		MatchPairIndex mid;
		mid.l = pMatch[i].id1;
		mid.r = pMatch[i].id2;
		matchRes.push_back(mid);
	}
    delete[] pMatch;
	return 1;
}
int CPanoMatch::Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch )
{

	return 0;
}


//////////////////////////////////////////////////////////////////////////
CKNNMatch::CKNNMatch()
{

}

CKNNMatch::~CKNNMatch()
{

}

int CKNNMatch::Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats )
{


	return 1;
}

int CKNNMatch::Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes)
{
	if(lImage.featPts.size()<8 || rImage.featPts.size()<8)
		return 0;

	//int64 t1 = getTickCount();
	unsigned long t1 = timeGetTime();

	int nDim = lImage.featPts[0].feat.size();

	//data format conversion
	int lPtNum = lImage.featPts.size();
	int rPtNum = rImage.featPts.size();
	
	Key_Point* lPts = (Key_Point*)malloc(lPtNum*sizeof(Key_Point));
	Key_Point* rPts = (Key_Point*)malloc(rPtNum*sizeof(Key_Point));

	int i = 0;
	for(i=0; i<lPtNum; i++)
	{
		FeatureConvert( lImage.featPts[i], lPts[i] );
	}
	for(i=0; i<rPtNum; i++)
	{
		FeatureConvert( rImage.featPts[i], rPts[i] );
	}

	MatchPoint* pMatch;
	int nMatch = 0;
	SiftPairMatch(lPts, lPtNum, rPts, rPtNum, &pMatch, &nMatch, nDim);
	
	//error removal
    MatchPoints matchSet;
    matchSet.mat_key_num = nMatch;
	matchSet.matpoint = pMatch;
    
	
	MatchPoints goodMatch = del_gross_error(matchSet, lImage.wd, lImage.ht);
 
	//save result
	for(i=0; i<goodMatch.mat_key_num; i++)
	{
		MatchPairIndex mid;
		mid.l = goodMatch.matpoint[i].id1;
		mid.r = goodMatch.matpoint[i].id2;
		matchRes.push_back(mid);
	}
	delete[] goodMatch.matpoint;
   
	delete[] matchSet.matpoint;

	//int64 t2 = getTickCount();
	unsigned long t2 = timeGetTime();

	//printf("match time: %lf \n", (double)(t2-t1) / getTickFrequency());
	printf("match time: %lf \n", (double)(t2-t1) / 1000 );


	return 1;
}


int CKNNMatch::Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch )
{
	if(lImage.featPts.size()<8 || rImage.featPts.size()<8)
		return 0;


	unsigned long t1 = timeGetTime();

	int nDim = lImage.featPts[0].feat.size();

	//data format conversion
	int lPtNum = lImage.featPts.size();
	int rPtNum = rImage.featPts.size();

	Key_Point* lPts = (Key_Point*)malloc(lPtNum*sizeof(Key_Point));
	Key_Point* rPts = (Key_Point*)malloc(rPtNum*sizeof(Key_Point));

	int i = 0;
	for(i=0; i<lPtNum; i++)
	{
		FeatureConvert( lImage.featPts[i], lPts[i] );
	}
	for(i=0; i<rPtNum; i++)
	{
		FeatureConvert( rImage.featPts[i], rPts[i] );
	}

	MatchPoint* pMatch;
	int nMatch = 0;
	SiftPairMatch(lPts, lPtNum, rPts, rPtNum, &pMatch, &nMatch, nDim);

	//error removal
	MatchPoints matchSet;
	matchSet.mat_key_num = nMatch;
	matchSet.matpoint = pMatch;
		
	
	MatchPoints goodMatch = del_gross_error(matchSet, lImage.wd, lImage.ht);
	//inlier ratio
	pairMatch.inlierRatio =  (double)(goodMatch.mat_key_num)/(double)(matchSet.mat_key_num);
	//save result
	for(i=0; i<goodMatch.mat_key_num; i++)
	{
		MatchPairIndex mid;
		mid.l = goodMatch.matpoint[i].id1;
		mid.r = goodMatch.matpoint[i].id2;
		pairMatch.matchs.push_back(mid);
	}
	delete[] goodMatch.matpoint;
	


	for(i=0; i<matchSet.mat_key_num; i++)
	{
		MatchPairIndex mid;
		mid.l = matchSet.matpoint[i].id1;
		mid.r = matchSet.matpoint[i].id2;
		pairMatch.matchs.push_back(mid);
	}

	delete[] matchSet.matpoint;

	unsigned long t2 = timeGetTime();
	printf("match time: %lf s \n", (double)(t2-t1) / 1000 );

	return 1;	
}
//////////////////////////////////////////////////////////////////////////




int ConvertFeatToUChar(vector<stPtFeature> feats, unsigned char** key)
{
	int nFeatPtNum = feats.size();
	int nFeatDim   = feats[0].feat.size();
	*key = (unsigned char*)malloc(nFeatPtNum*nFeatDim);
	
	int index = 0;
	for(int i=0; i<nFeatPtNum; i++)
	{
		for(int j=0; j<nFeatDim; j++)
		{
			(*key)[index] = feats[i].feat[j];
			index ++;
		}
	}

	return nFeatPtNum;
}


CSiftMatch::CSiftMatch()
{

}
CSiftMatch::~CSiftMatch()
{
}

int CSiftMatch::Match(vector<PtFeature> lImageFeats, vector<PtFeature> rImageFeats )
{
	return 1;
}

int CSiftMatch::Match(ImgFeature& lImage, ImgFeature& rImage, vector<MatchPairIndex>& matchRes)
{
	return 1;
}   



//two images matching
int CSiftMatch::Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch )
{
	unsigned long t1 = timeGetTime();

	//feature convert
	unsigned char* lKey = NULL;
	int nLeftKey  = ConvertFeatToUChar(lImage.featPts, &lKey);
	unsigned char* rKey = NULL;
	int nRightKey = ConvertFeatToUChar(rImage.featPts, &rKey);

	vector<int> lKeyIds;
	vector<int> rKeyIds;
	PairMatchUsingSiftKey(lKey, nLeftKey, rKey, nRightKey, lKeyIds, rKeyIds);

	//duplicate removal
	int nMatch = lKeyIds.size();
	vector<int> mask( nMatch, 0 );
	for(int j=0; j<nMatch; j++)
	{
		for(int i=j+1; i<nMatch; i++)
		{
			if( rKeyIds[j]==rKeyIds[i] )
			{
				mask[i]=1;
			}
		}
	}
	vector<int> lKeyIdsNew;
	vector<int> rKeyIdsNew;
	for(int j=0; j<nMatch; j++)
	{
		if(mask[j]==0)
		{
			lKeyIdsNew.push_back( lKeyIds[j] );
			rKeyIdsNew.push_back( rKeyIds[j] );
		}
	}
	lKeyIds = lKeyIdsNew;
	rKeyIds = rKeyIdsNew;

	//save result
	for(int i=0; i<lKeyIds.size(); i++)
	{
		MatchPairIndex mid;
		mid.l = lKeyIds[i];
		mid.r = rKeyIds[i];
		pairMatch.matchs.push_back(mid);
	}	
	
	//wrong match removal based on Homography
	double dH[9];
	vector<int> inliers = EstimateTransform( lImage.featPts, rImage.featPts, pairMatch.matchs,
											 MotionHomography, HOMOGRAPHY_ROUNDS, HOMOGRAPHY_THRESHOLD, dH);
  	
	
	//inlier ratio
	pairMatch.inlierRatio = (double)(inliers.size()) / (double)( pairMatch.matchs.size() );
	//pairMatch.inlierRatio =  (double)(lKeyIds.size())/(double)();

	vector<MatchPairIndex> rightMatch;
	for(int i=0; i<inliers.size(); i++)
	{
		rightMatch.push_back( pairMatch.matchs[ inliers[i] ] );
	}
	pairMatch.matchs = rightMatch;
	

	free(lKey);
	free(rKey);

	printf("Matching Point Number: %d \n", pairMatch.matchs.size());

	unsigned long t2 = timeGetTime();
	printf("matching time: %lf s \n", (double)(t2-t1)/1000.0 );

	return 1;
}


//two images matching
int CSiftMatch::Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch, double* m )
{
	//feature convert
	unsigned char* lKey = NULL;
	int nLeftKey  = ConvertFeatToUChar(lImage.featPts, &lKey);
	unsigned char* rKey = NULL;
	int nRightKey = ConvertFeatToUChar(rImage.featPts, &rKey);

	vector<int> lKeyIds;
	vector<int> rKeyIds;
	PairMatchUsingSiftKey(lKey, nLeftKey, rKey, nRightKey, lKeyIds, rKeyIds);

	//duplicate removal
	int nMatch = lKeyIds.size();
	vector<int> mask( nMatch, 0 );
	for(int j=0; j<nMatch; j++)
	{
		for(int i=j+1; i<nMatch; i++)
		{
			if( rKeyIds[j]==rKeyIds[i] )
			{
				mask[i]=1;
			}
		}
	}
	vector<int> lKeyIdsNew;
	vector<int> rKeyIdsNew;
	for(int j=0; j<nMatch; j++)
	{
		if(mask[j]==0)
		{
			lKeyIdsNew.push_back( lKeyIds[j] );
			rKeyIdsNew.push_back( rKeyIds[j] );
		}
	}
	lKeyIds = lKeyIdsNew;
	rKeyIds = rKeyIdsNew;

	//save result
	for(int i=0; i<lKeyIds.size(); i++)
	{
		MatchPairIndex mid;
		mid.l = lKeyIds[i];
		mid.r = rKeyIds[i];
		pairMatch.matchs.push_back(mid);
	}

	//wrong match removal based on Homography
	double dH[9];
	vector<int> inliers = EstimateTransform( lImage.featPts, rImage.featPts, pairMatch.matchs,
		MotionHomography, HOMOGRAPHY_ROUNDS, HOMOGRAPHY_THRESHOLD, dH);
	
	//save Homography
	for(int i=0; i<9; i++)
		m[i] = dH[i];
	
	//inlier ratio
	pairMatch.inlierRatio = (double)(inliers.size()) / (double)( pairMatch.matchs.size() );
	//pairMatch.inlierRatio =  (double)(lKeyIds.size())/(double)();

	vector<MatchPairIndex> rightMatch;
	for(int i=0; i<inliers.size(); i++)
	{
		rightMatch.push_back( pairMatch.matchs[ inliers[i] ] );
	}

	pairMatch.matchs = rightMatch;


	free(lKey);
	free(rKey);

	return 1;
}


//////////////////////////////////////////////////////////////////////////
CMyGenerateTrack::CMyGenerateTrack()
{

}

CMyGenerateTrack::~CMyGenerateTrack()
{

}

int CMyGenerateTrack::GenerateTracks( vector<PairMatchRes> pairMatchs, vector<TrackInfo>& tracks )
{
	int nTrack = 0;

	for(int i=0; i<pairMatchs.size(); i++)
	{
		int nMatch = pairMatchs[i].matchs.size();
		ImageKey lkey,rkey;
		ImageKey newKey;
		for(int k=0; k<nMatch; k++ )
		{			
			lkey.first = pairMatchs[i].lId;
			lkey.second = pairMatchs[i].matchs[k].l;
			rkey.first = pairMatchs[i].rId;
			rkey.second = pairMatchs[i].matchs[k].r;

			bool bIsNew = true;
			int  iKey = 0;
			int  nMatch = 0; //determine if the match pair has been added, 1:no, 2:yes
			for(int j=0; j<tracks.size(); j++)
			{
				//check if the current projection belong to existing track
				for(int m=0; m<tracks[j].views.size(); m++)
				{
					ImageKey key = tracks[j].views[m];
					if(key.first == lkey.first && key.second == lkey.second)
					{
						bIsNew = false;
						iKey = j;
						newKey = rkey;
						nMatch++;
						//break;
					}
					if(key.first == rkey.first && key.second == rkey.second)
					{
						bIsNew = false;
						iKey = j;
						newKey = lkey;
						nMatch++;
						//break;
					}
					//if(!bIsNew)
					//	break;
				}
				if(!bIsNew)
					break;
			}
			if(bIsNew)
			{
				TrackInfo track;
				track.id = nTrack;
				track.views.push_back(lkey);
				track.views.push_back(rkey);
				tracks.push_back(track);
				nTrack++;
			}
			else
			{
				//add existing track
				if(nMatch==1)
					tracks[iKey].views.push_back(newKey);
			}
		}
	}

	return 1;
}


int CMyGenerateTrack::GenerateTracks( vector<PairMatchRes> pairMatchs, vector<CImageDataBase*> imageData, vector<TrackInfo>& tracks )
{
	int nTrack = 0;

	for(int i=0; i<pairMatchs.size(); i++)
	{
		int nMatch = pairMatchs[i].matchs.size();
		ImageKey lkey,rkey;
		ImageKey newKey;
		for(int k=0; k<nMatch; k++ )
		{			
			lkey.first  = pairMatchs[i].lId;
			lkey.second = pairMatchs[i].matchs[k].l;
			rkey.first  = pairMatchs[i].rId;
			rkey.second = pairMatchs[i].matchs[k].r;

			bool bIsNew = true;
			int  iKey = 0;
			int  nMatch = 0; //determine if the match pair has been added, 1:no, 2:yes
			for(int j=0; j<tracks.size(); j++)
			{
				//check if the current projection belong to existing track
				for(int m=0; m<tracks[j].views.size(); m++)
				{
					ImageKey key = tracks[j].views[m];
					if(key.first == lkey.first && key.second == lkey.second)
					{
						bIsNew = false;
						iKey = j;
						newKey = rkey;
						nMatch++;
						//break;
					}
					if(key.first == rkey.first && key.second == rkey.second)
					{
						bIsNew = false;
						iKey = j;
						newKey = lkey;
						nMatch++;
						//break;
					}
					//if(!bIsNew)
					//	break;
				}
				if(!bIsNew)
					break;
			}
			if(bIsNew)
			{
				TrackInfo track;
				track.id = nTrack;
				track.extra = -1;
				track.views.push_back(lkey);
				track.views.push_back(rkey);
				tracks.push_back(track);
				
				//imageData[lkey.first]->SetTrackIdx(lkey.second, nTrack);
				//imageData[rkey.first]->SetTrackIdx(rkey.second, nTrack);
				
				nTrack++;
			}
			else
			{
				//add existing track
				if(nMatch==1)
				{
					tracks[iKey].views.push_back(newKey);
					//imageData[newKey.first]->SetTrackIdx(newKey.second, iKey);
				}
			}
		}
	}

	//prund bad tracks
	PruneBadTracks(tracks);

	//set track id for each projected point of image
	for(int i=0; i<tracks.size(); i++)
	{
		int nView = tracks[i].views.size();
		for(int m=0; m<nView; m++)
		{
			int imgId = tracks[i].views[m].first;
			int ptId  = tracks[i].views[m].second;
			imageData[imgId]->SetTrackIdx(ptId, i);
		}
	}

	return 1;
}


//////////////////////////////////////////////////////////////////////////
//functions based on tracks

int      GetFileRows1(char* file)
{
	char  sline[512];
	FILE* fp = NULL;
	int rows = 0;

	fp = fopen(file, "r");

	if(fp==NULL)
		return 0;

	while( fgets(sline, 512, fp)!=NULL )
		rows++;

	fclose(fp);

	return rows;
}

int ReadTracks(char* filename, vector<TrackInfo>& tracks)
{	
	int nRow = GetFileRows1(filename);
	tracks.resize(nRow);

	FILE* fp = fopen(filename, "r");
	if(fp==NULL) 
		return 0;

	for(int i=0; i<nRow; i++)
	{
		TrackInfo oneTrack;
		int t=0;
		int npt=0;
		fscanf(fp, "%d %d", &t, &npt);
		int imgId,ptId;
		ImageKey key;
		for(int k=0; k<npt; k++)
		{
			fscanf(fp, "%d_%d ", &imgId, &ptId);
			key.first = imgId;
			key.second = ptId;
			oneTrack.views.push_back(key);
		}
		tracks[i] = oneTrack;
	}
	fclose(fp);
	
	return 0;
}


int PrintTracks(vector<TrackInfo> tracks, char* filename)
{
	FILE* fp = fopen(filename, "w");
	if(fp==NULL)
		return 0;

	for(int i=0; i<tracks.size(); i++)
	{
		fprintf(fp, "%d  %d ", i, tracks[i].views.size());
		for(int j=0; j<tracks[i].views.size(); j++)
		{
			ImageKey key = tracks[i].views[j]; //key for image point
			int nImgId = key.first;
			int nPtId = key.second;

			fprintf(fp," %d_%d ", nImgId, nPtId);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 1;
}


void GetMatch(int imgId1, int imgId2, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq)
{
	int nTracks = tracks.size();
	for(int i=0; i<nTracks; i++)
	{
		//check each image key to find if there is the image pair
		int nFind = 0;
		int nViews = tracks[i].views.size();	
		TrackInfo oneTrack;
		for(int j=0; j<nViews; j++)
		{
			ImageKey key = tracks[i].views[j]; //key for image point
			int nImgId = key.first;
			int nPtId = key.second;

			/*
			if( nImgId==imgId1 || nImgId==imgId2)
			{
				nFind++;
				oneTrack.views.push_back(key);
			}
			*/
			if( nImgId==imgId1)
			{
				nFind++;
				key.first = 0;
				oneTrack.views.push_back(key);
			}
			if( nImgId==imgId2)
			{
				nFind++;
				key.first = 1;
				oneTrack.views.push_back(key);
			}
		}
		if(nFind==2)
		{
			oneTrack.extra = i;
			trackSeq.push_back(oneTrack);
		}
	}
}

int FindProjections(int imgId, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq)
{
	int nProj = 0;    
	for(int i=0; i<trackSeq.size(); i++)
	{
		int trackId = trackSeq[i].extra;
		int nView = tracks[trackId].views.size();
		for(int j=0; j<nView; j++)
		{
			ImageKey key = tracks[trackId].views[j];
			if(key.first == imgId)
			{
				nProj++;
			}
		}
	}
	return nProj;
}

int FindNewImage(int numCam, vector<int> cameraIDOrder, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq)
{
	int imgId = 0;
    
	vector<int> projHist;
    projHist.resize(numCam);
	for(int i=0; i<numCam; i++)
		projHist[i] = 0;

	for(int i=0; i<trackSeq.size(); i++)
	{
		int trackId = trackSeq[i].extra;

		int nView = tracks[trackId].views.size();
		for(int j=0; j<nView; j++)
		{
			ImageKey key = tracks[trackId].views[j];
			int nImgId = key.first;
			
			bool bIsAdded = false;
			for(int j=0; j<cameraIDOrder.size(); j++)
			{
				if( nImgId== cameraIDOrder[j] )
					bIsAdded = true;
			}
			if(bIsAdded)
				continue;

			projHist[ key.first ]++;
		}
	}

	int nMaxProj = 0;
	int nNewCamId = 0;
	for(int i=0; i<projHist.size(); i++)
	{
		if(nMaxProj<projHist[i])
		{
			nMaxProj = projHist[i];
			nNewCamId = i;
		}
	}

	return nNewCamId;
}


int UpdateTracks(int imgId, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq)
{
	for(int i=0; i<trackSeq.size(); i++)
	{
		int trackId = trackSeq[i].extra;
		int nView = tracks[trackId].views.size();
		for(int j=0; j<nView; j++)
		{
			ImageKey key = tracks[trackId].views[j];
			int nImgId = key.first;
			if(imgId == nImgId)
			{
				trackSeq[i].views.push_back(key);
			}
		}
	}
	return 1;
}

int UpdateTracks(int imgId,  int newId, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq, 
				 vector<ImgFeature>& imgFeatures, vector<Point3DDouble>& pt3, vector<Point2DDouble>& pt2)
{
	for(int i=0; i<trackSeq.size(); i++)
	{
		int trackId = trackSeq[i].extra;

		int nView = tracks[trackId].views.size();
		for(int j=0; j<nView; j++)
		{
			ImageKey key = tracks[trackId].views[j];
			int nImgId = key.first;
			int nPtId  = key.second;

			if(imgId == nImgId)
			{
				key.first = newId;
				trackSeq[i].views.push_back(key);
				Point3DDouble p3;
				Point2DDouble p2;
				p3 = trackSeq[i].grd;
				p2.p[0] = imgFeatures[imgId].featPts[nPtId].cx;
				p2.p[1] = imgFeatures[imgId].featPts[nPtId].cy;
				pt3.push_back(p3);
				pt2.push_back(p2);
			}
		}
	}

	return 1;
}


/*
   imgId:         image ID to be added
   cameraIDOrder: images having been added
   tracks:        all tracks
   trackSeq:      current tracks
*/
int AddNewTracks(int imgId, vector<int> cameraIDOrder, vector<CameraPara> cameras,
				 vector<ImgFeature>& imgFeatures, vector<TrackInfo>& tracks, vector<TrackInfo>& trackSeq)
{
	CTriangulateBase* pTri = new CTriangulateCV();

	int nAllTracks = tracks.size();
	vector<int> trackIndex;
	trackIndex.resize(nAllTracks);
	for(int i=0; i<trackIndex.size(); i++)
		trackIndex[i] = 0;
    
	//label those tracks whose 3D is recovered
	for(int i=0; i<trackSeq.size(); i++)
	{
		int tid = trackSeq[i].extra;
		trackIndex[tid] = 1;
	}

	//find the match (imgId, cid ) 
	for(int j=0; j<cameraIDOrder.size(); j++)
	{
		int cid = cameraIDOrder[j];
		for( int i=0; i<tracks.size(); i++)
		{
			if( trackIndex[i]>0 )
				continue;
			
			TrackInfo oneTrack;
			oneTrack.extra = i;
			int nView = tracks[i].views.size();
			int nHit = 0;
			for(int k=0; k<nView; k++)
			{
				ImageKey key = tracks[i].views[k];
				if(  key.first==cid || key.first==imgId )
				{
					oneTrack.views.push_back(key);
					nHit++;
				}
			}

			if(nHit==2)
			{
				//recover the 3D point
				vector<Point2DDouble> lpts,rpts;
				vector<Point3DDouble> gpts;
				lpts.resize(1);
				rpts.resize(1);
				int lImgId = oneTrack.views[0].first;
				int lPtId  = oneTrack.views[0].second;
				int rImgId = oneTrack.views[1].first;
				int rPtId  = oneTrack.views[1].second;
				lpts[0].p[0] = imgFeatures[lImgId].featPts[lPtId].cx;
				lpts[0].p[1] = imgFeatures[lImgId].featPts[lPtId].cy;
				rpts[0].p[0] = imgFeatures[rImgId].featPts[rPtId].cx;
				rpts[0].p[1] = imgFeatures[rImgId].featPts[rPtId].cy;

				pTri->Triangulate(lpts, rpts, cameras[lImgId], cameras[rImgId], gpts);
				oneTrack.grd = gpts[0];
				trackSeq.push_back(oneTrack);
			}
		}
	}

	delete pTri;
	return 1;
}


void PruneBadTracks(vector<TrackInfo>& trackSeq)
{
	vector<TrackInfo> trackSeqRight;

	for(int i=0; i<trackSeq.size(); i++)
	{
		bool bIsBad = false;
		int nView = trackSeq[i].views.size();
		for(int m=0; m<nView; m++)
			for(int n=m+1; n<nView; n++)
			{
				if( trackSeq[i].views[m].first == trackSeq[i].views[n].first )
				{
					bIsBad = true;
					break;
				}
				if(bIsBad)
					break;
			}
			if(!bIsBad)
			{
				trackSeqRight.push_back( trackSeq[i] );
			}
	}

	trackSeq = trackSeqRight;
}

