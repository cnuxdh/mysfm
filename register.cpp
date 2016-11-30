
#ifdef WIN32
#include "windows.h"
//#include "time.h"
#include "mmsystem.h"
#endif

//siftlib
#include "siftmatch.h"
#include "siftlib/ransac.h"

#include "register.hpp"
#include "relativepose.hpp"

#include <iostream>
#include <fstream>
#include <string>


//coredll
#include "main.h"

#include "baselib.h"


//matrix
//#include "matrix/matrix.h"


//imagelib
//#include "homography.h"
//#include "horn.h"
//#include "matrix.h"
//#include "tps.h"
//#include "vector.h"


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

		p[0] = k1[matches[i].l].x;
		p[1] = k1[matches[i].l].y;
		p[2] = 1.0;

		double q[3];
		dll_matrix_product(3, 3, 3, 1, M, p, q);

		double qx = q[0] / q[2];
		double qy = q[1] / q[2];

		double dx = qx - k2[matches[i].r].x;
		double dy = qy - k2[matches[i].r].y;

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
		dll_matrix_product(3, 3, 3, 1, M, l, rp);

		rp[0] /= rp[2];
		rp[1] /= rp[2];

		double dx = rp[0] - r[0];
		double dy = rp[1] - r[1];

		error += dx * dx + dy * dy;
	}

	//printf("[LeastSquaresFit] Residual error (before) is %0.3e\n", error);    


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
			dll_align_horn((int) inliers.size(), r_pts, l_pts, R, T, Tout, &scale, weight);
			memcpy(M, Tout, 9 * sizeof(double));
			break;
		}

	case MotionHomography: 
		{
			dll_align_homography((int) inliers.size(), r_pts, l_pts, M, 1);
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
		dll_matrix_product(3, 3, 3, 1, M, l, rp);

		rp[0] /= rp[2];
		rp[1] /= rp[2];

		double dx = rp[0] - r[0];
		double dy = rp[1] - r[1];

		error += dx * dx + dy * dy;
	}

	//printf("[LeastSquaresFit] Residual error (after) is %0.3e\n", error);    

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

			Vx(l_pts[i]) = k1[idx1].x;
			Vy(l_pts[i]) = k1[idx1].y;
			Vz(l_pts[i]) = 1.0;

			Vx(r_pts[i]) = k2[idx2].x;
			Vy(r_pts[i]) = k2[idx2].y;
			Vz(r_pts[i]) = 1.0;

			weight[i] = 1.0;
		}

		double Mcurr[9];

		switch (mm)
		{
		case MotionRigid: 
			{
				double R[9], T[9], Tout[9], scale;
				dll_align_horn(min_matches, r_pts, l_pts, R, T, Tout, &scale, weight);
				memcpy(Mcurr, Tout, 9 * sizeof(double));
				break;
			}

		case MotionHomography: 
			{
				dll_align_homography(min_matches, r_pts, l_pts, Mcurr, 1);
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


vector<int> EstimateFMatrix(  const vector<Point2DDouble>& pl, 
							  const vector<Point2DDouble>& pr,
							  int num_trials, 
							  double threshold,
							  double* fm)
{
	bool essential = false;

	int num_pts = (int) pl.size();

	/* num_pts should be greater than a threshold */
	if (num_pts < 16) 
	{
		std::vector<int> inliers;
		return inliers;
	}

	v3_t *k1_pts = new v3_t[num_pts];
	v3_t *k2_pts = new v3_t[num_pts];

	v3_t *k1_pts_in = new v3_t[num_pts];
	v3_t *k2_pts_in = new v3_t[num_pts];

	for (int i = 0; i < pl.size(); i++) 
	{
		k1_pts[i] = dll_v3_new(pl[i].p[0], pl[i].p[1], 1.0);
		k2_pts[i] = dll_v3_new(pr[i].p[0], pr[i].p[1], 1.0);
	}

	//double Fm[9];

	double F[9];	
	dll_estimate_fmatrix_ransac_matches(num_pts, k2_pts, k1_pts, 
		num_trials, threshold, 0.95, (essential ? 1 : 0), F);

	/* Find the inliers */
	std::vector<int> inliers;

	for (int i = 0; i < num_pts; i++) 
	{
		double dist = dll_fmatrix_compute_residual(F, k2_pts[i], k1_pts[i]);
		if (dist < threshold) 
		{
			inliers.push_back(i);
		}
	}

	/* Re-estimate using inliers */
	int num_inliers = (int) inliers.size();

	for (int i = 0; i < num_inliers; i++) 
	{
		k1_pts_in[i] = k1_pts[inliers[i]]; // v3_new(k1[idx1]->m_x, k1[idx1]->m_y, 1.0);
		k2_pts_in[i] = k2_pts[inliers[i]]; // v3_new(k2[idx2]->m_x, k2[idx2]->m_y, 1.0);
	}

	// printf("[1] num_inliers = %d\n", num_inliers);

#if 0
	double F0[9];
	double e1[3], e2[3];
	estimate_fmatrix_linear(num_inliers, k2_pts_in, k1_pts_in, F0, e1, e2);

	inliers.clear();
	for (int i = 0; i < num_pts; i++) {
		double dist = fmatrix_compute_residual(F0, k2_pts[i], k1_pts[i]);
		if (dist < threshold) {
			inliers.push_back(i);
		}
	}
	num_inliers = inliers.size();
	// printf("[2] num_inliers = %d\n", num_inliers);

	// matrix_print(3, 3, F0);
#else
	double F0[9];
	memcpy(F0, F, sizeof(double) * 9);
#endif

	if (!essential) 
	{
		/* Refine using NLLS */
		for (int i = 0; i < num_inliers; i++) 
		{
			k1_pts_in[i] = k1_pts[inliers[i]];
			k2_pts_in[i] = k2_pts[inliers[i]];
		}

		dll_refine_fmatrix_nonlinear_matches(num_inliers, k2_pts_in, k1_pts_in, 
			F0, F);
	}
	else 
	{
		memcpy(F, F0, sizeof(double) * 9);
	}

#if 0
	if (essential) 
	{
		/* Compute the SVD of F */
		double U[9], S[3], VT[9];
		dgesvd_driver(3, 3, F, U, S, VT);
		double E0[9] = { 1.0, 0.0, 0.0,
			0.0, 1.0, 0.0,
			0.0, 0.0, 0.0 };

		double tmp[9];
		matrix_product(3, 3, 3, 3, U, E0, tmp);
		matrix_product(3, 3, 3, 3, tmp, VT, F);
	}
#endif

	inliers.clear();
	for (int i = 0; i < num_pts; i++) 
	{
		double dist = dll_fmatrix_compute_residual(F, k2_pts[i], k1_pts[i]);
		if (dist < threshold) 
		{
			inliers.push_back(i);
		}
	}
	num_inliers = (int) inliers.size();

	delete [] k1_pts;
	delete [] k2_pts;
	delete [] k1_pts_in;
	delete [] k2_pts_in;


	if(fm!=NULL)
	{
		//memcpy(fm, F, sizeof(9*sizeof(double)));
		for(int i=0; i<9; i++)
			fm[i] = F[i];
	}

	return inliers;
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

#ifdef _WIN32
	//int64 t1 = getTickCount();
	//unsigned long t1 = timeGetTime();
#endif

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
	free(lPts);
	free(rPts);	

#ifdef _WIN32
	//int64 t2 = getTickCount();
	//unsigned long t2 = timeGetTime();
	//printf("match time: %lf \n", (double)(t2-t1) / getTickFrequency());
	//printf("match time: %lf \n", (double)(t2-t1) / 1000 );
#endif


	return 1;
}


int CKNNMatch::Match(ImgFeature& lImage, ImgFeature& rImage, PairMatchRes& pairMatch )
{
	if(lImage.featPts.size()<8 || rImage.featPts.size()<8)
		return 0;

#ifdef _WIN32
	//unsigned long t1 = timeGetTime();
#endif

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
		
	/*
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
	*/
	
	
	for(i=0; i<matchSet.mat_key_num; i++)
	{
		MatchPairIndex mid;
		mid.l = matchSet.matpoint[i].id1;
		mid.r = matchSet.matpoint[i].id2;
		pairMatch.matchs.push_back(mid);
	}
	

	delete[] matchSet.matpoint;
	free(lPts);
	free(rPts);

#ifdef _WIN32
	//unsigned long t2 = timeGetTime();
	//printf("match time: %lf s \n", (double)(t2-t1) / 1000 );
#endif

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
	#ifdef _WIN32
	//unsigned long t1 = timeGetTime();
	#endif	
	
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
	//printf("before duplicate removal: %d \n", nMatch);
	vector<int> mask( nMatch, 0 );
	for(int j=0; j<nMatch; j++)
	{
		if(mask[j]>0)
			continue;

		for(int i=j+1; i<nMatch; i++)
		{
			if( (lKeyIds[j]==lKeyIds[i]) || (rKeyIds[j]==rKeyIds[i]) )
			{
				mask[j]=1;
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
	//printf("after duplicate removal: %d \n", lKeyIds.size());

	//save result
	for(int i=0; i<lKeyIds.size(); i++)
	{
		MatchPairIndex mid;
		mid.l = lKeyIds[i];
		mid.r = rKeyIds[i];
		pairMatch.matchs.push_back(mid);
	}	
	printf("direct matching number : %d \n", pairMatch.matchs.size());
	

	//wrong match removal based on Homography
	/*
	if(1)
	{
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
	}*/

	if(1)
	{
		vector<Point2DDouble> lpts;
		vector<Point2DDouble> rpts;
		for(int i=0; i<pairMatch.matchs.size(); i++)
		{
			int lid = pairMatch.matchs[i].l;
			int rid = pairMatch.matchs[i].r;

			Point2DDouble lp,rp;
			lp.p[0] = lImage.featPts[lid].x;
			lp.p[1] = lImage.featPts[lid].y;
			
			rp.p[0] = rImage.featPts[rid].x;
			rp.p[1] = rImage.featPts[rid].y;

			lpts.push_back(lp);
			rpts.push_back(rp);
		}
	
		//for debug
		if(0)
		{
			char filename[256];
			sprintf(filename, "c:\\temp\\match-%d-%d.txt", pairMatch.lId, pairMatch.rId);
			FILE* fp = fopen(filename, "w");
			for(int ti=0; ti<lpts.size(); ti++)
			{
				fprintf(fp, "%.4lf %.4lf  %.4lf %.4lf \n", lpts[ti].p[0], lpts[ti].p[1],
					rpts[ti].p[0], rpts[ti].p[1]);
			}
			fclose(fp);
		}

		//from "Modeling the World from Internet Photo Collections"
		int threshold = 0.008*max( lImage.ht, lImage.wd );

		vector<int> inliers = EstimateFMatrix(lpts, rpts, 2048, threshold);

		vector<MatchPairIndex> inlierMatch;
		for(int i=0; i<inliers.size(); i++)
		{
			inlierMatch.push_back( pairMatch.matchs[ inliers[i] ] );
		}

		pairMatch.matchs = inlierMatch;
		pairMatch.inlierRatio = (double)(inlierMatch.size()) / double(lpts.size()); 
	}
	
	free(lKey);
	free(rKey);

	printf("Matching Point Number: %d \n", pairMatch.matchs.size());

#ifdef _WIN32
	//unsigned long t2 = timeGetTime();
	//printf("matching time: %lf s \n", (double)(t2-t1)/1000.0 );
#endif

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

int CMyGenerateTrack::GenerateTracks(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks )
{

	return 0;
}

int CMyGenerateTrack::GenerateTracks( vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks )
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
				track.extra = -1;
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


int CMyGenerateTrack::GenerateTracks( vector<PairMatchRes>& pairMatchs, vector<CImageDataBase*> imageData, vector<TrackInfo>& tracks )
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


CFastGenerateTrack::CFastGenerateTrack()
{

}

CFastGenerateTrack::~CFastGenerateTrack()
{

}

int CFastGenerateTrack::GenerateTracks(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks )
{
	int nImage = imgFeatures.size();

	//clear the "extra" value
	for(int i=0; i<nImage; i++)
	{
		int nFeat = imgFeatures[i].featPts.size();
		for(int j=0; j<nFeat; j++)
		{
			imgFeatures[i].featPts[j].extra = -1; //must be initialize as -1
		}
	}
	
	int nMatchPair = pairMatchs.size();

	tracks.clear();

	int nTrackIndex = 0;
	for(int i=0; i<nMatchPair; i++)
	{
		int nLeftImage  = pairMatchs[i].lId; //left image index
		int nRightImage = pairMatchs[i].rId; //right image index
		int nMatchNum   = pairMatchs[i].matchs.size();

		printf("match pair: %d - %d \n", nLeftImage, nRightImage);
		
		//process each match point pair
		for(int j=0; j<nMatchNum; j++)
		{
			int nLeftFeatIndex  = pairMatchs[i].matchs[j].l;
			int nRightFeatIndex = pairMatchs[i].matchs[j].r;

			int lTrackIndex = imgFeatures[nLeftImage].featPts[nLeftFeatIndex].extra;
			int rTrackIndex = imgFeatures[nRightImage].featPts[nRightFeatIndex].extra;

			TrackInfo tp;
			tp.extra = -1;

			ImageKey ik;
			
			//create a new track point
			if( lTrackIndex==-1 && rTrackIndex==-1 )
			{
				ik.first = nLeftImage;
				ik.second = nLeftFeatIndex;
				tp.views.push_back(ik);

				ik.first = nRightImage;
				ik.second = nRightFeatIndex;
				tp.views.push_back(ik);

				tp.id = nTrackIndex;

				//bridge the image point and track
				imgFeatures[nLeftImage].featPts[nLeftFeatIndex].extra   = nTrackIndex;
				imgFeatures[nRightImage].featPts[nRightFeatIndex].extra = nTrackIndex;

				tracks.push_back(tp);

				nTrackIndex ++;
			}
			///////////////////////////////////////////////////////////////////////
			
			//adding the match point into existing track point
			if(lTrackIndex != -1 && rTrackIndex == -1)
			{
				ik.first  = nRightImage;
				ik.second = nRightFeatIndex;
				tp.views.push_back(ik);
				
				//adding the right image point into the track
				tracks[lTrackIndex].views.push_back(ik);
				imgFeatures[nRightImage].featPts[nRightFeatIndex].extra = lTrackIndex;
			}
			if(lTrackIndex == -1 && rTrackIndex != -1)
			{
				ik.first  = nLeftImage;
				ik.second = nLeftFeatIndex;
				tp.views.push_back(ik);

				//adding the left image point into the track
				tracks[rTrackIndex].views.push_back(ik);
				imgFeatures[nLeftImage].featPts[nLeftFeatIndex].extra = rTrackIndex;
			}
			////////////////////////////////////////////////////////////////////////////
		}
	}

	//remove the track with multiple projections in the same image
	//vector<TrackInfo> newTracks;
	for(int i=0; i<tracks.size(); i++)
	{
		vector<int> vecNumberFreq;
		vecNumberFreq.resize(nImage, 0);
		bool bIsRemove = false;

		for(int j=0; j<tracks[i].views.size(); j++)
		{
			int nImageIndex = tracks[i].views[j].first;
			vecNumberFreq[nImageIndex]++;

			if(vecNumberFreq[nImageIndex]>1)
			{
				bIsRemove = true;
				break;
			}
		}
		//if(!bIsRemove)
		//	newTracks.push_back(tracks[i]);

		if(bIsRemove)
		{
			for(int k=0; k<vecNumberFreq.size(); k++)
				vecNumberFreq[k] = 0;

			ImageKeyVector newView;
			for(int j=0; j<tracks[i].views.size(); j++)
			{
				int nImageIndex = tracks[i].views[j].first;
				vecNumberFreq[nImageIndex]++;

				ImageKey ik = tracks[i].views[j];

				if(vecNumberFreq[nImageIndex]==1)
				{
					newView.push_back(ik);
				}
			}
			tracks[i].views = newView;
		}
	}

	//tracks = newTracks;

	return 0;
}


int CFastGenerateTrack::GenerateTracks( vector<PairMatchRes>& pairMatchs, vector<TrackInfo>& tracks )
{
	
	return 0;
}
int CFastGenerateTrack::GenerateTracks( vector<PairMatchRes>& pairMatchs, vector<CImageDataBase*> imageData, vector<TrackInfo>& tracks )
{
	return 0;
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
			int nPtId  = key.second;

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
				//key.first = 0;
				oneTrack.views.push_back(key);
			}
			if( nImgId==imgId2)
			{
				nFind++;
				//key.first = 1;
				oneTrack.views.push_back(key);
			}
		}

		if(nFind==2)
		{
			//save the new track index in the original tracks
			int currentTrackIndex = trackSeq.size();
			tracks[i].extra = currentTrackIndex;

			oneTrack.extra = i; //save the new track point index
			oneTrack.valid = 1; //initialized as 1
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


int LoadMatchFile( char* filename, vector<PairMatchRes>& pairMatchs)
{

	FILE *f = fopen(filename, "r");

	if (f == NULL) 
	{
		printf("[LoadMatchFile] Error opening file %s for reading\n", filename);	
		return -1;
	}

	char buf[256];
	while (fgets(buf, 256, f)) 
	{
		PairMatchRes pair;
		
		// Read the images
		int i1, i2;
		sscanf(buf, "%d %d\n", &i1, &i2);

		pair.lId = i1;
		pair.rId = i2;

		// Read the number of matches
		int nMatches;
		fscanf(f, "%d\n", &nMatches);

		// Read the matches
		for (int i = 0; i < nMatches; i++) 
		{
			int k1, k2;
			fscanf(f, "%d %d\n", &k1, &k2);

			MatchPairIndex mi;
			mi.l = k1;
			mi.r = k2;
			pair.matchs.push_back(mi);
		}
		pairMatchs.push_back(pair);
	}
	fclose(f);

	return 0;
}
