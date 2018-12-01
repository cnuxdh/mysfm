
#include "warp.hpp"
#include "CalcAngle.h"
#include "Matrix.h"
#include "CommonFuncs.h"

#include"feature.hpp"
#include"sift.hpp"
#include"register.hpp"


#include <opencv\cv.h>
#include <opencv\highgui.h>

#include <opencv2/core/core.hpp> 
#include <opencv2/highgui/highgui.hpp> 
#include <opencv2/imgproc/imgproc.hpp> 
#include <opencv2/features2d/features2d.hpp>
#include "opencv2/nonfree/nonfree.hpp"


using namespace cv;


/*
int  AffineRansac(MyPointF* desPt, MyPointF* srcPt, int npt, double *ap)
{
    int ransac_rounds = 500;
    int maxInliers = 0;
    int DistanceThrshold = 16;

#define RANSAC_PT_NUM 4

    for (int round = 0; round < ransac_rounds; round++)
    {
        int indices[RANSAC_PT_NUM];
        ChooseRandIndex(npt, RANSAC_PT_NUM, indices);

        MyPointF dpts[RANSAC_PT_NUM];
        MyPointF spts[RANSAC_PT_NUM];
        for (int k = 0; k < RANSAC_PT_NUM; k++)
        {
            int randi = indices[k];
            dpts[k].x = desPt[randi].x;
            dpts[k].y = desPt[randi].y;
            spts[k].x = srcPt[randi].x;
            spts[k].y = srcPt[randi].y;
        }

        double p[6];
        SimilarityTransform(spts, dpts, RANSAC_PT_NUM, p);
        //float p[6];
        //AffineTransformFloat(spts, dpts, 3, p);

        //error
        double error = 0;
        int inliers = 0;
        //vInliers.clear();
        for (int i = 0; i < npt; i++)
        {
            double tx, ty;
            double dx = desPt[i].x;
            double dy = desPt[i].y;
            double ox = srcPt[i].x;
            double oy = srcPt[i].y;

            //tx = AFFINE_X(ox, oy, p);
            //ty = AFFINE_Y(ox, oy, p);

            tx = SIMILARITY_X(ox, oy, p);
            ty = SIMILARITY_Y(ox, oy, p);

            double len = sqrt((tx - dx)*(tx - dx) + (ty - dy)*(ty - dy));
            if (len < DistanceThrshold)
            {
                //vInliers.push_back(i);
                inliers++;
            }
            error += len;
        }
        error /= double(npt);

        if (maxInliers < inliers)
        {
            maxInliers = inliers;

            for (int k = 0; k < 6; k++) {
                ap[k] = p[k];
            }
        }
    }

    return maxInliers;

}
*/

//find homography using opencv
void  HomographyRansac(MyPointF* desPt, MyPointF* srcPt, int npt, double *hp) 
{
    //-- Localize the object
    std::vector<Point2f> despts;
    std::vector<Point2f> srcpts;

    for (int i = 0; i < npt; i++)
    {
        Point2f dp, sp;

        dp.x = desPt[i].x;
        dp.y = desPt[i].y;
        sp.x = srcPt[i].x;
        sp.y = srcPt[i].y;

        despts.push_back(dp);
        srcpts.push_back(sp);
    }
    
    
    Mat H = findHomography(srcpts, despts, CV_LMEDS, 8);
    int index = 0;
    for (int j = 0; j < 3; j++) {
        for (int i = 0; i < 3; i++) {
            hp[index] = H.at<double>(j, i);
            index++;
        }
    }
    

    /*
    Mat A = getAffineTransform(srcpts, despts);
    int index = 0;
    for (int j = 0; j < 2; j++) {
        for (int i = 0; i < 3; i++) {
            hp[index] = A.at<float>(j, i);
            index++;
        }
    }*/

    //return H;
}

void CalcWrapParas(char* dstfile, char*  srcfile, double* hp)
//void CalcWrapParas(IplImage* pDst, IplImage* pSrc, double* hp)
{
    //char* cfile = (char*)dstfile.data();
    
    CFeatureBase* pFeatDetect = new CSIFTFloat();

    ImgFeature dstfeats;
    ImgFeature srcfeats;
    pFeatDetect->Detect(dstfile, dstfeats, 640);
    pFeatDetect->Detect(srcfile, srcfeats, 640);

    CMatchBase* pMatch = new CSiftMatch();
    PairMatchRes mr;
    pMatch->Match(dstfeats, srcfeats, mr);
      
    printf("%d \n", mr.matchs.size());

    int nmatch = mr.matchs.size();
    MyPointF* pSrcPts = new MyPointF[nmatch];
    MyPointF* pDstPts = new MyPointF[nmatch];
    for (int i = 0; i < nmatch; i++) {
        int di = mr.matchs[i].l;
        int si = mr.matchs[i].r;
        pSrcPts[i].x = srcfeats.GetFeatPt(si).x;
        pSrcPts[i].y = srcfeats.GetFeatPt(si).y;
        pDstPts[i].x = dstfeats.GetFeatPt(di).x;
        pDstPts[i].y = dstfeats.GetFeatPt(di).y;
    }
    
    HomographyRansac(pDstPts, pSrcPts, nmatch, hp);


    delete pFeatDetect;
    delete pMatch;

    //Mat mSrc(pSrc);
    //Mat mDst(pDst);
    //vector<KeyPoint> kpSrc, kpDst; //keypoints
    //Mat dpSrc, dpDst; //descrptors
  
    /*
    ORB orb(2000);
    //resize the frame
    //double scale = (double)(mPanoPatch.rows) / (double)( mFrame.rows);
    //Size dsize = Size(mFrame.cols*scale, mFrame.rows*scale);
    //Mat resizeFrame = Mat(dsize, CV_8U);
    //resize(mFrame, resizeFrame, dsize);
    //imwrite("c:\\temp\\src.jpg", resizeFrame);
    //imwrite("c:\\temp\\dst.jpg", mFrame);
    //orb.detect(panoPatch, kpPanoPatch );
    //orb.detect(frame, kpFrame);
    orb(mSrc, Mat(), kpSrc, dpSrc);
    orb(mDst, Mat(), kpDst, dpDst);
    //BruteForceMatcher<HammingLUT> matcher;
    BFMatcher matcher(NORM_HAMMING);
    vector<DMatch> matches;
    matcher.match(dpDst, dpSrc, matches);
    */

    /*
    // detecting keypoints
    SurfFeatureDetector detector(400);
    detector.detect(mSrc, kpSrc);
    detector.detect(mDst, kpDst);

    // computing descriptors
    SurfDescriptorExtractor extractor;
    Mat descriptors1, descriptors2;
    extractor.compute(mSrc, kpSrc, dpSrc);
    extractor.compute(mDst, kpDst, dpDst);

    // matching descriptors
    BFMatcher matcher(NORM_L2);
    vector<DMatch> matches;
    matcher.match(dpDst, dpSrc, matches);

    printf("%d %d %d \n", dpDst.rows, dpSrc.rows, matches.size());
    */

    /*
    double max_dist = 0; double min_dist = 100;
    //Quick calculation of max and min distances between keypoints
    for (int i = 0; i < dpDst.rows; i++)
    {
        double dist = matches[i].distance;
        if (dist < min_dist) min_dist = dist;
        if (dist > max_dist) max_dist = dist;
    }
    printf("-- Max dist : %f \n", max_dist);
    printf("-- Min dist : %f \n", min_dist);
    //Draw only "good" matches (i.e. whose distance is less than 0.6*max_dist )
    // PS.- radiusMatch can also be used here.
    std::vector< DMatch > good_matches;
    for (int i = 0; i < dpDst.rows; i++)
    {
        if (matches[i].distance < 0.6*max_dist)
        {
            good_matches.push_back(matches[i]);
        }
    }
    printf("good match: %d \n", good_matches.size());
    Mat img_matches;
    drawMatches(mDst, kpDst, mSrc, kpSrc,
        good_matches, img_matches, Scalar::all(-1), Scalar::all(-1),
        vector<char>(), DrawMatchesFlags::NOT_DRAW_SINGLE_POINTS);
   
    imwrite("c:\\temp\\surf-match.jpg", img_matches);

    int nmatch = good_matches.size();
    MyPointF* pSrcPts = new MyPointF[nmatch];
    MyPointF* pDstPts = new MyPointF[nmatch];
    for (int i = 0; i < nmatch; i++) {
        int si = good_matches[i].queryIdx;
        int di = good_matches[i].trainIdx;
        pSrcPts[i].x = kpSrc[si].pt.x;
        pSrcPts[i].y = kpSrc[si].pt.y;
        pDstPts[i].x = kpDst[di].pt.x;
        pDstPts[i].y = kpDst[di].pt.y;
    }
    

    
    //double ap[9];
    //AffineRansac(pDstPts, pSrcPts, nmatch, ap);
    Mat hm = HomographyRansac(pDstPts, pSrcPts, nmatch, hp);    
    */


}


//generate the homography matrix: H = R - T.N'/d, x1=H.x2
//R,T:  rotation and translation, R(X - T)
//N',d: N'=(nx,ny,nz) is the plane normal, plane is described by: nx.X+ny.Y+nz.Z + distance=0
int GenerateHomographyMatrix(double omiga, double kapa, double phi, 
							 double tx, double ty, double tz,
							 double nx, double ny, double nz, 
							 double distance, double focal,
							 double* H )
{
	
	double R[9];
	double T[3];
	double pn[3];
	double K[9];

	for(int i=0; i<9; i++)
		K[i] = 0;
	K[0] = K[4] = 1;
	K[8] = -1.0/focal;

	GenerateRMatrixDirect(omiga, kapa, phi, R);

	T[0] = tx;
	T[1] = ty;
	T[2] = tz;

	double newT[3];
	mult(R, T, newT, 3, 3, 1);
	newT[0] = -newT[0];
	newT[1] = -newT[1];
	newT[2] = -newT[2];

	pn[0] = nx;
	pn[1] = ny;
	pn[2] = nz;
	
	double tn[9];

	mult(newT, pn, tn, 3, 1, 3);

	for(int i=0; i<9; i++)
	{
		H[i] = R[i] - tn[i]/distance;
	}

	double inversK[9];
	for(int i=0; i<9; i++)
		inversK[i] = K[i];
	invers_matrix(inversK, 3);
	
	double m1[9];
	mult(H, inversK, m1, 3, 3, 3);
	mult(K, m1, H, 3, 3, 3);

	invers_matrix(H, 3);

	return 0;
}



int HomographyWarp(unsigned char* pSrc, int sht, int swd, 
				   unsigned char* pDst, int& dht, int& dwd, 
				   double* H)
{
	for(int j=0; j<dht; j++)
		for(int i=0; i<dwd; i++)
		{
			double dx = double(i);
			double dy = double(j);

			//normalize the coordinates
			dx = dx - dwd*0.5;
			dy = dht*0.5 - dy;
						
			//homography transform
			double denominator = dx*H[6] + dy*H[7] + H[8];

			double sx = (dx*H[0] + dy*H[1] + H[2]) / denominator;
			double sy = (dx*H[3] + dy*H[4] + H[5]) / denominator;

			int ix = (int)(sx); 
			int iy = (int)(sy);

			ix += swd*0.5;
			iy = sht*0.5 - iy;
			
			if(ix>=0 && ix<swd && iy>=0 && iy<sht)
			{
				pDst[j*dwd+i] = pSrc[iy*swd+ix];
			}
		}

	return 0;
}


//IplImage* benchImage = cvLoadImage(benchfile);
//IplImage* videoImage = cvLoadImage(videofile);

/*
double hp[9];
CalcWrapParas(benchfile, videofile, hp);


Mat hm(3, 3, CV_64F);
int index = 0;
for (int j = 0; j < 3; j++) {
for (int i = 0; i < 3; i++) {
hm.at<double>(j, i) = hp[index];
index++;
}
}
*/


/*
Mat am(2, 3, CV_32FC1);
int index = 0;
for (int j = 0; j < 3; j++) {
for (int i = 0; i < 3; i++) {
am.at<double>(j, i) = hp[index];
index++;
}
}
*/