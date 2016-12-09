
#include"stdio.h"

#include"cvInvoke.hpp"
#include"sift.hpp"
#include"dataBase.hpp"
#include"register.hpp"
#include"relativepose.hpp"
#include "CalcAngle.h"


//corelib
#include"commonfile.h"



// function for feature detection of multiple images, for input files
// Input
//	filenames : the image files
//  nFile :
// Output 
//	outpath: image feature files saved into the output path
int DetectFileFeaturePts(char** filenames, int nFile, char* outpath)
{
	//retrive the file title
	
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	CPointFeatureBase* pFeatureData = new CSiftFeatureDataBinary();
	
	for(int i=0; i<nFile; i++)
	{
		printf("image: %s \n", filenames[i]);
		
		ImgFeature feats;
		pFeatDetect->Detect(filenames[i], feats);
		
		//save the feat into the file
		//printf("save the feature points into the file.... \n");
		char* title;
		GetTitleName(filenames[i], &title);
		printf("title: %s \n", title);
		
		char featPath[512];
		sprintf(featPath, "%s/%s.dat", outpath, title);
		printf("%s \n", featPath);
		
		pFeatureData->Write(featPath, feats);		
	}	
	delete pFeatDetect;
	
	return 0;
}


int DLL_EXPORT DetectFileFeaturePts(char* filenames, ImgFeature& imgFeatures, int maxHt)
{
	CFeatureBase* pFeatDetect = new CSIFTFloat();

	pFeatDetect->Detect(filenames, imgFeatures, maxHt);

	delete pFeatDetect;

	return 0;
}



//interface for memory
int DetectFileFeaturePts(char** filenames, int nFile, vector<ImgFeature>& imgFeatures, int maxHt)
{

	printf("[DetectFileFeaturePts] ... \n");

	CFeatureBase* pFeatDetect = new CSIFTFloat();
	//CPointFeatureBase* pFeatureData = new CSiftFeatureDataBinary();

	for(int i=0; i<nFile; i++)
	{
		printf("image: %s \n", filenames[i]);
		
		ImgFeature feats;
		pFeatDetect->Detect(filenames[i], feats, maxHt);
		imgFeatures.push_back(feats);
	}
	
	delete pFeatDetect;
	return 0;
}


//matching between image files
int MatchImageFiles(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& matchRes)
{
	printf("[MatchImageFiles] ... \n");

	int nImageNum = imgFeatures.size();
		
	CMatchBase* pMatch = new CSiftMatch();  //new CKNNMatch();

	for(int i=0; i<nImageNum; i++)
		for(int j=i+1; j<nImageNum; j++)
		{

			PairMatchRes mr;
			mr.lId = i;
			mr.rId = j;

			pMatch->Match(imgFeatures[i], imgFeatures[j], mr);

			//printf("%d-%d  %d %lf \n", i, j, mr.matchs.size(), mr.inlierRatio );
			printf("%d-%d  %d  \n", i, j, mr.matchs.size() );
			
			matchRes.push_back(mr);			
		}
	
	delete pMatch;

	return 0;
}


DLL_EXPORT int dll_EstimatePose5Point_Pano( vector<Point3DDouble>& pl, 
	vector<Point3DDouble>& pr,
	double radius,
	int num_trials, double threshold, 
	double *R, double *t, vector<double>& residual)
{
	return EstimatePose5Point_Pano(pl, pr, radius, num_trials, threshold, R, t, residual);
}


DLL_EXPORT int dll_EstimatePose( vector<Point2DDouble> lPts, vector<Point2DDouble> rPts,
	CameraPara& cam1, CameraPara& cam2 )
{

	CRelativePoseBase* pRP = new CEstimatePose5Point(); 
	pRP->EstimatePose(lPts, rPts, cam1, cam2);
	delete pRP;

	return 0;
}

DLL_EXPORT Point3DDouble dll_TriangulatePt(Point2DDouble p, Point2DDouble q, 
	double *R0, double *t0, 
	double *R1, double *t1, double *error)
{
	return  TriangulatePt(p, q, R0, t0, R1, t1, error);
}



DLL_EXPORT void dll_Triangulate(vector<Point2DDouble> pts, vector<CameraPara> cams, Point3DDouble& gps,
	bool explicit_camera_centers,double& ferror)
{

	CTriangulateBase* pTriangulate = new CTriangulateCV();
	pTriangulate->Triangulate(pts, cams, gps, explicit_camera_centers, ferror);
	delete pTriangulate;
}


void dll_GenerateRMatrix(double omiga, double phi, double kapa, double* R)
{
	GenerateRMatrix(omiga, phi, kapa, R);
}

void dll_GrdToImg(double gx, double gy, double gz, double* ix, double* iy, 
	double* R, double* Ts, double f, double x0, double y0,
	int ht, int wd)
{
	GrdToImg(gx, gy, gz, ix, iy, R, Ts, f, x0, y0, ht, wd);
}

