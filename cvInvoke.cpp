
#include"stdio.h"

#include"cvInvoke.hpp"
#include"sift.hpp"
#include"dataBase.hpp"
#include "register.hpp"


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

//interface for memory
int DetectFileFeaturePts(char** filenames, int nFile, vector<ImgFeature>& imgFeatures)
{

	printf("[DetectFileFeaturePts] ... \n");

	CFeatureBase* pFeatDetect = new CSIFTFloat();
	//CPointFeatureBase* pFeatureData = new CSiftFeatureDataBinary();

	for(int i=0; i<nFile; i++)
	{
		printf("image: %s \n", filenames[i]);
		
		ImgFeature feats;
		pFeatDetect->Detect(filenames[i], feats);
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

	CMatchBase* pMatch = new CKNNMatch();

	for(int i=0; i<nImageNum; i++)
		for(int j=i+1; j<nImageNum; j++)
		{
			PairMatchRes mr;
			mr.lId = i;
			mr.rId = j;

			pMatch->Match(imgFeatures[i], imgFeatures[j], mr);

			printf("%d-%d  %d %lf \n", i, j, mr.matchs.size(), mr.inlierRatio );
			
			matchRes.push_back(mr);
		}

	delete pMatch;

	return 0;
}

