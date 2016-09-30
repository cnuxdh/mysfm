
#include"stdio.h"





#include"cvInvoke.hpp"
#include"sift.hpp"
#include"dataBase.hpp"


//corelib
#include"commonfile.h"


//multiple cpu cores




//core function for feature detection of multiple images
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

