
#include"stdio.h"





#include"cvInvoke.hpp"



#include"sift.hpp"




//multiple cpu cores




//core function for feature detection of multiple images
int DetectFileFeaturePts(char** filenames, int nFile, char* outpath)
{
	
	//retrive the file title
	
	CFeatureBase* pFeatDetect = new CSIFTFloat();
	
	for(int i=0; i<nFile; i++)
	{
		printf("image: %s \n", filenames[i]);
		
		ImgFeature feats;
		pFeatDetect->Detect(filenames[i], feats);
		
		//save the feat into the file
		
		
	}
	
	delete pFeatDetect;
	
	
	return 0;
}

