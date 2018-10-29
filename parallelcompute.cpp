

#include"parallelcompute.hpp"



//feature detection class
CDetectFeatPts::CDetectFeatPts()
{

}

CDetectFeatPts::~CDetectFeatPts()
{

}

int CDetectFeatPts::run(char** filenames, int nfile, double& dProgress)
{
	double dStep = 100.0 / (double)nfile;

	for (int i = 0; i<nfile; i++)
	{
		CImageDataBase* pImageData = new CImageFeature();
		pImageData->Load(filenames[i]);
		pImageData->DetectPtFeature(SIFT_FLOAT_FEATURE);
		vecImageDataPointer.push_back(pImageData);

		dProgress += dStep;
	}
	
	return 0;
}




//batch class
int CGenerateSparsePts::run(char** filenames, int nfile, double& dProgress)
{



	return 0;
}