

#ifndef PARALLEL_COMPUTING
#define PARALLEL_COMPUTING

#include"imagedata.hpp"


//feature detection
class DLL_EXPORT CDetectFeatPts
{
public:
	CDetectFeatPts();
	~CDetectFeatPts();

	int run(char** filenames, int nfile, double& dProgress);

private:
	//save the feature 
	vector<CImageDataBase*> vecImageDataPointer;
};


//match



//bundle adjust


//batch
class CGenerateSparsePts
{
public:
	CGenerateSparsePts();
	~CGenerateSparsePts();
	
	int run(char** filenames, int nfile, double& dProgress);

private:
		
};















#endif
