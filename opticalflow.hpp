
#ifndef OPTICALFLOW_HEADER
#define OPTICALFLOW_HEADER

#include <vector>
using namespace std;




_declspec(dllexport) int SimpleFlow(char* lfile, char* rfile, char* outfile, 
									int& ht, int& wd, 
									vector<double>& tx,  
									vector<double>& ty);


//void HarrisDetect(IplImage* img, vector<int>& featx, vector<int>& featy, int number);

void HarrisDetect(char* filename, vector<int>& featx, vector<int>& featy, int number);



#endif