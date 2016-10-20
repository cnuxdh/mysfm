
#ifndef CV_INVOKE_H
#define CV_INVOKE_H

#include <vector>
using namespace std;


#include"sift.hpp"
#include"dataBase.hpp"
#include"register.hpp"



int DetectFileFeaturePts(char** filenames, int nFile, char* outpath);

int DetectFileFeaturePts(char** filenames, int nFile, vector<ImgFeature>& imgFeatures);

int MatchImageFiles(vector<ImgFeature>& imgFeatures, vector<PairMatchRes>& matchRes);



#endif





