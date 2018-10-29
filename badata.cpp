
#include "badata.hpp"


ImgFeature::ImgFeature()
{
	featPts.clear();
}

ImgFeature::~ImgFeature()
{

}

Point2DDouble ImgFeature::GetCenteredPt(int ptIndex)
{
	Point2DDouble cp;

	cp.p[0] = featPts[ptIndex].cx;
	cp.p[1] = featPts[ptIndex].cy;

	return cp;
}
Point2DDouble ImgFeature::GetTopLeftPt(int ptIndex)
{
	Point2DDouble cp;

	cp.p[0] = featPts[ptIndex].x;
	cp.p[1] = featPts[ptIndex].y;

	return cp;
}

int ImgFeature::SaveImgFeaturePts(string filepath)
{
	char *cfile = (char*)filepath.c_str();

	FILE* fp = fopen(cfile, "wb");
	if (fp == NULL)
		return 0;

	//dim of image
	fwrite(&ht, sizeof(int), 1, fp);
	fwrite(&wd, sizeof(int), 1, fp);

	//number of feature points
	int numFeat = featPts.size();
	fwrite(&numFeat, sizeof(int), 1, fp);
	
	//the dim of feature 
	int ndim = 128;
	fwrite(&ndim, sizeof(int), 1, fp);
	//write feature vector one by one
	for (int i = 0; i<numFeat; i++)
	{
		//y x scale orientation
		float cy, cx, scale, ori;
		cx = featPts[i].cx;
		cy = featPts[i].cy;
		scale = featPts[i].scl;
		ori   = featPts[i].ori;

		fwrite(&cx, sizeof(float), 1, fp);
		fwrite(&cy, sizeof(float), 1, fp);
		fwrite(&scale, sizeof(float), 1, fp);
		fwrite(&ori, sizeof(float), 1, fp);

		//feature vector
		float feat[128];
		for (int k = 0; k < 128; k++)
			feat[k] = featPts[i].feat[k];
		fwrite(feat, sizeof(float), 128, fp);
	}
	fclose(fp);

	return 0;
}

int ImgFeature::ReadImgFeaturePts(string filepath)
{
	int num = 0;

	char *cfile = (char*)filepath.c_str();

	int dim = 128; //dimension by feature detector

	FILE* fp = fopen(cfile, "rb");
	fread(&ht, sizeof(int), 1, fp);
	fread(&wd, sizeof(int), 1, fp);
	fread(&num, sizeof(int), 1, fp); //feature point number
	fread(&dim, sizeof(int), 1, fp); //sift feature length
	float* pbuffer = new float[num*(dim+4)];
	fread(pbuffer, sizeof(float), num*(dim + 4), fp);
	fclose(fp);


	featPts.resize(num);
	for (int i = 0; i < num; i++)
	{
		featPts[i].cx  = pbuffer[i*(dim + 4)];
		featPts[i].cy  = pbuffer[i*(dim + 4)+1];
		featPts[i].scl = pbuffer[i*(dim + 4)+2];
		featPts[i].ori = pbuffer[i*(dim + 4)+3];

		featPts[i].extra = -1;
		featPts[i].id = i;

		featPts[i].feat.resize(dim);
		for (int ki = 0; ki < dim; ki++)
			featPts[i].feat[ki] = pbuffer[i*(dim + 4) + 4 + ki];
	}
	delete pbuffer;

	//
	//int resizeDim = 128;
	//for (int i = 0; i <num; i++)
	//{
	//	float cx, cy, scl, ori;
	//	float pd[128];

	//	fread(&cy,  sizeof(float), 1, fp);
	//	fread(&cx,  sizeof(float), 1, fp);
	//	fread(&scl, sizeof(float), 1, fp);
	//	fread(&ori, sizeof(float), 1, fp);
	//	fread(pd,   sizeof(float), dim, fp);

	//	featPts[i].cx = cx;
	//	featPts[i].cy = cy;
	//	featPts[i].scl = scl;
	//	featPts[i].ori = ori;
	//	
	//	featPts[i].extra = -1;
	//	featPts[i].id = i;

	//	////resize the input feature vector, written by xiedonghai, 2015.10.16
	//	//int nstep = dim / resizeDim;
	//	////float sp[128];
	//	//for (int k = 0; k<resizeDim; k++)
	//	//{
	//	//	double sumValue = 0;
	//	//	for (int n = k*nstep; n<(k + 1)*nstep; n++)
	//	//	{
	//	//		sumValue += pd[n];
	//	//	}
	//	//	float sp = sumValue / nstep;
	//	//	featPts[i].feat.push_back(sp);
	//	//}

	//	featPts[i].feat.resize(dim);
	//	for (int ki = 0; ki < dim; ki++)
	//		featPts[i].feat[ki] = pd[ki];

	//}

	//fclose(fp);

	return 0;
}





TrackInfo::TrackInfo()
{
	views.clear();
	ctrl = 0;
}

TrackInfo::~TrackInfo()
{

}

int GetGoodTracks(vector<TrackInfo> srcTracks, vector<TrackInfo>& goodTracks)
{
	goodTracks.clear();
	for (int i = 0; i < srcTracks.size(); i++)
	{
		if (srcTracks[i].GetImageKeySum() == 0)
			continue;

		if (srcTracks[i].derror < 32)
		{
			goodTracks.push_back(srcTracks[i]);
		}
	}

	
	//analyze the height
	vector<double> zbuffer;
	for (int i = 0; i < goodTracks.size(); i++)
	{
		zbuffer.push_back(goodTracks[i].GetGround().p[2]);
	}
	std::sort(zbuffer.begin(), zbuffer.end());
	int half = zbuffer.size()*0.5;
	int crop = zbuffer.size()*0.1;
	double meanz = zbuffer[half];
	double maxz = zbuffer[zbuffer.size() - crop];
	double minz = zbuffer[crop];
	double zthresh = (maxz - minz) * 2;
	
	vector<TrackInfo> tracks = goodTracks;
	goodTracks.clear();
	for (int i = 0; i < tracks.size(); i++)
	{
		double z = tracks[i].GetGround().p[2];
		if (fabs(z - meanz) < zthresh)
			goodTracks.push_back(tracks[i]);
	}

	return 0;
}
