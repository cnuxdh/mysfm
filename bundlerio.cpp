

#include "bundlerio.hpp"

//matrix
#include"vector.h"




/* reading bundler *.out file, including camera pose information, ground point and its projections
*/
int ReadBundlerOutFile(char* filename, vector<stPOS>& camParas, vector<stTrack>& tracks )
{
	char strLine[256];
	int numCamera, numPt;
	double f,k1,k2;
	double R[9];
	double T[3];
	int i,j;

	FILE* fp =fopen(filename, "r");
	fgets(strLine,256, fp);
	printf("%s\n", strLine);
	fscanf(fp, "%d %d", &numCamera, &numPt);
	printf("%d %d\n", numCamera, numPt);

	stPOS cam;
	for(i=0; i<numCamera; i++)
	{
		fscanf(fp, "%lf %lf %lf", &cam.f, &cam.k1, &cam.k2);
		for(j=0; j<9; j++)	
			fscanf(fp, "%lf", &(cam.R[j]) );
		for(j=0; j<3; j++)	
			fscanf(fp, "%lf", &(cam.T[j]) );
	
		camParas.push_back(cam);
	}

	//read 3d point
	double rgb[3];
	int nPt;
	int nImageIndex, nPtIndex;
	double ix,iy;
	
	for( i=0; i<numPt; i++)
	{
		stTrack singleTrack;

		fscanf(fp, "%lf %lf %lf", &singleTrack.x, &singleTrack.y, &singleTrack.z );
		fscanf(fp, "%d %d %d", &singleTrack.r, &singleTrack.g, &singleTrack.b );
		fscanf(fp, "%d", &nPt);
		
		for(int k=0; k<nPt; k++)
		{
			POINT2 pt;				
			fscanf(fp,"%d %d %lf %lf", &nImageIndex, &nPtIndex, &pt.x, &pt.y);
			singleTrack.imgpt.push_back(pt);
			singleTrack.imgid.push_back(nImageIndex);
			singleTrack.ptid.push_back(nPtIndex);
		}		
		tracks.push_back(singleTrack);
	}
	fclose(fp);

	return 1;
}






//////////////////////////////////////////////////////////////////////////
//3D model 
//ply mode file header
/*
static char ply_header[] = 
"ply\n"
"format ascii 1.0\n"
"element face 0\n"
"property list uchar int vertex_indices\n"
"element vertex %d\n"
"property float x\n"
"property float y\n"
"property float z\n"
"property uchar diffuse_red\n"
"property uchar diffuse_green\n"
"property uchar diffuse_blue\n"
"end_header\n";
*/

static char ply_header[] = 
"ply\n"
"format ascii 1.0\n"
"element face 0\n"
"property list uchar int vertex_indices\n"
"element vertex %d\n"
"property float x\n"
"property float y\n"
"property float z\n"
"property uchar diffuse_red\n"
"property uchar diffuse_green\n"
"property uchar diffuse_blue\n"
"end_header\n";


CPlyModel::CPlyModel()
{

}

CPlyModel::~CPlyModel()
{

}

int CPlyModel::Save(char *modelFile, vector<Point3DDouble> pts)
{
	
	int r=255,g=0,b=0;

	int num_good_pts = pts.size();

	FILE* fp = fopen(modelFile, "w");

	fprintf(fp, ply_header, num_good_pts);
    for(int i=0; i<num_good_pts; i++)
	{
		fprintf(fp, "%0.6e %0.6e %0.6e %d %d %d\n", 
			Vx(pts[i]), Vy(pts[i]), Vz(pts[i]),
			r,g,b);
	}
	fclose(fp);


	return 1;
}

int CPlyModel::Save(char* modelFile, vector<Point3DDouble> pts, vector<Point3DDouble> colors)
{

	int r,g,b;
	int num_good_pts = pts.size();

	FILE* fp = fopen(modelFile, "w");

	fprintf(fp, ply_header, num_good_pts);
	for(int i=0; i<num_good_pts; i++)
	{
		r = Vx(colors[i]);
		g = Vy(colors[i]);
		b = Vz(colors[i]);

		//fprintf(fp, "%0.6e %0.6e %0.6e %d %d %d\n", 
		fprintf(fp, "%6.3lf %6.3lf %6.3lf %d %d %d\n", 
			Vx(pts[i]), Vy(pts[i]), Vz(pts[i]),
			r,g,b);
	}
	fclose(fp);

	return 1;
}