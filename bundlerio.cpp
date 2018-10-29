
#include "defines.hpp"
#include "bundlerio.hpp"
#include "badata.hpp"
#include "funcs.hpp"

#include "Matrix.h"

#include"baselib.h"

//#include "commondata.h"

//matrix
#include"vector.h"


//#include <iostream>
//#include <string>
//using namespace std;

int WritePMVSCamFile(char* file, CameraPara cam)
{
	//generate the projection matrix for PMVS and save it	  

	double focalLen = cam.focalLen;
	double outHt = cam.rows;
	double outWd = cam.cols;

	printf("ht: %lf wd: %lf \n", outHt, outWd);

	double K[9] =
	{ -focalLen, 0.0, 0.5 * outWd - 0.5,
	0.0, focalLen, 0.5 * outHt - 0.5,
	0.0, 0.0, 1.0 };

	double R[9];
	double T[3];
	memcpy(R, cam.R, sizeof(double) * 9);
	memcpy(T, cam.T, sizeof(double) * 3);

	//for R(x-T) to RX+T
	double it[3];
	/*it[0] = T[0];
	it[1] = T[1];
	it[2] = T[2];*/	
	mult(R, T, it, 3, 3, 1);
	for (int m = 0; m < 3; m++)
	{
		it[m] = -it[m];
	}
	
	//generate P
	double Ptmp[12] =
	{ R[0], R[1], R[2], it[0],
	  R[3], R[4], R[5], it[1],
	  R[6], R[7], R[8], it[2] };

	double P[12];
	dll_matrix_product(3, 3, 3, 4, K, Ptmp, P);
	dll_matrix_scale(3, 4, P, -1.0, P);

	FILE* f = fopen(file, "w");
	fprintf(f, "CONTOUR\n");
	fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[0], P[1], P[2], P[3]);
	fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[4], P[5], P[6], P[7]);
	fprintf(f, "%0.6f %0.6f %0.6f %0.6f\n", P[8], P[9], P[10], P[11]);
	//fprintf(f, "%lf %lf %lf \n", focalLen, outHt, outWd);
	////append the R and T for other applications
	//fprintf(f, "%lf %lf %lf \n", T[0], T[1], T[2]);
	//fprintf(f, "%lf %lf %lf \n", R[0], R[1], R[2]);
	//fprintf(f, "%lf %lf %lf \n", R[3], R[4], R[5]);
	//fprintf(f, "%lf %lf %lf \n", R[6], R[7], R[8]);
	fclose(f);


	return 0;
}

int SaveTracks(char* filepath, vector<TrackInfo>& tracks)
{
	FILE* fp = fopen(filepath, "w");

	for (int i = 0; i < tracks.size(); i++)
	{
		int nview = tracks[i].GetImageKeySum();
		Point3DDouble p3 = tracks[i].GetGround();
		fprintf(fp, "%lf %lf %lf, ", p3.p[0], p3.p[1], p3.p[2]);
		for (int k = 0; k < nview; k++)
		{
			ImageKey ik = tracks[i].GetImageKey(k);
			fprintf(fp, " %d-%d ", ik.first, ik.second);
		}
		fprintf(fp, "\n");
	}
	fclose(fp);

	return 0;
}

int SaveTracksToPly(char* filepath, vector<TrackInfo>& trackSeq, const vector<CameraPara>& cameras){

	//save the ba results: camera position, track points
	vector<Point3DDouble> goodGrds;
	vector<Point3DDouble> colors;
	printf("output the optimized track points.... \n");
	for (int i = 0; i<trackSeq.size(); i++)
	{
		if (trackSeq[i].GetImageKeySum() == 0)
			continue;

		if( trackSeq[i].derror<32 )
		{
			goodGrds.push_back(trackSeq[i].grd);

			Point3DDouble ptColor;
			ptColor.p[0] = 255;
			ptColor.p[1] = 0;
			ptColor.p[2] = 0;
			colors.push_back(ptColor);

		}
	}

	//save camera positions
	for (int i = 0; i<cameras.size(); i++)
	{
		int id = i;

		if (!cameras[i].bIsAddedIntoNet)
			continue;

		Point3DDouble cp;
		cp.p[0] = cameras[id].xs;
		cp.p[1] = cameras[id].ys;
		cp.p[2] = cameras[id].zs;

		Point3DDouble camColor;
		camColor.p[0] = 0;
		camColor.p[1] = 255;
		camColor.p[2] = 0;
		colors.push_back(camColor);

		goodGrds.push_back(cp);
	}

	WritePMVSPly(filepath, goodGrds, colors);

	return 0;
}

int SaveTracksToPly(char* filepath, vector<TrackInfo>& trackSeq,
					vector<int> cameraIDOrder, const vector<CameraPara>& cameras)
{	

	////calculate the average distance of image pair
	//double sumdis = 0;
	//int nsum = 0;
	//for(int i=0; i<cameraIDOrder.size()-1; i++)
	//{
	//	int id = cameraIDOrder[i];
	//	Point3DDouble cp1;
	//	cp1.p[0] = cameras[id].T[0];
	//	cp1.p[1] = cameras[id].T[1];
	//	cp1.p[2] = cameras[id].T[2];
	//	
	//	id = cameraIDOrder[i+1];
	//	Point3DDouble cp2;
	//	cp2.p[0] = cameras[id].T[0];
	//	cp2.p[1] = cameras[id].T[1];
	//	cp2.p[2] = cameras[id].T[2];

	//	double dis = distanceVec(cp1, cp2);

	//	sumdis += dis;
	//	nsum ++;
	//}
	//sumdis /= (double)(nsum);

	//double disThreshold = sumdis*20;
	//printf("distance threshold: %lf \n", disThreshold);

	//save the ba results: camera position, track points
	vector<Point3DDouble> goodGrds;
	vector<Point3DDouble> colors;
	printf("output the optimized track points.... \n");
	for(int i=0; i<trackSeq.size(); i++)
	{
		//if( trackSeq[i].valid<1 )
		//	continue;

		if( trackSeq[i].GetImageKeySum() == 0 )
			continue;

		//double gx = trackSeq[i].grd.p[0];
		//double gy = trackSeq[i].grd.p[1];
		//double gz = trackSeq[i].grd.p[2];
		////calculate the minimal distance between the 3D point and camera center
		//double minDis = 100000;
		//for(int k=0; k<cameraIDOrder.size(); k++)
		//{
		//	int id = cameraIDOrder[k];

		//	Point3DDouble cp;
		//	cp.p[0] = cameras[id].T[0];
		//	cp.p[1] = cameras[id].T[1];
		//	cp.p[2] = cameras[id].T[2];

		//	double dis = distanceVec(cp,  trackSeq[i].grd);
		//	if(minDis>dis)
		//		minDis = dis;
		//}

		////remove the 3D point far away from the cameras
		//if(minDis>disThreshold)
		//	continue;

		//if( trackSeq[i].derror<32 )
		{
			goodGrds.push_back( trackSeq[i].grd );

			Point3DDouble ptColor;
			ptColor.p[0] = 255;
			ptColor.p[1] = 0;
			ptColor.p[2] = 0;
			colors.push_back(ptColor);

			/*for(int j=0; j<trackSeq[i].views.size(); j++)
			{
				printf("%d %d ", trackSeq[i].views[j].first, trackSeq[i].views[j].second);
			}
			printf("\n");*/
		}
	}
	for(int i=0; i<cameraIDOrder.size(); i++)
	{
		int id = cameraIDOrder[i];

		Point3DDouble cp;
		cp.p[0] = cameras[id].T[0];
		cp.p[1] = cameras[id].T[1];
		cp.p[2] = cameras[id].T[2];

		Point3DDouble camColor;
		camColor.p[0] = 0;
		camColor.p[1] = 255;
		camColor.p[2] = 0;
		colors.push_back(camColor);

		goodGrds.push_back( cp );
	}

	WritePMVSPly(filepath, goodGrds, colors);
	
	return 0;
}




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
			Point2DDouble pt;				
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



int WriteBundlerOutFile(char* filepath, vector<CameraPara>& camParas )
{
	FILE *f = fopen(filepath, "w");
	if (f == NULL) 
	{
		printf("Error opening file %s for writing\n", filepath);
		return -1;
	}

	int num_images = camParas.size();
	int num_visible_points = 0; 

	/* Print version number */
	fprintf(f, "# Bundle file v0.3\n");
	/* Print number of cameras and points */
	fprintf(f, "%d %d\n", num_images, num_visible_points);

	for(int i=0; i<num_images; i++)
	{
		fprintf(f, "%0.10e %0.10e %0.10e \n",	camParas[i].focalLen, camParas[i].k1, camParas[i].k2);

		fprintf(f, "%0.10e %0.10e %0.10e\n", camParas[i].R[0], camParas[i].R[1],camParas[i].R[2]);
		fprintf(f, "%0.10e %0.10e %0.10e\n", camParas[i].R[3], camParas[i].R[4],camParas[i].R[5]);
		fprintf(f, "%0.10e %0.10e %0.10e\n", camParas[i].R[6], camParas[i].R[7],camParas[i].R[8]);

		//double t[3];
		//matrix_product(3, 3, 3, 1, camParas[i].R, camParas[i].t, t);
		//matrix_scale(3, 1, t, -1.0, t);
		//fprintf(f, "%0.10e %0.10e %0.10e\n", t[0], t[1], t[2]);
		fprintf(f, "%0.10e %0.10e %0.10e\n", camParas[i].T[0], camParas[i].T[1], camParas[i].T[2]);
	}
	fclose(f);

	return 0;
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

int ReadPMVSPly(char* filename, vector<stTrack>& tracks)
{
	char sline[256];

	FILE* fp = fopen(filename, "r");

	//skip the first two lines
	for(int i=0; i<3; i++)
		fgets(sline, 256, fp);

	//reading the number of points
	char tc[256];
	int  nPt;
	//fscanf(fp, "%s %s %d", tc, tc, &nPt);
	sscanf(sline, "%s %s %d", tc, tc, &nPt);

	//skip the following 10 lines
	for(int i=0; i<10; i++)
		fgets(sline, 256, fp);

	//reading the track points
	tracks.resize(nPt);
	for(int i=0; i<nPt; i++)
	{
		fscanf(fp, "%lf %lf %lf  %lf %lf %lf  %d %d %d ",
			&(tracks[i].x), &(tracks[i].y), &(tracks[i].z), 
			&(tracks[i].nx), &(tracks[i].ny), &(tracks[i].nz), 
			&(tracks[i].r), &(tracks[i].g), &(tracks[i].b) );
	}
	
	fclose(fp);
	
	return 0;
}

int WritePMVSPly(char* filename, vector<stTrack>& tracks)
{
	FILE* fp = fopen(filename, "w");
	fprintf(fp, "ply \n");
	fprintf(fp, "format ascii 1.0 \n");
	fprintf(fp, "element vertex %d \n", tracks.size());
	fprintf(fp, "property float x \n");
	fprintf(fp, "property float y \n");
	fprintf(fp, "property float z \n");
	fprintf(fp, "property float nx \n");
	fprintf(fp, "property float ny \n");
	fprintf(fp, "property float nz \n");
	fprintf(fp, "property uchar diffuse_red \n");
	fprintf(fp, "property uchar diffuse_green \n");
	fprintf(fp, "property uchar diffuse_blue \n");
	fprintf(fp, "end_header \n");

	for(int i=0; i<tracks.size(); i++)
	{
		fprintf(fp, "%lf %lf %lf  %lf %lf %lf %d %d %d \n",
			tracks[i].x,  tracks[i].y, tracks[i].z,
			tracks[i].nx, tracks[i].ny, tracks[i].nz,
			tracks[i].r,  tracks[i].g, tracks[i].b );
	}

	fclose(fp);

	return 0;
}


DLL_EXPORT int WritePMVSPly(char* filename, 
	double* px, double* py,  double* pz, int nPt)
{

	FILE* fp = fopen(filename, "w");
	fprintf(fp, "ply \n");
	fprintf(fp, "format ascii 1.0 \n");
	fprintf(fp, "element vertex %d \n", nPt);
	fprintf(fp, "property float x \n");
	fprintf(fp, "property float y \n");
	fprintf(fp, "property float z \n");
	fprintf(fp, "property float nx \n");
	fprintf(fp, "property float ny \n");
	fprintf(fp, "property float nz \n");
	fprintf(fp, "property uchar diffuse_red \n");
	fprintf(fp, "property uchar diffuse_green \n");
	fprintf(fp, "property uchar diffuse_blue \n");
	fprintf(fp, "end_header \n");

	for(int i=0; i<nPt; i++)
	{
		fprintf(fp, "%lf %lf %lf  %lf %lf %lf %d %d %d \n",
			px[i],  py[i], pz[i],
			1,1,1,255,0,0);
	}

	fclose(fp);


	return 0;
}

DLL_EXPORT int WritePMVSPly(char* filename, const vector<Point3DDouble>& gpts)
{

	FILE* fp = fopen(filename, "w");
	fprintf(fp, "ply \n");
	fprintf(fp, "format ascii 1.0 \n");
	fprintf(fp, "element vertex %d \n", gpts.size());
	fprintf(fp, "property float x \n");
	fprintf(fp, "property float y \n");
	fprintf(fp, "property float z \n");
	fprintf(fp, "property float nx \n");
	fprintf(fp, "property float ny \n");
	fprintf(fp, "property float nz \n");
	fprintf(fp, "property uchar diffuse_red \n");
	fprintf(fp, "property uchar diffuse_green \n");
	fprintf(fp, "property uchar diffuse_blue \n");
	fprintf(fp, "end_header \n");

	for(int i=0; i<gpts.size(); i++)
	{
		fprintf(fp, "%lf %lf %lf  %lf %lf %lf %d %d %d \n",
			gpts[i].p[0],  gpts[i].p[1], gpts[i].p[2],
			1.0,1.0,1.0,255,0,0);
	}

	fclose(fp);

	return 0;
}

DLL_EXPORT int WritePMVSPly(char* filename, const vector<Point3DDouble>& gpts, 
	const vector<Point3DDouble>& colors)
{

	FILE* fp = fopen(filename, "w");

	if (fp == NULL){
		printf("failed to open %s ! \n", filename);
		return -1;
	}

	fprintf(fp, "ply \n");
	fprintf(fp, "format ascii 1.0 \n");
	fprintf(fp, "element vertex %d \n", gpts.size());
	fprintf(fp, "property float x \n");
	fprintf(fp, "property float y \n");
	fprintf(fp, "property float z \n");
	fprintf(fp, "property float nx \n");
	fprintf(fp, "property float ny \n");
	fprintf(fp, "property float nz \n");
	fprintf(fp, "property uchar diffuse_red \n");
	fprintf(fp, "property uchar diffuse_green \n");
	fprintf(fp, "property uchar diffuse_blue \n");
	fprintf(fp, "end_header \n");

	for(int i=0; i<gpts.size(); i++)
	{
		fprintf(fp, "%lf %lf %lf  %lf %lf %lf %d %d %d \n",
			gpts[i].p[0],  gpts[i].p[1], gpts[i].p[2],
			1.0,1.0,1.0, int(colors[i].p[0]), int(colors[i].p[1]) , int(colors[i].p[2]));
	}

	fclose(fp);


	return 0;
}


int ReadPMVSPly(char* filename, stTrack** tracks, int* nTrack)
{
	char sline[256];

	FILE* fp = fopen(filename, "r");

	//skip the first two lines
	for(int i=0; i<3; i++)
		fgets(sline, 256, fp);

	//reading the number of points
	char tc[256];
	int  nPt;
	//fscanf(fp, "%s %s %d", tc, tc, &nPt);
	sscanf(sline, "%s %s %d", tc, tc, &nPt);

	//skip the following 10 lines
	for(int i=0; i<10; i++)
		fgets(sline, 256, fp);

	//reading the track points
	//tracks.resize(nPt);
	*tracks = (stTrack*)malloc(nPt*sizeof(stTrack));
	for(int i=0; i<nPt; i++)
	{
		fscanf(fp, "%lf %lf %lf  %lf %lf %lf  %d %d %d ",
			&((*tracks)[i].x), &((*tracks)[i].y), &((*tracks)[i].z), 
			&((*tracks)[i].nx), &((*tracks)[i].ny), &((*tracks)[i].nz), 
			&((*tracks)[i].r), &((*tracks)[i].g), &((*tracks)[i].b) );
	}

	fclose(fp);

	*nTrack = nPt;

	return 0;



	return 0;
}

int WritePMVSPly(char* filename, stTrack* tracks, int nTrack)
{

	FILE* fp = fopen(filename, "w");
	fprintf(fp, "ply \n");
	fprintf(fp, "format ascii 1.0 \n");
	fprintf(fp, "element vertex %d \n", nTrack);
	fprintf(fp, "property float x \n");
	fprintf(fp, "property float y \n");
	fprintf(fp, "property float z \n");
	fprintf(fp, "property float nx \n");
	fprintf(fp, "property float ny \n");
	fprintf(fp, "property float nz \n");
	fprintf(fp, "property uchar diffuse_red \n");
	fprintf(fp, "property uchar diffuse_green \n");
	fprintf(fp, "property uchar diffuse_blue \n");
	fprintf(fp, "end_header \n");

	for(int i=0; i<nTrack; i++)
	{
		fprintf(fp, "%lf %lf %lf  %lf %lf %lf %d %d %d \n",
			tracks[i].x,  tracks[i].y, tracks[i].z,
			tracks[i].nx, tracks[i].ny, tracks[i].nz,
			tracks[i].r,  tracks[i].g, tracks[i].b );
	}

	fclose(fp);

	return 0;
}
