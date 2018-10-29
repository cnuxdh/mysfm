
#include "absOri.hpp"
#include "baselib.h"



//matrixlib
//#include "matrix/matrix.h"


//coredll
//#include "OrthoImage.h"



//corelib
#include "commondata.h"
#include "commonfile.h"
#include "CommonFuncs.h"
#include "LatLong-UTMconversion.h"


//generate the coordinate for plane projection according to the known Z axis
int GenerateProjectAxis(POINT3D za, POINT3D& xa, POINT3D& ya)
{
	
	POINT3D vz;
	vz.x = 0;
	vz.y = 0;
	vz.z = 1;
	
	POINT3D va;
	va.x = za.x;
	va.y = za.y;
	va.z = 0;
	
	xa = OuterProduct(vz, va);
	//xa.x = -xa.x;
	//xa.y = -xa.y;
	//xa.z = -xa.z;

	ya = OuterProduct(za, xa);
	
	return 0;
}


//generate the rotation matrix between two coordinates
int RotationAlign(vector<POINT3D> srcPts, vector<POINT3D> dstPts, double* rotationMatrix)
{	
	//calculate the R
	double sumC[9];
	memset(sumC, 0, sizeof(double)*9);
	for(int j=0; j<srcPts.size(); j++)
	{
		double gp[3]; //ground point
		double fp[3]; //free point

		fp[0] = srcPts[j].x;
		fp[1] = srcPts[j].y;
		fp[2] = srcPts[j].z;
		gp[0] = dstPts[j].x;
		gp[1] = dstPts[j].y;
		gp[2] = dstPts[j].z;

		double mc[9];
		dll_matrix_product(3, 1, 1, 3, gp, fp, mc);
		for(int k=0; k<9; k++)
			sumC[k] += mc[k];    
	}
	//svd
	double U[9];
	double S[3];
	double VT[9];
	dll_dgesvd_driver(3, 3, sumC, U, S, VT);
	dll_matrix_product(3, 3, 3, 3, U, VT, rotationMatrix );
	
	return 0;
}


/*
using only three points to calculate the absolute orientation, 
written by Donghai,Xie, 2015.10.10
the function can be invoked in RANSAC 
inputs:
	freePts,grdPts: free and ground control points
*/
int AbsOriP3( vector<POINT3D> freePts, vector<POINT3D> grdPts, 
			  stAbsPOS& absPosParams)
{
	if( freePts.size() != grdPts.size() )
	{
		return -1;
	}
	if( freePts.size() != 3)
	{
		printf("size must be 3 ! \n");
		return -1;
	}
	
	//calculate the centroid
	double fcx=0,fcy=0,fcz=0;
	double gcx=0,gcy=0,gcz=0;
	for(int j=0; j<freePts.size(); j++)
	{
		double xj = freePts[j].x;
		double yj = freePts[j].y;
		double zj = freePts[j].z;
		double gxj = grdPts[j].x;
		double gyj = grdPts[j].y;
		double gzj = grdPts[j].z;

		fcx += xj;		fcy += yj;		fcz += zj;
		gcx += gxj;		gcy += gyj;		gcz += gzj;
	}
    POINT3D freeCenterPt,grdCenterPt;	
	freeCenterPt.x = fcx / double(freePts.size());
	freeCenterPt.y = fcy / double(freePts.size());
	freeCenterPt.z = fcz / double(freePts.size());
	grdCenterPt.x  = gcx / double(freePts.size());
	grdCenterPt.y  = gcy / double(freePts.size());
	grdCenterPt.z  = gcz / double(freePts.size());

	//
	for(int j=0; j<freePts.size(); j++)
	{
		freePts[j].x = freePts[j].x-freeCenterPt.x;
		freePts[j].y = freePts[j].y-freeCenterPt.y;
		freePts[j].z = freePts[j].z-freeCenterPt.z;
		grdPts[j].x = grdPts[j].x-grdCenterPt.x;
		grdPts[j].y = grdPts[j].y-grdCenterPt.y;
		grdPts[j].z = grdPts[j].z-grdCenterPt.z;
	}

	//calculate the scale
	double sumScale = 0;
	int    ns = 0;
	for(int j=0; j<freePts.size(); j++)
	{
		double fx = freePts[j].x;
		double fy = freePts[j].y;
		double fz = freePts[j].z;
		double gx = grdPts[j].x;
		double gy = grdPts[j].y;
		double gz = grdPts[j].z;
		double gLen = sqrt(  gx*gx + gy*gy + gz*gz );
		double fLen = sqrt(  fx*fx + fy*fy + fz*fz );

		if(fLen!=0)
		{
			sumScale += gLen / fLen;
			ns ++;
		}
	}
	absPosParams.scale = sumScale / (double)(ns);

	//calculate the R
	double sumC[9];
	memset(sumC, 0, sizeof(double)*9);
	for(int j=0; j<freePts.size(); j++)
	{
		double gp[3]; //ground point
		double fp[3]; //free point

		fp[0] = freePts[j].x*absPosParams.scale;
		fp[1] = freePts[j].y*absPosParams.scale;
		fp[2] = freePts[j].z*absPosParams.scale;

		gp[0] = grdPts[j].x;
		gp[1] = grdPts[j].y;
		gp[2] = grdPts[j].z;

		double mc[9];
		//dll_matrix_product(3, 1, 1, 3, gp, fp, mc);
		dll_matrix_product(3, 1, 1, 3, fp, gp, mc); //revised by xiedonghai, 2018.5.27
		for(int k=0; k<9; k++)
			sumC[k] += mc[k];    
	}
	//svd
	double U[9];
	double S[3];
	double VT[9];
	double mS[9];
	dll_dgesvd_driver(3, 3, sumC, U, S, VT);

	double Ut[9];
	double V[9];
	dll_matrix_transpose(3, 3, U, Ut);
	dll_matrix_transpose(3, 3, VT, V);
	double m2[9];
	dll_matrix_product(3, 3, 3, 3, V, Ut, m2 );

	//adding constraints , https://en.wikipedia.org/wiki/Kabsch_algorithm
	double det = dll_matrix_determint3(m2);
	memset(mS, 0, 9 * sizeof(double));
	mS[0] = 1;  mS[4] = 1; mS[8] = det;
	double m1[9];
	//dll_matrix_product(3, 3, 3, 3, U, mS, m1);
	//dll_matrix_product(3, 3, 3, 3, m1, VT, absPosParams.R);
	dll_matrix_product(3, 3, 3, 3, V, mS, m1);
	dll_matrix_product(3, 3, 3, 3, m1, Ut, absPosParams.R);



	/*
	//printf("determinant of R: %lf \n", det);
	if (det < 0){
		for (int i = 0; i < 9; i++)
			absPosParams.R[i] *= -1;
	}
    det = dll_matrix_determint3(absPosParams.R);
	printf("determinant of R: %lf \n", det);
	*/

	//calculate the T
	double gpc[3]; //center of the ground points
	double fpc[3]; //center of the free points
	gpc[0]=grdCenterPt.x;	gpc[1]=grdCenterPt.y;	gpc[2]=grdCenterPt.z;
	fpc[0]=freeCenterPt.x;	fpc[1]=freeCenterPt.y;	fpc[2]=freeCenterPt.z;	
	double tc[3];
	dll_matrix_product(3, 3, 3, 1, absPosParams.R, fpc, tc);
	for(int i=0; i<3; i++)
		absPosParams.T[i] = gpc[i] - absPosParams.scale*tc[i];
	
	return 0;
}

int GetZoneNumber1(vector<CameraPara> camParas)
{
	double slon = 0;
	int n = 0;
	for (int i = 0; i<camParas.size(); i++)
	{
		if (camParas[i].focalLen != 0)
		{
			slon += camParas[i].lon;
			n++;
		}
	}

	slon /= (double)(n);
	int zoneNumber = int((slon + 180) / 6) + 1;

	return zoneNumber;
}

int GetZoneNumber1(vector<stPOS> camParas)
{
	double slon = 0;
	int n = 0;
	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].f != 0)
		{
			slon += camParas[i].lon;
			n ++;
		}
	}

	slon /= (double)(n);
	int zoneNumber = int((slon + 180)/6) + 1;

	return zoneNumber;
}

double AbsOriOrthogonal(stAbsPOS& absPosParams, vector<stPOS>& camParas, vector<stTrack>& tracks)
{

	return 0;
}


double AbsOriOrthogonal(stAbsPOS& absPosParams, vector<CameraPara>& camParas)
{
	//collect valid cameras
	vector<int> validCameraIndex;
	for (int i = 0; i<camParas.size(); i++)
	{
		if (camParas[i].focalLen != 0)
		{
			validCameraIndex.push_back(i);
		}
	}

	if (validCameraIndex.size() < 3)
	{
		printf("[AbsOriOrthogonal]: the number of camera is less than 3 ! \n");
		return -1;
	}


	// form R(X-T) to RX+T
	for (int i = 0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];

		camParas[ci].xs = camParas[ci].T[0];
		camParas[ci].ys = camParas[ci].T[1];
		camParas[ci].zs = camParas[ci].T[2];

		//dll_matrix_invert(3, camParas[ci].R, R1);
		double t1[3];
		dll_matrix_product(3, 3, 3, 1, camParas[ci].R, camParas[ci].T, t1);

		camParas[ci].T[0] = -t1[0];
		camParas[ci].T[1] = -t1[1];
		camParas[ci].T[2] = -t1[2];
	}

	//convert from lon,lat to ground coordinate
	int zoneNumber = GetZoneNumber1(camParas);
	printf("zonenumber: %d \n", zoneNumber);
	for (int i = 0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];
		/*double lat = camParas[ci].lat;
		double lon = camParas[ci].lon;
		double gx, gy;
		LLtoUTM(23, lat, lon, gy, gx, zoneNumber);
		camParas[ci].gx = gx;
		camParas[ci].gy = gy;
		camParas[ci].gz = camParas[ci].altitude;*/
		printf("%lf %lf %lf \n", camParas[ci].gx, camParas[ci].gy, camParas[ci].gz);
	}

	//RANSAC absolute orientation
	int ransac_rounds = 1000;
	int minWrongNumber = 10000000;
	int DistanceThrshold = 16;
	for (int round = 0; round < ransac_rounds; round++)
	{
		int indices[3];
		ChooseRandIndex(validCameraIndex.size(), 3, indices);

		vector<POINT3D> freePts;
		vector<POINT3D> grdPts;
		freePts.resize(3);
		grdPts.resize(3);
		for (int i = 0; i<3; i++)
		{
			int ri = validCameraIndex[indices[i]];

			freePts[i].x = camParas[ri].xs;
			freePts[i].y = camParas[ri].ys;
			freePts[i].z = camParas[ri].zs;

			grdPts[i].x = camParas[ri].gx;
			grdPts[i].y = camParas[ri].gy;
			grdPts[i].z = camParas[ri].gz;
		}

		stAbsPOS absPara;
		AbsOriP3(freePts, grdPts, absPara);

		//calculate the errors
		int nWrongNumber = 0;
		for (int i = 0; i<validCameraIndex.size(); i++)
		{
			int ci = validCameraIndex[i];

			double fP[3];
			fP[0] = camParas[ci].xs;
			fP[1] = camParas[ci].ys;
			fP[2] = camParas[ci].zs;

			double tp[3];
			dll_matrix_product(3, 3, 3, 1, absPara.R, fP, tp);
			for (int k = 0; k<3; k++)
				fP[k] = absPara.scale*tp[k] + absPara.T[k];

			double gP[3];
			gP[0] = camParas[ci].gx;
			gP[1] = camParas[ci].gy;
			gP[2] = camParas[ci].gz;

			double len = 0;
			for (int k = 0; k<3; k++)
				len += (fP[k] - gP[k])*(fP[k] - gP[k]);
			len = sqrt(len);
			if (len>DistanceThrshold)
				nWrongNumber++;
		}

		if (minWrongNumber>nWrongNumber)
		{
			minWrongNumber = nWrongNumber;
			absPosParams = absPara;
		}
	}
	printf("minimux wrong number: %d \n", minWrongNumber);


	//transform each camera parameters from current coordinate to new coordinate
	double Rg[9];
	dll_matrix_invert(3, absPosParams.R, Rg);
	double Tg[3];
	memcpy(Tg, absPosParams.T, 3 * sizeof(double));

	double sumErr = 0;
	vector<double> vecError;
	vecError.resize(validCameraIndex.size());
	for (int i = 0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];

		//new R
		double newR[9];
		dll_matrix_product(3, 3, 3, 3, camParas[ci].R, Rg, newR);

		//new T
		double newT[3];
		dll_matrix_product(3, 3, 3, 1, newR, Tg, newT);
		for (int k = 0; k<3; k++)
			newT[k] = absPosParams.scale*camParas[ci].T[k] - newT[k];

		memcpy(camParas[ci].R, newR, sizeof(double)* 9);
		memcpy(camParas[ci].T, newT, sizeof(double)* 3);

		//rotation angle
		camParas[ci].pitch = atan(camParas[ci].R[5] / camParas[ci].R[8]) / PI * 180;
		camParas[ci].roll = asin(-camParas[ci].R[2]) / PI * 180;
		camParas[ci].yaw = atan(camParas[ci].R[1] / camParas[ci].R[0]) / PI * 180;
		//printf("rotation angle: %lf %lf %lf \n", camParas[ci].pitch, camParas[ci].roll, camParas[ci].yaw);

		//convert form RX+T to R( X - (-inv(R)*T) )
		double t1[3];
		double R1[9];
		dll_matrix_invert(3, camParas[ci].R, R1);
		dll_matrix_product(3, 3, 3, 1, R1, camParas[ci].T, t1);

		camParas[ci].xs = -t1[0];
		camParas[ci].ys = -t1[1];
		camParas[ci].zs = -t1[2];

		//the error
		double distance = sqrt((camParas[ci].xs - camParas[ci].gx)*(camParas[ci].xs - camParas[ci].gx) +
			(camParas[ci].ys - camParas[ci].gy)*(camParas[ci].ys - camParas[ci].gy) +
			(camParas[ci].zs - camParas[ci].gz)*(camParas[ci].zs - camParas[ci].gz));

		vecError[i] = distance;
		sumErr += distance;
		printf("abs distance: %lf \n", distance);
	}

	//printf("\n");
	sumErr /= validCameraIndex.size();

	printf("Absolute orientation error.... \n");
	for (int i = 0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];
		printf("error: %lf  rotation angle: %lf %lf %lf \n", vecError[i],
			camParas[ci].pitch, camParas[ci].roll, camParas[ci].yaw);

		//remove the camera with large rotation angle 
		//if( fabs(camParas[ci].pitch)>10 || fabs(camParas[ci].roll)>10 )
		//	camParas[ci].f = 0;
	}
	printf("\n\n");

	return sumErr;
}

double AbsOriOrthogonal(stAbsPOS& absPosParams, vector<CameraPara>& camParas, vector<TrackInfo>& tracks)
{
	//collect valid cameras
	vector<int> validCameraIndex;
	for(int i=0; i<camParas.size(); i++)
	{
		if(camParas[i].bIsAddedIntoNet)
		{
			validCameraIndex.push_back(i);
		}
	}

	if( validCameraIndex.size() < 3 )
	{
		printf("[AbsOriOrthogonal]: the number of camera is less than 3 ! \n");
		return -1;
	}


	// form R(X-T) to RX+T
	printf("input for absolute orientation ..... \n");
	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];		

		camParas[ci].xs = camParas[ci].T[0];
		camParas[ci].ys = camParas[ci].T[1];
		camParas[ci].zs = camParas[ci].T[2];

		//dll_matrix_invert(3, camParas[ci].R, R1);
		double t1[3];
		dll_matrix_product(3, 3, 3, 1, camParas[ci].R, camParas[ci].T, t1);
		camParas[ci].T[0] = -t1[0];
		camParas[ci].T[1] = -t1[1];
		camParas[ci].T[2] = -t1[2];

		printf("%lf %lf %lf - %lf %lf %lf \n", 
			camParas[ci].xs, camParas[ci].ys, camParas[ci].zs,
			camParas[ci].gx, camParas[ci].gy, camParas[ci].gz);
	}

	/*
	//convert from lon,lat to ground coordinate
	int zoneNumber =GetZoneNumber1(camParas);
	printf("zonenumber: %d \n", zoneNumber);
	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];
		double lat = camParas[ci].lat;
		double lon = camParas[ci].lon;
		double gx,gy;
		LLtoUTM(23, lat, lon, gy, gx, zoneNumber);
		camParas[ci].gx = gx;
		camParas[ci].gy = gy;
		camParas[ci].gz = camParas[ci].altitude;
		printf("%lf %lf %lf \n", gx, gy, camParas[ci].altitude);
	}
	*/

	//RANSAC absolute orientation
	int ransac_rounds  = 1000;
	int minWrongNumber = 10000000;
	int DistanceThrshold = 16;
	for(int round = 0; round < ransac_rounds; round++)
	{
		int indices[3];
		ChooseRandIndex(validCameraIndex.size(), 3, indices);
		
		vector<POINT3D> freePts;
		vector<POINT3D> grdPts;
		freePts.resize(3);
		grdPts.resize(3);
		for(int i=0; i<3; i++)
		{
			int ri = validCameraIndex[ indices[i] ];

			freePts[i].x = camParas[ri].xs;
			freePts[i].y = camParas[ri].ys;
			freePts[i].z = camParas[ri].zs;

			grdPts[i].x = camParas[ri].gx;
			grdPts[i].y = camParas[ri].gy;
			grdPts[i].z = camParas[ri].gz;
		}

		stAbsPOS absPara;
		AbsOriP3(freePts, grdPts, absPara);
		

		//calculate the errors
		int nWrongNumber     = 0; 
        for(int i=0; i<validCameraIndex.size(); i++)
		{
			int ci = validCameraIndex[i];

			double fP[3];
			fP[0] = camParas[ci].xs;
			fP[1] = camParas[ci].ys;
			fP[2] = camParas[ci].zs;
			
			double tp[3];
			dll_matrix_product(3, 3, 3, 1, absPara.R, fP, tp);
			for(int k=0; k<3; k++)
				fP[k] = absPara.scale*tp[k] + absPara.T[k];

			double gP[3];
			gP[0] = camParas[ci].gx;
			gP[1] = camParas[ci].gy;
			gP[2] = camParas[ci].gz;

			double len = 0;
			for(int k=0; k<3; k++)
				len += (fP[k]-gP[k])*(fP[k]-gP[k]);
			len = sqrt(len);
			if(len>DistanceThrshold)
				nWrongNumber++;
		}
		
		if( minWrongNumber>nWrongNumber )
		{
			minWrongNumber = nWrongNumber;
			absPosParams = absPara;
		}
	}	

	printf("minimux wrong number: %d \n", minWrongNumber);

	printf("absolute orientation results: \n");
	printf("scale: %lf \n", absPosParams.scale);
	printf(" rotation.... \n");
	for (int j = 0; j < 3; j++){
		for (int i = 0; i < 3; i++){
			printf("%lf ", absPosParams.R[j * 3 + i]);
		}
		printf("\n");
	}
	
	

	printf("\n translation.... \n");
	for (int j = 0; j < 3; j++){
		printf("%lf ", absPosParams.T[j]);
	}
	printf("\n\n");
	
	

	//transform each camera parameters from current coordinate to new coordinate
	double Rg[9];
	dll_matrix_invert(3, absPosParams.R, Rg);
	double Tg[3];
	memcpy(Tg, absPosParams.T, 3*sizeof(double) );

	double sumErr = 0;
	vector<double> vecError;
	vecError.resize(validCameraIndex.size());
	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];

		//new R
		double newR[9];
		dll_matrix_product(3, 3, 3, 3, camParas[ci].R, Rg, newR);

		//new T
		double newT[3];
		dll_matrix_product(3, 3, 3, 1, newR, Tg, newT);
		for(int k=0; k<3; k++)
			newT[k] = absPosParams.scale*camParas[ci].T[k] - newT[k];

		memcpy( camParas[ci].R, newR, sizeof(double)*9 );
		memcpy( camParas[ci].T, newT, sizeof(double)*3 );

		//rotation angle
		camParas[ci].pitch = atan(  camParas[ci].R[5]/camParas[ci].R[8] )/PI*180; 
		camParas[ci].roll  = asin( -camParas[ci].R[2] )/PI*180;
		camParas[ci].yaw   = atan(  camParas[ci].R[1]/camParas[ci].R[0])/PI*180;
		//printf("rotation angle: %lf %lf %lf \n", camParas[ci].pitch, camParas[ci].roll, camParas[ci].yaw);

		//convert form RX+T to R( X - (-inv(R)*T) )
		double t1[3];
		double R1[9];
		dll_matrix_invert(3, camParas[ci].R, R1);
		dll_matrix_product(3, 3, 3, 1, R1, camParas[ci].T, t1);

		//camParas[ci].T[0] = -t1[0];
		//camParas[ci].T[1] = -t1[1];
		//camParas[ci].T[2] = -t1[2];
		camParas[ci].xs = -t1[0]; 
		camParas[ci].ys = -t1[1]; 
		camParas[ci].zs = -t1[2];

		//the error
		double distance = sqrt( (camParas[ci].xs-camParas[ci].gx)*(camParas[ci].xs-camParas[ci].gx) +
								(camParas[ci].ys-camParas[ci].gy)*(camParas[ci].ys-camParas[ci].gy) + 
								(camParas[ci].zs-camParas[ci].gz)*(camParas[ci].zs-camParas[ci].gz) );

		vecError[i] = distance;
		sumErr += distance;
		
		//printf("absolute orientation distance: %lf \n", distance);
	}

	//printf("\n");
	sumErr /= validCameraIndex.size();

	printf("Absolute orientation error.... \n");
	for(int i=0; i<validCameraIndex.size(); i++)
	{
		int ci = validCameraIndex[i];
		printf("error: %lf  rotation angle: %lf %lf %lf \n", vecError[i],
			camParas[ci].pitch, camParas[ci].roll, camParas[ci].yaw);

		//remove the camera with large rotation angle 
		//if( fabs(camParas[ci].pitch)>10 || fabs(camParas[ci].roll)>10 )
		//	camParas[ci].f = 0;
	}
	printf("\n\n");

	//transform all tracks 
	for(int i=0; i<tracks.size(); i++)
	{
		double fP[3];
		fP[0] = tracks[i].grd.p[0];
		fP[1] = tracks[i].grd.p[1];
		fP[2] = tracks[i].grd.p[2];

		double tp[3];
		dll_matrix_product(3, 3, 3, 1, absPosParams.R, fP, tp);
		for(int k=0; k<3; k++)
			fP[k] = absPosParams.scale*tp[k] + absPosParams.T[k];

		tracks[i].grd.p[0] = fP[0];
		tracks[i].grd.p[1] = fP[1];
		tracks[i].grd.p[2] = fP[2];
	}

	return sumErr;
}