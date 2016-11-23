
#include "warp.hpp"

#include "CalcAngle.h"

#include "Matrix.h"


//generate the homography matrix: H = R - T.N'/d, x1=H.x2
//R,T:  rotation and translation, R(X - T)
//N',d: N'=(nx,ny,nz) is the plane normal, plane is described by: nx.X+ny.Y+nz.Z + distance=0
int GenerateHomographyMatrix(double omiga, double kapa, double phi, 
							 double tx, double ty, double tz,
							 double nx, double ny, double nz, 
							 double distance, double focal,
							 double* H )
{
	
	double R[9];
	double T[3];
	double pn[3];
	double K[9];

	for(int i=0; i<9; i++)
		K[i] = 0;
	K[0] = K[4] = 1;
	K[8] = -1.0/focal;

	GenerateRMatrixDirect(omiga, kapa, phi, R);

	T[0] = tx;
	T[1] = ty;
	T[2] = tz;

	double newT[3];
	mult(R, T, newT, 3, 3, 1);
	newT[0] = -newT[0];
	newT[1] = -newT[1];
	newT[2] = -newT[2];

	pn[0] = nx;
	pn[1] = ny;
	pn[2] = nz;
	
	double tn[9];

	mult(newT, pn, tn, 3, 1, 3);

	for(int i=0; i<9; i++)
	{
		H[i] = R[i] - tn[i]/distance;
	}

	double inversK[9];
	for(int i=0; i<9; i++)
		inversK[i] = K[i];
	invers_matrix(inversK, 3);
	
	double m1[9];
	mult(H, inversK, m1, 3, 3, 3);
	mult(K, m1, H, 3, 3, 3);

	invers_matrix(H, 3);

	return 0;
}



int HomographyWarp(unsigned char* pSrc, int sht, int swd, 
				   unsigned char* pDst, int& dht, int& dwd, 
				   double* H)
{
	for(int j=0; j<dht; j++)
		for(int i=0; i<dwd; i++)
		{
			double dx = double(i);
			double dy = double(j);

			//normalize the coordinates
			dx = dx - dwd*0.5;
			dy = dht*0.5 - dy;
						
			//homography transform
			double denominator = dx*H[6] + dy*H[7] + H[8];

			double sx = (dx*H[0] + dy*H[1] + H[2]) / denominator;
			double sy = (dx*H[3] + dy*H[4] + H[5]) / denominator;

			int ix = (int)(sx); 
			int iy = (int)(sy);

			ix += swd*0.5;
			iy = sht*0.5 - iy;
			
			if(ix>=0 && ix<swd && iy>=0 && iy<sht)
			{
				pDst[j*dwd+i] = pSrc[iy*swd+ix];
			}
		}

	return 0;
}
