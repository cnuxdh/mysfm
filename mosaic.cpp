
#include "mosaic.hpp"

//#include "bundlerio.hpp"


//#include "opencv2/core/core_c.h"
#include "opencv2/core/core.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/imgproc/imgproc_c.h"
//#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/video/tracking.hpp"

#include "opencv2/stitching/detail/seam_finders.hpp"
#include "opencv2/stitching/detail/blenders.hpp"
#include "opencv2/stitching/detail/exposure_compensate.hpp"
#include "opencv2/stitching/detail/util.hpp"
#include "opencv2/stitching/warpers.hpp"


//gdal
#include "gdal_priv.h"
#include "ogr_spatialref.h"


//coredll
#include "geotiff.h"
//#include "MosaicFuncs.h"
//#include "CommonFuncs.h"


//corelib
/*
#include "corelib/imagefunc.h"
#include "Corelib/image.h"
#include "Corelib/commonfile.h"
*/
#include"commondata.h"


using namespace cv;
using namespace cv::detail;
using namespace std;


void CalculateVarianceAndMean(int* hist, int nHist, double& mean, double& std)
{
	double sum = 0;
	for(int i=1; i<nHist; i++)
	{
		sum += hist[i];
	}

	mean = 0;
	for(int i=1; i<nHist; i++)
	{
		double ratio = (double)(hist[i])/(double)(sum);
		mean += ratio*(double)(i);
	}

	std = 0;
	for(int i=1; i<nHist; i++)
	{
		double ratio = (double)(hist[i])/(double)(sum-1);
		std += ( (double)(i)-mean )*((double)(i)-mean ) * ratio;
	}
	std = sqrt(std);
}


void AutoLevelRange(int* hist, int nHist, double ratio, int& nMin, int& nMax)
{
	int sum = 0;
	int delSum = 0;
	int lg,rg;
	int s = 0;

	for(int i=0; i<nHist; i++)
		sum += hist[i];
	
	delSum = sum*ratio;
	lg = 0;
	s = 0;
	while(s<delSum)
	{
		s += hist[lg];
		lg ++;
	}

	s = 0;
	rg = 255;  //considering the white cloud
	while(s<delSum)
	{
		s += hist[rg];
		rg --;
	}

	nMax = rg;
	nMin = lg;
}


/* function: auto level
   default: 1%
*/
void   AutoLevel(unsigned char* image, int ht, int wd)
{
	int i,j;
	int index;
	int hist[256];
	int imgSize;
	float delRatio = 0.02;
	int   delSum;
	int   s;
	int   lg,rg;
	float scale;

	memset(hist, 0, sizeof(int)*256);

	imgSize = ht*wd;
	int sum = 0;
	//1: 统计灰度分布的直方图
	for(i=0; i<imgSize; i++)
	{
		if(image[i]>0 )
		{
			hist[ image[i] ] ++;
			sum++;
		}
	}

	//2: 搜索最黑与最白的边界
	delSum = sum*delRatio;
	lg = 0;
	s = 0;
	while(s<delSum)
	{
		s += hist[lg];
		lg ++;
	}

	s = 0;
	rg = 255;  //considering the white cloud
	while(s<delSum)
	{
		s += hist[rg];
		rg --;
	}

	//3: 根据搜索的结果来对图像进行拉伸
	scale = 1 / (float)(rg-lg);
	for(i=0; i<imgSize; i++)
	{
		if(image[i]<=lg)
			image[i] = 0;
		else if(image[i]>=rg)
			image[i] = 255;
		else
			image[i] = (image[i]-lg)*scale*255;
	}
}


//read image of type unsigned int (16bit) with one band
int GetUInt16GrayImage(char* filename, Mat& image, double ratio)
{
	GDALAllRegister();

	stGeoInfo geoinfo;
	GetGeoInformation(filename, geoinfo);

	GDALDataset* poDataset = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	int ht = geoinfo.ht;
	int wd = geoinfo.wd;	
	int zht = ht*ratio;
	int zwd = wd*ratio;
	unsigned short* pBuffer = (unsigned short*)malloc(zht*zwd*sizeof(unsigned short));
	GDALRasterBand* poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Read, 0, 0, wd, ht, pBuffer, zwd, zht, GDT_UInt16, 0, 0);
	GDALClose( (GDALDatasetH) poDataset );
	

	int nMin = 1000000;
	int nMax = 0;
	for(int i=0; i<zht*zwd; i++)
	{
		nMin = min(nMin, pBuffer[i]);
		nMax = max(nMax, pBuffer[i]);
	}
	
	printf("%d %d \n", nMin, nMax);

	//convert to BYTE type
	unsigned char* pDispBuffer = (unsigned char*)malloc(zht*zwd*sizeof(unsigned char));
	memset(pDispBuffer, 0, zht*zwd*sizeof(unsigned char));
	for(int i=0; i<zht*zwd; i++)
	{
		pDispBuffer[i] = (double)(pBuffer[i]-nMin) / (double)(nMax-nMin) * 255;
	}

	//auto level algorithm
	AutoLevel(pDispBuffer, zht, zwd);
	
	//convert to Mat format
	image.create(zht, zwd, CV_8UC3);
	image.setTo(Scalar::all(0));
	uchar* pImage = image.data;
	for(int j=0; j<zht; j++)
		for(int i=0; i<zwd; i++)
		{
			pImage[j*3*zwd+3*i]   = pDispBuffer[j*zwd+i];
			pImage[j*3*zwd+3*i+1] = pDispBuffer[j*zwd+i];
			pImage[j*3*zwd+3*i+2] = pDispBuffer[j*zwd+i];
		}

	free(pBuffer);
	free(pDispBuffer);
	return 1;
}


//read image of type unsigned int (16bit) with one band
int GetBYTEGrayImage(char* filename, Mat& image, double ratio)
{
	GDALAllRegister();

	stGeoInfo geoinfo;
	GetGeoInformation(filename, geoinfo);

	GDALDataset* poDataset = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	int ht = geoinfo.ht;
	int wd = geoinfo.wd;	
	int zht = ht*ratio;
	int zwd = wd*ratio;
	unsigned char* pBuffer = (unsigned char*)malloc(zht*zwd*sizeof(unsigned char));
	GDALRasterBand* poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Read, 0, 0, wd, ht, pBuffer, zwd, zht, GDT_Byte, 0, 0);
	GDALClose( (GDALDatasetH) poDataset );

	//convert to Mat format
	image.create(zht, zwd, CV_8UC3);
	image.setTo(Scalar::all(0));
	uchar* pImage = image.data;
	for(int j=0; j<zht; j++)
		for(int i=0; i<zwd; i++)
		{
			pImage[j*3*zwd+3*i]   = pBuffer[j*zwd+i];
			pImage[j*3*zwd+3*i+1] = pBuffer[j*zwd+i];
			pImage[j*3*zwd+3*i+2] = pBuffer[j*zwd+i];
		}

	free(pBuffer);
	return 1;
}

int GetImage(char* filename, Mat& image, double ratio)
{
	GDALAllRegister();

	stGeoInfo geoinfo;
	GetGeoInformation(filename, geoinfo);

	GDALDataset* poDataset = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	int ht = geoinfo.ht;
	int wd = geoinfo.wd;	
	int zht = ht*ratio;
	int zwd = wd*ratio;
	unsigned char* r = (unsigned char*)malloc(zht*zwd);
	unsigned char* g = (unsigned char*)malloc(zht*zwd);
	unsigned char* b = (unsigned char*)malloc(zht*zwd);

	GDALRasterBand* poBand = poDataset->GetRasterBand(3);
	poBand->RasterIO(GF_Read, 0, 0, wd, ht, r, zwd, zht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(2);
	poBand->RasterIO(GF_Read, 0, 0, wd, ht, g, zwd, zht, GDT_Byte, 0, 0);
	poBand = poDataset->GetRasterBand(1);
	poBand->RasterIO(GF_Read, 0, 0, wd, ht, b, zwd, zht, GDT_Byte, 0, 0);
	GDALClose( (GDALDatasetH) poDataset );
		
	//AutoLevel(r, zht, zwd);
	//AutoLevel(g, zht, zwd);
	//AutoLevel(b, zht, zwd);
	
	image.create(zht, zwd, CV_8UC3);
	image.setTo(Scalar::all(0));
	uchar* pBuffer = image.data;
	for(int j=0; j<zht; j++)
		for(int i=0; i<zwd; i++)
		{
			pBuffer[j*3*zwd+3*i]   = r[j*zwd+i];
			pBuffer[j*3*zwd+3*i+1] = g[j*zwd+i];
			pBuffer[j*3*zwd+3*i+2] = b[j*zwd+i];
		}

	free(r);
	free(g);
	free(b);

	return 1;
}


int GenerateMask(char* filename, Mat& mask, double ratio)
{
	stGeoInfo geoinfo;
	GetGeoInformation(filename, geoinfo);

	GDALDataset* poDataset = (GDALDataset *) GDALOpen( filename, GA_ReadOnly );
	GDALRasterBand* poBand = poDataset->GetRasterBand(1);
	int ht = geoinfo.ht;
	int wd = geoinfo.wd;	

	int zht = ht*ratio;
	int zwd = wd*ratio;
	unsigned char* pBench = (unsigned char*)malloc(zht*zwd);
	poBand->RasterIO(GF_Read, 0, 0, wd, ht, pBench, zwd, zht, GDT_Byte, 0, 0);
	GDALClose( (GDALDatasetH) poDataset );


	mask.create(zht, zwd, CV_8U);
	mask.setTo(Scalar::all(0));
	uchar* pBuffer = mask.data;
	for(int j=0; j<zht; j++)
		for(int i=0; i<zwd; i++)
		{
			//pBuffer[j*zwd+i] = pBench[j*zwd+i];
			if( pBench[j*zwd+i]>0 )
				pBuffer[j*zwd+i] = 255;
		}
		free(pBench);

		int dilation_type = MORPH_ELLIPSE;
		int dilation_size = 5;
		Mat element = getStructuringElement( dilation_type,
			Size( 2*dilation_size + 1, 2*dilation_size+1 ),
			Point( dilation_size, dilation_size ) );
		dilate(mask, mask, element);
		erode(mask, mask, element);

		return 1;
}

//fill the black part of the image based on interpolation
void ImageFill(Mat& img)
{


}




CVoroniSeamFinder::CVoroniSeamFinder()
{

}

CVoroniSeamFinder::~CVoroniSeamFinder()
{

}

/*
 mask: image mask to show the forground
 corner1,corner2: the position of left-top at the same coordinate 
*/
int CVoroniSeamFinder::FindPair(unsigned char* mask1, int ht1, int wd1, 
								unsigned char* mask2, int ht2, int wd2, 
								iPoint corner1, iPoint corner2)
{
	//generate masks
	Mat m1,m2;
	m1.create(ht1, wd1, CV_8UC1);
	m2.create(ht2, wd2, CV_8UC1);

	uchar* p1 = m1.data;
	uchar* p2 = m2.data;

    for(int j=0; j<ht1; j++)
		for(int i=0; i<wd1; i++)
			p1[j*wd1+i] = mask1[j*wd1+i];

	for(int j=0; j<ht2; j++)
		for(int i=0; i<wd2; i++)
			p2[j*wd2+i] = mask2[j*wd2+i];
	
	//data conversion
	vector<Mat> masks;
	masks.push_back(m1);
	masks.push_back(m2);
	
	vector<Point> corners;
	Point pt1,pt2;
	pt1.x = corner1.x;
	pt1.y = corner1.y;
	pt2.x = corner2.x;
	pt2.y = corner2.y;
	corners.push_back(pt1);
	corners.push_back(pt2);

	//the size of image is only needed 
	vector<Mat> srcs;
	Mat s1,s2;
	s1.create(ht1, wd1, CV_8UC1);
	s2.create(ht2, wd2, CV_8UC1);
	srcs.push_back(s1);
	srcs.push_back(s2);
    
	//seam finder 
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	//Ptr<SeamFinder> seam_finder = new detail::GraphCutSeamFinder();
	seam_finder->find(srcs, corners, masks);

	//reset the mask
	p1 = masks[0].data;
	p2 = masks[1].data;
	for(int j=0; j<ht1; j++)
		for(int i=0; i<wd1; i++)
			mask1[j*wd1+i] = p1[j*wd1+i];
	for(int j=0; j<ht2; j++)
		for(int i=0; i<wd2; i++)
			mask2[j*wd2+i] = p2[j*wd2+i];

	return 1;
}

int CVoroniSeamFinder::FindPair(MatrixByte mask1, MatrixByte mask2, iPoint corner1, iPoint corner2)
{
	int ht1 = mask1.ht;
	int wd1 = mask1.wd;
	int ht2 = mask2.ht;
	int wd2 = mask2.wd;

	unsigned char* pMask1 = mask1.pbuffer;
	unsigned char* pMask2 = mask2.pbuffer;

	//generate masks
	Mat m1,m2;
	m1.create(ht1, wd1, CV_8UC1);
	m2.create(ht2, wd2, CV_8UC1);

	uchar* p1 = m1.data;
	uchar* p2 = m2.data;

	for(int j=0; j<ht1; j++)
		for(int i=0; i<wd1; i++)
			p1[j*wd1+i] = pMask1[j*wd1+i];

	for(int j=0; j<ht2; j++)
		for(int i=0; i<wd2; i++)
			p2[j*wd2+i] = pMask2[j*wd2+i];

	//data conversion
	vector<Mat> masks;
	masks.push_back(m1);
	masks.push_back(m2);

	vector<Point> corners;
	Point pt1,pt2;
	pt1.x = corner1.x;
	pt1.y = corner1.y;
	pt2.x = corner2.x;
	pt2.y = corner2.y;
	corners.push_back(pt1);
	corners.push_back(pt2);

	//the size of image is only needed 
	vector<Mat> srcs;
	Mat s1,s2;
	s1.create(ht1, wd1, CV_8UC1);
	s2.create(ht2, wd2, CV_8UC1);
	srcs.push_back(s1);
	srcs.push_back(s2);

	//seam finder 
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	seam_finder->find(srcs, corners, masks);

	//reset the mask
	p1 = masks[0].data;
	p2 = masks[1].data;
	for(int j=0; j<ht1; j++)
		for(int i=0; i<wd1; i++)
			pMask1[j*wd1+i] = p1[j*wd1+i];
	for(int j=0; j<ht2; j++)
		for(int i=0; i<wd2; i++)
			pMask2[j*wd2+i] = p2[j*wd2+i];

	return 1;
}

CMultiBandBlender::CMultiBandBlender()
{

}

CMultiBandBlender::~CMultiBandBlender()
{

}

int CMultiBandBlender::Prepair(vector<MatrixByte>& images, vector<iPoint> corners)
{


	return 1;
}

int CMultiBandBlender::Feed( MatrixByte image, MatrixByte mask, iPoint corner)
{

	return 1;
}

int CMultiBandBlender::Blend(vector<MatrixByte> images, vector<MatrixByte> masks, vector<iPoint> corners, 
							 MatrixByte& dst, MatrixByte& dstMask)
{
	//data conversion
	//mask 
	vector<Mat> masks_warped;
	masks_warped.resize(masks.size());
	for(int i=0; i<masks.size(); i++)
	{
		int ht = masks[i].ht;
		int wd = masks[i].wd;

		masks_warped[i].create(ht, wd, CV_8U);
		uchar* pdata = masks_warped[i].data;

		for(int m=0; m<ht; m++)
			for(int n=0; n<wd; n++)
			{
				pdata[m*wd+n] = masks[i].pbuffer[m*wd+n];
			}
	}

	//image conversion
	vector<Mat> images_warped;
	images_warped.resize(images.size());
	for(int i=0; i<masks.size(); i++)
	{
		int ht = images[i].ht;
		int wd = images[i].wd;

		images_warped[i].create(ht, wd, CV_8U);
		uchar* pdata = images_warped[i].data;

		for(int m=0; m<ht; m++)
			for(int n=0; n<wd; n++)
			{
				pdata[m*wd+n] = images[i].pbuffer[m*wd+n];
			}
	}

	vector<Point> vps;
	vector<Size>  sizes;
	vps.resize(corners.size());
	for(int i=0; i<corners.size(); i++)
	{
		vps[i].x = corners[i].x;
		vps[i].y = corners[i].y;
	}
	sizes.resize(images.size());
	for(int i=0; i<corners.size(); i++)
	{
		sizes[i].width  = images[i].wd;
		sizes[i].height = images[i].ht;
	}

	//1. define the blender
	Ptr<Blender> blender;
	int blend_type = Blender::MULTI_BAND;
	bool try_gpu   = false;
	blender = Blender::createDefault(blend_type, try_gpu);
   
	//get the area of mosaic
	Size   dst_sz = resultRoi(vps, sizes).size();

	//calculate the area of blend
	double blend_strength = 5.0;
	float  blend_width = sqrt(static_cast<float>(dst_sz.area())) * blend_strength / 100.f;
	
	//prepare
	MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
	mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));	
	blender->prepare(vps, sizes);

	//2. add the images one by one	
	Mat dilated_mask, seam_mask, mask, mask_warped,img_warped_s;
	int num_images = images.size();
	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		int ht = images_warped[img_idx].rows;
		int wd = images_warped[img_idx].cols;
		mask_warped.create(ht, wd, CV_8U);

		dilate(masks_warped[img_idx], dilated_mask, Mat());
		resize(dilated_mask, seam_mask, mask_warped.size());
		mask_warped = seam_mask & mask_warped;
	   
		//imwrite("d:\\image.jpg", images_warped[img_idx]);
		//imwrite("d:\\mask.jpg", mask_warped);
		images_warped[img_idx].convertTo(img_warped_s, CV_16S);
		blender->feed(img_warped_s, mask_warped, vps[img_idx]);
		//blender->feed(images_warped[img_idx], mask_warped, vps[img_idx]);
	}
	
	Mat result, result_mask;
	blender->blend(result, result_mask);

	imwrite("d:\\blend.jpg", result);

	return 1;
}

/*
int CMultiBandBlender::GeoBlend( vector<Mat> images, char* outFilePath )
{
	printf("Begin blend ... \n");

	vector<Mat>   masks;
	vector<Point> pts;
	vector<Size>  sizes;

	images.resize(2);
	masks.resize(2);
	pts.resize(2);
	sizes.resize(2);

	printf("%d \n", images[0].type());

	//prepare images masks
	for (int i = 0; i < images.size(); ++i)
	{
		masks[i].create(images[i].size(), CV_8U); //must be CV_8U
		masks[i].setTo(Scalar::all(255));
	}

	pts[0].x = 0;	pts[0].y = 0;
	pts[1].x = 10;	pts[1].y = 60;

	for(int i=0; i<2; i++)
	{
		sizes[i].width = images[i].cols;
		sizes[i].height = images[i].rows;
	}

	//1. seam finder
	//seam finder 
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	seam_finder->find(images, pts, masks);

	//2. define the blender
	Ptr<Blender> blender;
	int blend_type = Blender::MULTI_BAND;
	bool try_gpu   = false;
	blender = Blender::createDefault(blend_type, try_gpu);

	//get the area of mosaic
	Size   dst_sz = resultRoi(pts, sizes).size();
	//calculate the area of blend
	double blend_strength = 60.0;
	float  blend_width = sqrt(static_cast<float>(dst_sz.area())) * blend_strength / 100.f;

	//prepare
	MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
	mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));	
	blender->prepare(pts, sizes);

	//2.1 add the images one by one	
	Mat dilated_mask, seam_mask, mask, mask_warped,img_warped_s;
	int num_images = images.size();
	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		int ht = images[img_idx].rows;
		int wd = images[img_idx].cols;
		mask_warped.create(images[img_idx].size(), CV_8U);
		mask_warped.setTo(Scalar::all(255));
		printf("%d \n", mask_warped.type());

		dilate(masks[img_idx], dilated_mask, Mat());
		resize(dilated_mask, seam_mask, mask_warped.size());
		mask_warped = seam_mask & mask_warped;

		//imwrite("d:\\image.jpg", images_warped[img_idx]);
		//imwrite("d:\\mask.jpg", mask_warped);
		images[img_idx].convertTo(img_warped_s, CV_16S);
		blender->feed(img_warped_s, mask_warped, pts[img_idx]);
		//blender->feed(images_warped[img_idx], mask_warped, vps[img_idx]);
	}

	//2.2 do blend
	Mat result, result_mask;
	blender->blend(result, result_mask);

	//output
	imwrite("d:\\blend.jpg", result);


	printf("Finished and Out is  'blend.jpg' ");


	return 1;
}
*/


/*
if the background pixel is near the image foreground, then fill it using the image pixel, written by xdh, 2014.12.31
*/
void SpreadImageBorder(Mat& img, Mat& dst)
{
	int ht = img.rows;
	int wd = img.cols;

	for(int j=0; j<ht; j++)
		for(int i=0; i<wd; i++)
		{
			int l = max(0, i-2);
			int r = min(wd-1, i+2);
			int t = max(0, j-2);
			int b = min(ht-1, j+2);			

			int num = 0;
			if(img.channels()==1)
			{
				if( img.at<uchar>(j,i)!=0 )
					continue;

				int value = 0;
				int sum = 0;
				for(int m=t; m<=b; m++)
					for(int n=l; n<=r; n++)
					{
						value = img.at<uchar>(m,n);
						if(value>0)
						{
							sum += value;
							num ++;
						}
					}
				if(num>0)
				{
					dst.at<uchar>(j,i) = sum / num;
				}
			}
			else if(img.channels()==3)
			{
				int vr,vg,vb;

				vr = img.at<cv::Vec3b>(j,i)[0];
				vg = img.at<cv::Vec3b>(j,i)[1];
				vb = img.at<cv::Vec3b>(j,i)[2];

				//for black pixel
				if( (vr+vg+vb)<36 )
				{
					int sr = 0;
					int sg = 0;
					int sb = 0;
					for(int m=t; m<=b; m++)
						for(int n=l; n<=r; n++)
						{
							vr = img.at<cv::Vec3b>(m,n)[0];
							vg = img.at<cv::Vec3b>(m,n)[1];
							vb = img.at<cv::Vec3b>(m,n)[2];

							if(vr>5 && vg>5 && vb>5)
							{
								num ++;
								sr += vr;
								sg += vg;
								sb += vb;
							}
						}
					if(num>0)
					{
						dst.at<cv::Vec3b>(j,i)[0] = sr / num;
						dst.at<cv::Vec3b>(j,i)[1] = sg / num;
						dst.at<cv::Vec3b>(j,i)[2] = sb / num;
					}
				}
			}
		}

}


int CMultiBandBlender::Blend(char** filenames, int nFile, double blendWeight, char* outFilePath)
{
	printf("Begin blend ... \n");

	char* imgfile1  = "F:\\data\\Dodging\\redApple.jpg";
	char* imgfile2  = "F:\\data\\Dodging\\greenApple.jpg";

	//char* imgfile1  = "F:\\data\\mosaic\\lena.jpg";
	//char* imgfile2  = "F:\\data\\mosaic\\lena-blue.jpg";
	//char* maskfile1 = "d:\\data\\mosaic\\mask.jpg";
	//char* maskfile2 = "d:\\data\\mosaic\\mask.jpg";
	//char* imgfile1  = "F:\\data\\mosaic\\2013.jpg";
	//char* imgfile2  = "F:\\data\\mosaic\\2014.jpg";
	//char* imgfile1  = "F:\\data\\Images\\lena.jpg";
	//char* imgfile2  = "F:\\data\\Images\\lena-change.jpg";

	//char* imgfile1  = "F:\\Data\\Dodging\\rgb\\2013-482-1131.jpg";
	//char* imgfile2  = "F:\\Data\\Dodging\\rgb\\2014-482-1132.jpg";

	//char* imgfile1  = "F:\\Data\\Dodging\\rgb\\2014-483-1130.jpg";
	//char* imgfile2  = "F:\\Data\\Dodging\\rgb\\2014-483-1131.jpg";
	

	vector<Mat>   images;
	vector<Mat>   masks;
	vector<Point> pts;
	vector<Size>  sizes;

	images.resize(2);
	masks.resize(2);
	pts.resize(2);
	sizes.resize(2);

	//read images
	//images[0] = imread(imgfile1); 
	//images[1] = imread(imgfile2); 
	//printf("%d \n", images[0].type());

	double ratio = 1;
	GetImage(imgfile1, images[0], ratio);
	GetImage(imgfile2, images[1], ratio);
	GenerateMask(imgfile1, masks[0], ratio);
	GenerateMask(imgfile2, masks[1], ratio);

	//printf(" step: %d  elemsize: %d  type: %d \n", images[0].step, images[0].elemSize(), images[0].type());
    

	//loading the filled image
	//GetImage("F:\\Data\\Images\\lena-fill.jpg", images[0], ratio);
    
	/*
	//enlarge the image border
	for(int i=0; i<2; i++)
	{
		Mat dst = images[i].clone();
		SpreadImageBorder(images[i], dst);
		images[i] = dst.clone();
	
		char file[256];
		sprintf(file, "d:\\image_%d.jpg", i);
		imwrite(file, images[i]);
	}
	*/

	int xoff = images[0].cols*ratio*0.5;
	int yoff = 0*ratio;
	pts[0].x = 0;	pts[0].y = 0;
	pts[1].x = xoff;	pts[1].y = yoff;
	for(int i=0; i<2; i++)
	{
		sizes[i].width = images[i].cols;
		sizes[i].height = images[i].rows;
	}

	//1. seam finder
	//seam finder 
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	seam_finder->find(images, pts, masks);

	//2. define the blender
	Ptr<Blender> blender;
	int blend_type = Blender::MULTI_BAND;
	bool try_gpu   = false;
	blender = Blender::createDefault(blend_type, try_gpu);

	//get the area of mosaic
	Size   dst_sz = resultRoi(pts, sizes).size();
	//calculate the area of blend
	double blend_strength = blendWeight+1;
	float  blend_width = sqrt(static_cast<float>(dst_sz.area())) * blend_strength / 100.f;

	//prepare
	MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
	mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));	
	blender->prepare(pts, sizes);

	//2.1 add the images one by one	
	Mat dilated_mask, seam_mask, mask, mask_warped,img_warped_s;
	int num_images = images.size();
	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		int ht = images[img_idx].rows;
		int wd = images[img_idx].cols;
		
		//generate mask
		mask_warped.create(images[img_idx].size(), CV_8U);
		mask_warped.setTo(Scalar::all(255));
		printf("%d \n", mask_warped.type());
		dilate(masks[img_idx], dilated_mask, Mat());
		resize(dilated_mask, seam_mask, mask_warped.size());
		mask_warped = seam_mask & mask_warped;
		
		//prepare image
		images[img_idx].convertTo(img_warped_s, CV_16S);

		char file[256];
		sprintf(file, "d:\\image_%d.jpg", img_idx);
		imwrite(file, img_warped_s);
		sprintf(file, "d:\\mask_%d.jpg", img_idx);
		imwrite(file, mask_warped);
		
		blender->feed(img_warped_s, mask_warped, pts[img_idx]);
		//blender->feed(images_warped[img_idx], mask_warped, vps[img_idx]);
	}

	//2.2 do blend
	Mat result, result_mask;
	blender->blend(result, result_mask);

	//output
	imwrite(outFilePath, result);


	printf("Finished and Out is  'blend.jpg' ");

	return 1;
}

int CMultiBandBlender::GeoBlend( vector<string> filenames, string outFilePath )
{

	printf("Begin blend ... \n");

	char* imgfile1  = "d:\\data\\mosaic\\lena2.jpg";
	char* imgfile2  = "d:\\data\\mosaic\\IMG_0105.jpg";
	//char* maskfile1 = "d:\\data\\mosaic\\mask.jpg";
	//char* maskfile2 = "d:\\data\\mosaic\\mask.jpg";

	vector<Mat>   images;
	vector<Mat>   masks;
	vector<Point> pts;
	vector<Size>  sizes;

	images.resize(2);
	masks.resize(2);
	pts.resize(2);
	sizes.resize(2);

	//read images
	images[0] = imread(imgfile1); 
	images[1] = imread(imgfile2); 

	printf("%d \n", images[0].type());

	//prepare images masks
	for (int i = 0; i < images.size(); ++i)
	{
		masks[i].create(images[i].size(), CV_8U); //must be CV_8U
		masks[i].setTo(Scalar::all(255));
	}

	pts[0].x = 0;	pts[0].y = 0;
	pts[1].x = 10;	pts[1].y = 60;

	for(int i=0; i<2; i++)
	{
		sizes[i].width = images[i].cols;
		sizes[i].height = images[i].rows;
	}

	//1. seam finder
	//seam finder 
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	seam_finder->find(images, pts, masks);

	//2. define the blender
	Ptr<Blender> blender;
	int blend_type = Blender::MULTI_BAND;
	bool try_gpu   = false;
	blender = Blender::createDefault(blend_type, try_gpu);
	
	//get the area of mosaic
	Size   dst_sz = resultRoi(pts, sizes).size();
	//calculate the area of blend
	double blend_strength = 60.0;
	float  blend_width = sqrt(static_cast<float>(dst_sz.area())) * blend_strength / 100.f;

	//prepare
	MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
	mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));	
	blender->prepare(pts, sizes);

	//2.1 add the images one by one	
	Mat dilated_mask, seam_mask, mask, mask_warped,img_warped_s;
	int num_images = images.size();
	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		int ht = images[img_idx].rows;
		int wd = images[img_idx].cols;
		mask_warped.create(images[img_idx].size(), CV_8U);
		mask_warped.setTo(Scalar::all(255));
		printf("%d \n", mask_warped.type());

		dilate(masks[img_idx], dilated_mask, Mat());
		resize(dilated_mask, seam_mask, mask_warped.size());
		mask_warped = seam_mask & mask_warped;

		//imwrite("d:\\image.jpg", images_warped[img_idx]);
		//imwrite("d:\\mask.jpg", mask_warped);
		images[img_idx].convertTo(img_warped_s, CV_16S);
		blender->feed(img_warped_s, mask_warped, pts[img_idx]);
		//blender->feed(images_warped[img_idx], mask_warped, vps[img_idx]);
	}

	//2.2 do blend
	Mat result, result_mask;
	blender->blend(result, result_mask);

	//output
	imwrite("d:\\blend.jpg", result);


	printf("Finished and Out is  'blend.jpg' ");

	return 1;
}

int CMultiBandBlender::GeoBlendSingle( char** filenames, int nFile, char* outFilePath )
{	
	vector<Mat> images;
	vector<Mat> masks;
	masks.resize(nFile);
	images.resize(nFile);
	vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);

	GDALAllRegister();

	double ratio = 1.0;
	//1. generate masks
	for(int i=0; i<nFile; i++)
	{
		GenerateMask(filenames[i], masks[i], ratio);
		//imwrite("d:\\smallMask.jpg", masks[i]);
		//GetUInt16GrayImage(filenames[i], images[i], ratio);				
		GetBYTEGrayImage( filenames[i], images[i], ratio);
		GetGeoInformation( filenames[i], geoArray[i]);
	}

	//2. find seam
	double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);

		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

    vector<Point> pts;
	pts.resize(nFile);
	for(int i=0; i<nFile; i++)
	{
		int x = (geoArray[i].left - minx) / maxrx * ratio;
		int y = (maxy - geoArray[i].top) / maxry * ratio;
		pts[i].x = x;
		pts[i].y = y;
	} 
	vector<Size> sizes;
	sizes.resize(images.size());
	for(int i=0; i<nFile; i++)
	{
		sizes[i].width  = images[i].cols;
		sizes[i].height = images[i].rows;
	}
	
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	seam_finder->find(images, pts, masks);
	
	/*
	for(int i=0; i<nFile; i++)
	{
		char file[256];
		sprintf(file, "d:\\mask_%d.jpg", i);
		imwrite(file, masks[i]);
	}
	*/
	

	//3. define the blender
	Ptr<Blender> blender;
	int blend_type = Blender::MULTI_BAND;
	bool try_gpu   = false;
	blender = Blender::createDefault(blend_type, try_gpu);
	
	//get the area of mosaic
	Size   dst_sz = resultRoi(pts, sizes).size();
	
	//calculate the area of blend
	double blend_strength = 5;
	float  blend_width = sqrt(static_cast<float>(dst_sz.area())) * blend_strength / 100.f;
	if(blend_width<10) blend_width = 10;
	if(blend_width>100) blend_width = 100;
	printf("blend width: %f \n", blend_width);
	
	//prepare
	MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
	mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));	
	blender->prepare(pts, sizes);
	//3.1 add the images one by one	
	Mat dilated_mask, seam_mask, mask, mask_warped,img_warped_s;
	int num_images = images.size();
	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		int ht = images[img_idx].rows;
		int wd = images[img_idx].cols;
		
		//mask_warped = masks[img_idx].clone();
		mask_warped.create(images[img_idx].size(), CV_8U);
		mask_warped.setTo(Scalar::all(255));
		//printf("%d \n", mask_warped.type());
		
		dilate(masks[img_idx], dilated_mask, Mat());
		resize(dilated_mask, seam_mask, mask_warped.size());
		mask_warped = seam_mask & mask_warped;
		
		images[img_idx].convertTo(img_warped_s, CV_16S);
		blender->feed(img_warped_s, mask_warped, pts[img_idx]);

		//img_warped_s.release();
		//mask_warped.release();
	}

	//3.2 do blend
	Mat result, result_mask;
	blender->blend(result, result_mask);
	Mat finalImage;
	result.convertTo(finalImage, CV_8U);
    
	int mht = result.rows;
	int mwd = result.cols;
	unsigned char* pbuffer = (unsigned char*)malloc(mht*mwd);
	memset(pbuffer, 0, mht*mwd);
	uchar* pdata = finalImage.data;	


	//output as tiff file
	char **papszOptions  = NULL;	
	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* poDatasetMosaic = poDriver->Create(outFilePath, mwd, mht, 1, GDT_Byte, papszOptions );
	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = maxrx/ratio;
	geoTransform[5] = -maxry/ratio;
	poDatasetMosaic->SetGeoTransform(geoTransform);
	poDatasetMosaic->SetProjection(geoArray[0].projectRef);
	
	GDALRasterBand* poBand = poDatasetMosaic->GetRasterBand( 1 );
	for(int j=0; j<mht; j++)
		for(int i=0; i<mwd; i++)
		{
			pbuffer[j*mwd+i] = pdata[j*3*mwd+3*i];
		}	
	poBand->RasterIO(GF_Write,0,0,mwd,mht,pbuffer,mwd,mht,GDT_Byte,0,0);
	
    GDALClose( (GDALDatasetH) poDatasetMosaic );
    
	free(pbuffer);

	return 1;
}

int CMultiBandBlender::GeoBlend( char** filenames, int nFile, char* outFilePath)
{
	vector<Mat> images;
	vector<Mat> masks;
	masks.resize(nFile);
	images.resize(nFile);

	vector<stGeoInfo> geoArray;
	geoArray.resize(nFile);

	GDALAllRegister();

	double ratio = 1.0;
	printf("Generateing Masks... \n");
	//1. generate masks
	for(int i=0; i<nFile; i++)
	{
		GenerateMask(filenames[i], masks[i], ratio);
		//imwrite("d:\\smallMask.jpg", masks[i]);
		GetImage(filenames[i], images[i], ratio);				
		GetGeoInformation( filenames[i], geoArray[i]);
	}

	/*
	for(int i=0; i<nFile; i++)
	{
		char file[256];
		sprintf(file, "d:\\initmask_%d.jpg", i);
		imwrite(file, masks[i]);
	}*/
	    
	//2. find seam
	double minx,maxx,miny,maxy;
	double maxrx,maxry;
	double rx,ry;
	minx = 5000000000;
	maxx = -5000000000;
	miny = 5000000000;
	maxy = -5000000000;
	maxrx = 0;
	maxry = 0;
	int zoneNumber = geoArray[0].zoneNumber;
	for(int i=0; i<nFile; i++)
	{	
		if(zoneNumber!=geoArray[i].zoneNumber)
			continue;
		
		rx = geoArray[i].dx;
		ry = fabs(geoArray[i].dy);
		maxrx = max(rx, maxrx);
		maxry = max(ry, maxry);

		minx = min(minx, geoArray[i].left);
		maxx = max(maxx, geoArray[i].left + geoArray[i].wd*rx);
		miny = min(miny, geoArray[i].top - geoArray[i].ht*ry);
		maxy = max(maxy, geoArray[i].top);
	}
	printf(" left:%lf right:%lf  top:%lf bottom:%lf \n", minx, maxx, maxy, miny);

    vector<Point> pts;
	pts.resize(nFile);
	for(int i=0; i<nFile; i++)
	{
		int x = (geoArray[i].left - minx) / maxrx * ratio;
		int y = (maxy - geoArray[i].top) / maxry * ratio;
		pts[i].x = x;
		pts[i].y = y;
		printf("%d %d \n", x, y);
	} 

	vector<Size> sizes;
	sizes.resize(images.size());
	for(int i=0; i<nFile; i++)
	{
		sizes[i].width  = images[i].cols;
		sizes[i].height = images[i].rows;
	}
	
	printf("Finding seam... \n");
	Ptr<SeamFinder> seam_finder = new detail::VoronoiSeamFinder();
	seam_finder->find(images, pts, masks);
	
	/*
	for(int i=0; i<nFile; i++)
	{
		char file[256];
		sprintf(file, "d:\\mask_%d.jpg", i);
		imwrite(file, masks[i]);
	}*/
	

	//3. define the blender
	Ptr<Blender> blender;
	int blend_type = Blender::MULTI_BAND;
	bool try_gpu   = false;
	blender = Blender::createDefault(blend_type, try_gpu);
	
	//get the area of mosaic
	Size   dst_sz = resultRoi(pts, sizes).size();
	
	//calculate the area of blend
	double blend_strength = 10;
	float  blend_width = sqrt(static_cast<float>(dst_sz.area())) * blend_strength / 100.f;
	if(blend_width<1) blend_width = 10;
	printf("blend width: %f \n", blend_width);
	
	//blend_width = 2;

	//prepare
	MultiBandBlender* mb = dynamic_cast<MultiBandBlender*>(static_cast<Blender*>(blender));
	mb->setNumBands(static_cast<int>(ceil(log(blend_width)/log(2.)) - 1.));	
	blender->prepare(pts, sizes);
	//3.1 add the images one by one	
	Mat dilated_mask, seam_mask, mask, mask_warped,img_warped_s;
	int num_images = images.size();
	for (int img_idx = 0; img_idx < num_images; img_idx++)
	{
		int ht = images[img_idx].rows;
		int wd = images[img_idx].cols;
		
		//mask_warped = masks[img_idx].clone();
		mask_warped.create(images[img_idx].size(), CV_8U);
		mask_warped.setTo(Scalar::all(255));
		//printf("%d \n", mask_warped.type());
		
		dilate(masks[img_idx], dilated_mask, Mat());
		resize(dilated_mask, seam_mask, mask_warped.size());
		mask_warped = seam_mask & mask_warped;
		
		images[img_idx].convertTo(img_warped_s, CV_16S);
		blender->feed(img_warped_s, mask_warped, pts[img_idx]);

		//img_warped_s.release();
		//mask_warped.release();
	}

	//3.2 do blend
	Mat result, result_mask;
	blender->blend(result, result_mask);
	Mat finalImage;
	result.convertTo(finalImage, CV_8U);
    
	int mht = result.rows;
	int mwd = result.cols;
	unsigned char* pbuffer = (unsigned char*)malloc(mht*mwd);
	memset(pbuffer, 0, mht*mwd);
	uchar* pdata = finalImage.data;	

	printf("Output... \n");
	//output as tiff file
	char **papszOptions  = NULL;	
	GDALDriver* poDriver = GetGDALDriverManager()->GetDriverByName("GTiff");
	GDALDataset* poDatasetMosaic = poDriver->Create(outFilePath, mwd, mht, 3, GDT_Byte, papszOptions );
	double  geoTransform[6];
	memset(geoTransform, 0, sizeof(double)*6);
	geoTransform[0] = minx;
	geoTransform[3] = maxy;
	geoTransform[1] = maxrx/ratio;
	geoTransform[5] = -maxry/ratio;
	poDatasetMosaic->SetGeoTransform(geoTransform);
	poDatasetMosaic->SetProjection(geoArray[0].projectRef);
	
	GDALRasterBand* poBand = poDatasetMosaic->GetRasterBand( 1 );
	for(int j=0; j<mht; j++)
		for(int i=0; i<mwd; i++)
		{
			pbuffer[j*mwd+i] = pdata[j*3*mwd+3*i];
		}	
	poBand->RasterIO(GF_Write,0,0,mwd,mht,pbuffer,mwd,mht,GDT_Byte,0,0);
	
	poBand = poDatasetMosaic->GetRasterBand( 2 );
	for(int j=0; j<mht; j++)
		for(int i=0; i<mwd; i++)
		{
			pbuffer[j*mwd+i] = pdata[j*3*mwd+3*i+1];
		}	
	poBand->RasterIO(GF_Write,0,0,mwd,mht,pbuffer,mwd,mht,GDT_Byte,0,0);

	poBand = poDatasetMosaic->GetRasterBand( 3 );
	for(int j=0; j<mht; j++)
		for(int i=0; i<mwd; i++)
		{
			pbuffer[j*mwd+i] = pdata[j*3*mwd+3*i+2];
		}	
	poBand->RasterIO(GF_Write,0,0,mwd,mht,pbuffer,mwd,mht,GDT_Byte,0,0);

    GDALClose( (GDALDatasetH) poDatasetMosaic );
    
	free(pbuffer);

	return 1;
}


// calculate the whole translation based on optical flow
int OptialcFlowMotion(IplImage* pFirst, IplImage* pSecond, double& tx, double& ty)
{
	int nMaxCount = 10000;

	int ht = pFirst->height;
	int wd = pFirst->width;
	int scanwd = pFirst->widthStep;

	IplImage* curPYR  = cvCreateImage( cvSize(wd, ht), 8, 1 );
	IplImage* nextPYR = cvCreateImage( cvSize(wd, ht), 8, 1 );

	//1. find good features (harris feature points)
	IplImage* eig = cvCreateImage( cvGetSize(pFirst), 32, 1 );
	IplImage* temp = cvCreateImage( cvGetSize(pFirst), 32, 1 );
	double quality = 0.01;
	double min_distance = 10;
	int win_size = 10;
	CvPoint2D32f* points = (CvPoint2D32f*)malloc(nMaxCount*sizeof(CvPoint2D32f));
	CvPoint2D32f* npoints = (CvPoint2D32f*)malloc(nMaxCount*sizeof(CvPoint2D32f));        
	int count = nMaxCount;
	cvGoodFeaturesToTrack( pFirst, eig, temp, points, &count,
		quality, min_distance, 0, 3, 0, 0.04 );
	cvReleaseImage( &eig );
	cvReleaseImage( &temp );

	//2. calculate optical flow
	char* status = (char*)malloc(count);
	int Level = 3;
	cvCalcOpticalFlowPyrLK( pFirst, pSecond, curPYR, nextPYR, 
		points, npoints, 
		count, 
		cvSize(10,10), 
		Level, 
		status, 0, 
		cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,20,0.03),
		0);

	// calculate the optical flow and the mean translation		
	tx = 0;
	ty = 0;
	int num = 0;
	for(int k=0; k<count; k++)
	{		
		if( status[k] )
		{
			tx += (points[k].x - npoints[k].x);		
			ty += (points[k].y - npoints[k].y);
			num ++;
		}
	}
	tx /= (double)(num);
	ty /= (double)(num);

	//just for debug
	//if(tx>0) tx=0;
	//if(ty>0) ty=0;

	free(status);
	cvReleaseImage(&curPYR);
	cvReleaseImage(&nextPYR);
	free(points);
	free(npoints);

	return 1;
}



//calculate sparse optical motion, written using c style 
int KLTOptialcFlow(IplImage* pFirst, IplImage* pSecond, vector<MyPointF> motionVector)
{
	int nMaxCount = 10000;
	int ht = pFirst->height;
	int wd = pFirst->width;
	int scanwd = pFirst->widthStep;

	IplImage* curPYR  = cvCreateImage( cvSize(wd, ht), 8, 1 );
	IplImage* nextPYR = cvCreateImage( cvSize(wd, ht), 8, 1 );

	//1. find good features (harris feature points)
	IplImage* eig = cvCreateImage( cvGetSize(pFirst), 32, 1 );
	IplImage* temp = cvCreateImage( cvGetSize(pFirst), 32, 1 );
	double quality = 0.01;
	double min_distance = 10;
	int win_size = 10;
	CvPoint2D32f* points = (CvPoint2D32f*)malloc(nMaxCount*sizeof(CvPoint2D32f));
	CvPoint2D32f* npoints = (CvPoint2D32f*)malloc(nMaxCount*sizeof(CvPoint2D32f));        
	int count = nMaxCount;
	cvGoodFeaturesToTrack( pFirst, eig, temp, points, &count,
		quality, min_distance, 0, 3, 0, 0.04 );
	cvReleaseImage( &eig );
	cvReleaseImage( &temp );

	//2. calculate optical flow
	char* status = (char*)malloc(count);
	int Level = 3;
	cvCalcOpticalFlowPyrLK( pFirst, pSecond, curPYR, nextPYR, 
		points, npoints, 
		count, 
		cvSize(10,10), 
		3, 
		status, 0, 
		cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,3,0.03),
		0);
	
	// calculate the optical flow and the mean translation		
	double tx = 0;
	double ty = 0;
	MyPointF mt;  
	for(int k=0; k<count; k++)
	{
		tx += (points[k].x - npoints[k].x);
		ty += (points[k].y - npoints[k].y);
		if( status[k] )
		{
			mt.x = points[k].x - npoints[k].x;
			mt.y = points[k].y - npoints[k].y;
			motionVector.push_back(mt);
		}
	}
	tx /= (double)(count);
	ty /= (double)(count);

	free(status);
	cvReleaseImage(&curPYR);
	cvReleaseImage(&nextPYR);
	free(points);
	free(npoints);

	return 1;
}


/*  convert IplImage to gray image
*/
void   IplImageToColorImage(IplImage* pImage, unsigned char** r, unsigned char** g, unsigned char** b, int* ht, int* wd)
{
	if(pImage==NULL)
		return;

	int xoff = pImage->roi->xOffset;
	int yoff = pImage->roi->yOffset;

	*ht = pImage->roi->height;
	*wd = pImage->roi->width;
	*r = (unsigned char*)malloc( (*ht)*(*wd) );
	*g = (unsigned char*)malloc( (*ht)*(*wd) );
	*b = (unsigned char*)malloc( (*ht)*(*wd) );


	//8 bits
	if(pImage->nChannels == 1)
	{
		int scanwd = pImage->widthStep;
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*b)[j*(*wd)+i] = (unsigned char)( pImage->imageData[(j+yoff)*scanwd+(i+xoff)] );
				(*g)[j*(*wd)+i] = (unsigned char)( pImage->imageData[(j+yoff)*scanwd+(i+xoff)] );
				(*r)[j*(*wd)+i] = (unsigned char)( pImage->imageData[(j+yoff)*scanwd+(i+xoff)] );
			}
	}

	//24 bits
	if(pImage->nChannels == 3)
	{
		int scanwd = pImage->widthStep;
		for(int j=0; j<(*ht); j++)
			for(int i=0; i<(*wd); i++)
			{
				(*b)[j*(*wd)+i] = (unsigned char)( pImage->imageData[(j+yoff)*scanwd+(i+xoff)*3] );
				(*g)[j*(*wd)+i] = (unsigned char)( pImage->imageData[(j+yoff)*scanwd+(i+xoff)*3+1] );
				(*r)[j*(*wd)+i] = (unsigned char)( pImage->imageData[(j+yoff)*scanwd+(i+xoff)*3+2] );
			}		
	}
}



// mosaic based on optical flow estimation, written by xdh, 2013.12.18
int DirectMosaic(char** filenames, int nfile, char* mosaicFileName)
{
	CvRect rect;
	rect.x = 66;
	rect.y = 94;
	rect.width = 650;
	rect.height = 460;
    
	vector<Point2DDouble> motion;

    //calculate the optical flow between two frames
	vector<double> vtx;
	vector<double> vty;
	vtx.push_back(0);
	vty.push_back(0);

	for(int i=0; i<nfile-1; i++)
	{
		printf("%s \n", filenames[i]);
		printf("%s \n", filenames[i+1]);
		IplImage* pPrevious = cvLoadImage(filenames[i], 0);
		IplImage* pNext  = cvLoadImage(filenames[i+1], 0);

		//crop image
		cvSetImageROI(pPrevious, rect);
		cvSetImageROI(pNext, rect);
		IplImage* pFirst  = cvCloneImage(pPrevious);
		IplImage* pSecond = cvCloneImage(pNext);
        
		cvSaveImage("d:\\first.jpg", pFirst);
		cvSaveImage("d:\\second.jpg", pSecond);

		//detect feature points
		double tx,ty;
		OptialcFlowMotion(pFirst, pSecond, tx, ty);
		printf("%lf %lf \n", tx, ty);
		
		//optical flow
		vtx.push_back(tx);
		vty.push_back(-ty);
	
		cvReleaseImage(&pPrevious);
		cvReleaseImage(&pNext);
		cvReleaseImage(&pFirst);
		cvReleaseImage(&pSecond);
	}
	
	for(int i=1; i<vtx.size(); i++)
	{
		vtx[i] += vtx[i-1];
		vty[i] += vty[i-1];
	}


	//save as geotiff image
	unsigned char* r;
	unsigned char* g;
	unsigned char* b;
	int ht,wd;
	
	int nZoneNumber = 50;
	OGRSpatialReference oSRS;
	oSRS.SetUTM(nZoneNumber);
	oSRS.SetWellKnownGeogCS("WGS84");	
	char    *pszWKT =NULL;  
	oSRS.exportToWkt( &pszWKT );  

	for(int i=0; i<nfile; i++)
	{
		IplImage* pImage = cvLoadImage(filenames[i]);
		cvSetImageROI(pImage, rect);
		IplImage* pCrop  = cvCloneImage(pImage);
		cvSaveImage("d:\\crop.jpg", pCrop);

		IplImageToColorImage(pCrop, &r, &g, &b, &ht, &wd);
	
		stGeoInfo geo;
		geo.left = vtx[i];
		geo.top  = vty[i] + 10000;
		geo.dx   = 1;
		geo.dy   = -1;
		geo.zoneNumber = nZoneNumber;
		geo.projectRef = pszWKT;

		char tifFile[256];
		strcpy(tifFile, filenames[i]);
		strcpy(tifFile+strlen(tifFile)-3, "\0");
		strcat(tifFile, "tif");		
		GdalWriteImageColor(tifFile, r, g, b, ht, wd, geo);
	
		cvReleaseImage(&pImage);
		cvReleaseImage(&pCrop);
		free(r);
		free(g);
		free(b);
	}


	
	return 1;
}






