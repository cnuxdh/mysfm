#include"detect.hpp"
#include"time.h"


typedef struct stSegment
{
	int left, right;
	int row;
};

int DetectJYZLoss(Mat Img, vector<Rect>& lossRects)
{
	srand(time(NULL));
	int rindex = rand();

	//resize the image
	double rx = 640.0/double(Img.cols);
	double ry = rx;
	Mat resizeImg;
	resize(Img, resizeImg, Size(0,0), rx, ry);

	//blur the image
	blur(resizeImg, resizeImg, Size(3, 3));


	///////  0. color statistic of best hue for segmentation  




	///////  1. segmentation using Hue
	//hue
	int hmin = 150;
	int hmin_Max = 360;
	int hmax = 200;
	int hmax_Max = 360;
	//light  
	int lmin = 50;
	int lmin_Max = 255;
	int lmax = 255;
	int lmax_Max = 255;
	//sature  
	int smin = 0;
	int smin_Max = 255;
	int smax = 255;
	int smax_Max = 255;

	Mat bgr, hls;
	resizeImg.convertTo(bgr, CV_32FC3, 1.0 / 255, 0);
	cvtColor(bgr, hls, COLOR_BGR2HLS);

	//Mat dst = Mat::zeros(bgr.size(), CV_32FC3);
	Mat mask;
	inRange(hls, Scalar(hmin, lmin / float(lmin_Max), smin / float(smin_Max)),
		Scalar(hmax, lmax / float(lmax_Max), smax / float(smax_Max)),
		mask);

	char segfile[256];
	sprintf(segfile, "c:\\temp\\seg-color-%d.jpg", rindex);
	imwrite(segfile, mask);

	///////// 2. linear regression
	vector<Point> forePoints;
	for (int r = 0; r < bgr.rows; r++)
	{
		for (int c = 0; c < bgr.cols; c++)
		{
			if (mask.at<uchar>(r, c) == 255)
			{
				Point p;
				p.x = c;
				p.y = r;
				forePoints.push_back(p);
			}
		}
	}

	//forground area is small
	if (forePoints.size() < 100)
		return -1;

	cv::Vec4f line;
	cv::fitLine(forePoints,
		line,
		CV_DIST_HUBER,
		0,
		0.01,
		0.01);

	double cos_theta = line[0];
	double sin_theta = line[1];
	double x0 = line[2], y0 = line[3];
	double phi = atan2(sin_theta, cos_theta) + PI / 2.0;
	double rho = y0 * cos_theta - x0 * sin_theta;
	//std::cout << "phi = " << phi / PI * 180 << std::endl;
	//std::cout << "rho = " << rho << std::endl;
	//rotation
	double angle = phi / PI * 180 - 90;
	cv::Point2f center(mask.cols / 2, mask.rows / 2);
	cv::Mat  rot = cv::getRotationMatrix2D(center, angle, 1);
	cv::Rect bbox = cv::RotatedRect(center, mask.size(), angle).boundingRect();
	rot.at<double>(0, 2) += bbox.width / 2.0 - center.x;
	rot.at<double>(1, 2) += bbox.height / 2.0 - center.y;
	cv::Mat dst;
	cv::warpAffine(mask, dst, rot, bbox.size(), INTER_NEAREST);
	

	//srand(0);
	//int randid = rand();
	char outfile[256];
	sprintf(outfile, "c:\\temp\\rotate-%d.jpg", rindex);
	imwrite(outfile, dst);

	
	//for output
	//Mat rotateImg;
	//cv::warpAffine(resizeImg, rotateImg, rot, bbox.size());


	///////3.morphology
	Mat erodeStruct = getStructuringElement(MORPH_ELLIPSE, Size(5, 5));
	erode(dst, dst, erodeStruct);
	dilate(dst, dst, erodeStruct);
	Mat dilateStruct = getStructuringElement(MORPH_RECT, Size(2, 10));
	dilate(dst, dst, dilateStruct);

	char morphfile[256];
	sprintf(morphfile, "c:\\temp\\morph-%d.jpg", rindex);
	imwrite(morphfile, dst);

	///////4.scan the right edge point 
	vector<vector<int>> allCross;
	vector<vector<int>> crossMat;
	vector<int> segNumber;
	for (int r = 0; r < dst.rows; r++)
	{
		uchar *p = dst.ptr<uchar>(r);
		vector<int> lineCross;
		for (int c = 1; c < dst.cols - 1; c++)
		{
			if (p[c] == 0)
				continue;
			if (p[c - 1] == 0 && p[c]>0)
			{
				lineCross.push_back(c);
			}
		}
		allCross.push_back(lineCross);
		if (lineCross.size() < 3)
			continue;
		crossMat.push_back(lineCross);
		segNumber.push_back(lineCross.size());
	}
	sort(segNumber.begin(), segNumber.end());
	if (segNumber.size() < 1)
		return -1;
	int medianSegNumber = segNumber[segNumber.size() / 2];
	printf("median seg: %d \n", medianSegNumber);

	/////////5.calculate the distance of segment
	vector<vector<int>> crossDis;
	vector<int> segDis;
	for (int j = 0; j < crossMat.size(); j++)
	{
		if (crossMat[j].size()<(0.5*medianSegNumber) || crossMat[j].size() >(1.5*medianSegNumber))
			continue;

		vector<int> lineCrossDis;
		for (int i = 0; i < crossMat[j].size() - 1; i++)
		{
			int dis = crossMat[j][i + 1] - crossMat[j][i];
			lineCrossDis.push_back(dis);
			//printf("%d ", dis);
			segDis.push_back(dis);
		}
		//printf("\n");
		crossDis.push_back(lineCrossDis);
	}
	sort(segDis.begin(), segDis.end());
	int medianDis = segDis[segDis.size() / 2];
	double thresh = segDis[segDis.size() / 2] * 1.5;
	printf("medianDis  thresh: %d  %lf \n", medianDis, thresh);

	/////////6.generate histogram
	int nHist = dst.cols / medianDis + 1;
	vector<double> lossHist(nHist, 0);
	vector<stSegment> validSeg;
	for (int j = 0; j < allCross.size(); j++)
	{
		int ncp = allCross[j].size();

		//1.segment constraint
		if (ncp < (medianSegNumber*0.5) || ncp>(medianSegNumber*1.5))
			continue;

		for (int i = 0; i < ncp - 1; i++)
		{
			int dis = allCross[j][i + 1] - allCross[j][i];

			//2.distance constrait
			if (dis>thresh && dis<thresh * 3)
			{
				//in the loss area, background is more than foreground
				int fgSum = 0;
				int bgSum = 0;
				uchar *dp = dst.ptr<uchar>(j);
				for (int k = allCross[j][i]; k < allCross[j][i + 1]; k++)
				{
					if (dp[k]>0)
						fgSum++;
					else
						bgSum++;

				}

				//3.background & foreground constaint
				if (bgSum > fgSum && fgSum>4)
				{
					//generate histogram
					double hv = double(allCross[j][i] + allCross[j][i + 1])*0.5 / double(medianDis);
					int    ih = hv;
					if (ih >= (nHist - 1)) ih = nHist - 2;
					double dh = hv - ih;
					lossHist[ih] += (1 - dh);
					lossHist[ih + 1] += dh;

					stSegment seg;
					seg.left = allCross[j][i];
					seg.right = allCross[j][i + 1];
					seg.row = j;
					validSeg.push_back(seg);
				}
			}
		}
	}

	printf("histogram... \n");
	for (int i = 0; i < lossHist.size(); i++)
	{
		printf("%lf ", lossHist[i]);
	}
	printf("\n");

	////////7.determin and cluster the loss area by histogram
	for (int i = 0; i < lossHist.size() - 1; i++)
	{
		if ((lossHist[i] + lossHist[i + 1])>(medianDis*0.6))
		{
			Rect rec;
			int id = 0;
			for (int k = 0; k < validSeg.size(); k++)
			{
				double cx = double(validSeg[k].left + validSeg[k].right)*0.5;
				double dIndex = cx / double(medianDis);

				if (fabs(dIndex - i - 0.5) < 1)
				{
					if (id == 0){
						rec.x = validSeg[k].left;
						rec.y = validSeg[k].row;
						rec.width = validSeg[k].right - validSeg[k].left;
						rec.height = 1;
					}
					else{
						Rect curRec;
						curRec.x = validSeg[k].left;
						curRec.y = validSeg[k].row;
						curRec.width = validSeg[k].right - validSeg[k].left;
						curRec.height = 1;

						//cluster
						rec = rec | curRec;
					}
					id++;
				}
			}
			lossRects.push_back(rec);
		}
	}

	//printf("detect loss number: %d \n", lossRects.size());
	//for (int i = 0; i < lossRects.size(); i++)
	//{
	//	rectangle(rotateImg, lossRects[i], CV_RGB(255, 0, 0), 2);
	//}
	//imwrite("c:\\temp\\detecion.jpg", rotateImg);

	
	//printf("rotate the rect ...\n");
	//finaly, rotate the rect backward
	for (int i = 0; i < lossRects.size(); i++){
		/*Mat pt(2, 1, rot.type());
		pt.at<double>(0, 0) = lossRects[i].x;
		pt.at<double>(1, 0) = lossRects[i].y;
		Mat rpt = rot.inv()*pt;
		lossRects[i].x = rpt.at<double>(0, 0);
		lossRects[i].y = rpt.at<double>(1, 0);*/
		double cx = lossRects[i].x - dst.cols*0.5;
		double cy = lossRects[i].y - dst.rows*0.5;
		double ra = angle / 180.0*PI;
		double rx = cos(ra)*cx - sin(ra)*cy;
		double ry = sin(ra)*cx + cos(ra)*cy;
		lossRects[i].x = rx + mask.cols*0.5;
		lossRects[i].y = ry + mask.rows*0.5;
	}
	
	//resize the rect
	for (int i = 0; i < lossRects.size(); i++)
	{
		lossRects[i].x /= rx;
		lossRects[i].y /= ry;
		lossRects[i].width  /= rx;
		lossRects[i].height /= ry;
	}

	return 0;
}


