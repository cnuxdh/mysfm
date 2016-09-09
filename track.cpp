
//opencv
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/highgui/highgui.hpp>

using namespace cv;


#include "track.hpp"




/**
* Function to perform fast template matching with image pyramid
*/
void FastMatchTemplate(cv::Mat& srca,  // The reference image
					   cv::Mat& srcb,  // The template image
					   cv::Mat& dst,   // Template matching result
					   int maxlevel)   // Number of levels
{
	std::vector<cv::Mat> refs, tpls, results;

	// Build Gaussian pyramid
	cv::buildPyramid(srca, refs, maxlevel);
	cv::buildPyramid(srcb, tpls, maxlevel);

	cv::Mat ref, tpl, res;

	// Process each level
	for (int level = maxlevel; level >= 0; level--)
	{
		ref = refs[level];
		tpl = tpls[level];
		res = cv::Mat::zeros(ref.size() + cv::Size(1,1) - tpl.size(), CV_32FC1);

		if (level == maxlevel)
		{
			// On the smallest level, just perform regular template matching
			cv::matchTemplate(ref, tpl, res, CV_TM_CCORR_NORMED);
			//imwrite("d:\\tm.jpg", res);
		}
		else
		{
			// On the next layers, template matching is performed on pre-defined 
			// ROI areas.  We define the ROI using the template matching result 
			// from the previous layer.

			cv::Mat mask;
			cv::pyrUp(results.back(), mask);

			cv::Mat mask8u;
			mask.convertTo(mask8u, CV_8U);
			//imwrite("d:\\mask.jpg", mask);

			// Find matches from previous layer
			std::vector<std::vector<cv::Point> > contours;
			cv::findContours(mask8u, contours, CV_RETR_EXTERNAL, CV_CHAIN_APPROX_NONE);

			// Use the contours to define region of interest and 
			// perform template matching on the areas
			for (int i = 0; i < contours.size(); i++)
			{
				cv::Rect r = cv::boundingRect(contours[i]);
				cv::matchTemplate(
					ref(r + (tpl.size() - cv::Size(1,1))), 
					tpl, 
					res(r), 
					CV_TM_CCORR_NORMED
					);

			}
		}

		// Only keep good matches
		cv::threshold(res, res, 0.94,1., CV_THRESH_TOZERO);
        results.push_back(res);

		Mat img = res*255;
		imwrite("d:\\res.jpg", img);
	}
	res.copyTo(dst);
}


CFastTM::CFastTM()
{

}

CFastTM::~CFastTM()
{

}

int CFastTM::Match(char* srcFile, char* templateFile)
{

	cv::Mat ref = cv::imread(srcFile);
	cv::Mat tpl = cv::imread(templateFile);

	if (ref.empty() || tpl.empty())
		return -1;

	cv::Mat ref_gray, tpl_gray;
	cv::cvtColor(ref, ref_gray, CV_BGR2GRAY);
	cv::cvtColor(tpl, tpl_gray, CV_BGR2GRAY);

	cv::Mat dst;
	FastMatchTemplate(ref_gray, tpl_gray, dst, 2);

	while (true)
	{
		double minval, maxval;
		cv::Point minloc, maxloc;
		cv::minMaxLoc(dst, &minval, &maxval, &minloc, &maxloc);

		if (maxval >= 0.9)
		{
			cv::rectangle(ref, maxloc, 
				cv::Point(maxloc.x + tpl.cols, maxloc.y + tpl.rows), 
				CV_RGB(0,255,0), 2
				);

			cv::floodFill(
				dst, maxloc, 
				cv::Scalar(0), 0, 
				cv::Scalar(.1), 
				cv::Scalar(1.)
				);
		}
		else
			break;
	}

	//imwrite("d:\\res.jpg", ref);	
	namedWindow("result", CV_WINDOW_AUTOSIZE );  
	imshow("result", ref);

	return 1;
}

CBasicTM::CBasicTM()
{

}

CBasicTM::~CBasicTM()
{

}

int CBasicTM::Match(char* srcFile, char* templateFile)
{
	int match_method = 0;

	Mat img   = imread( srcFile, 1 );  
	Mat templ = imread( templateFile, 1 );  

	// 用于显示结果   
	Mat img_display;  
	img.copyTo( img_display );  

	Mat result;

	// 用于存储匹配结果的矩阵   
	int result_cols =  img.cols - templ.cols + 1;  
	int result_rows = img.rows - templ.rows + 1;  
	result.create( result_cols, result_rows, CV_32FC1 );  
    
	// 进行模板匹配   
	matchTemplate( img, templ, result, match_method );  
	

	// 归一化结果（方便显示结果）   
	normalize( result, result, 0, 1, NORM_MINMAX, -1, Mat() );  
    //result*=255;
	//imwrite("d:\\res.jpg", result);

	// 找到最佳匹配位置   
	double minVal;   
	double maxVal;   
	Point minLoc;   
	Point maxLoc;  
	Point matchLoc;  

	minMaxLoc( result, &minVal, &maxVal, &minLoc, &maxLoc, Mat() );   // 寻找result中的最大和最小值，以及它们所处的像素位置   

	// 使用SQDIFF和SQDIFF_NORMED方法时：值越小代表越相似   
	// 使用其他方法时：值越大代表越相似   
	if( match_method  == CV_TM_SQDIFF || match_method == CV_TM_SQDIFF_NORMED )  
	{  
		matchLoc = minLoc;   
	}  
	else  
	{   
		matchLoc = maxLoc;   
	}  

	// 显示匹配结果   
	rectangle( img_display, matchLoc, Point( matchLoc.x + templ.cols , matchLoc.y + templ.rows ), CV_RGB(255,0,0), 2, 8, 0 );  
	rectangle( result, matchLoc, Point( matchLoc.x + templ.cols , matchLoc.y + templ.rows ), Scalar::all(0), 2, 8, 0 );  

	char* image_window = "Source Image";  
	namedWindow( image_window, CV_WINDOW_AUTOSIZE );  
	imshow( image_window, img_display );  
	//imshow( result_window, result );  

	//imwrite("d:\\src.jpg", img_display);
	//imwrite("d:\\res.jpg", result);


	return 1;
}