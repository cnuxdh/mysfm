/* fgbg.cpp
 *
 * background model for foreground detection, 
 * created on: 2015.11.20
 *	   Author: Xie Donghai
 *
 */


#include "fgbg.hpp"


#include <opencv2/opencv.hpp>
#include <iostream>
using namespace cv;



int ConstructBgModel(char* videoFile)
{	
	initModule_video();
	setUseOptimized(true);
	setNumThreads(8);

	//string may be invalidate in Visual Studio 2005
	string algorithmName = "BackgroundSubtractor.GMG";
	
	Ptr<BackgroundSubtractorGMG> fgbg = Algorithm::create<BackgroundSubtractorGMG>("BackgroundSubtractor.GMG");
	if (fgbg.empty())
	{
		std::cerr << "Failed to create BackgroundSubtractor.GMG Algorithm." << std::endl;
		return -1;
	}
	fgbg->set("initializationFrames", 20);
	fgbg->set("decisionThreshold", 0.7);

	VideoCapture cap;
	cap.open(videoFile);
	
	if (!cap.isOpened())
	{
		std::cerr << "Cannot read video. Try moving video file to sample directory." << std::endl;
		return -1;
	}

	Mat frame, fgmask, segm;

	namedWindow("FG Segmentation", WINDOW_NORMAL);

	for (;;)
	{
		cap >> frame;

		if (frame.empty())
			break;

		(*fgbg)(frame, fgmask);

		frame.copyTo(segm);
		add(frame, Scalar(100, 100, 0), segm, fgmask);

		imshow("FG Segmentation", segm);

		int c = waitKey(30);
		if (c == 'q' || c == 'Q' || (c & 255) == 27)
			break;
	}

	return 0;
}

