
#include"lane.hpp"



void makeFromVid(string path)
{
    Mat frame;
    VideoCapture cap(path); // open the video file for reading

    if (!cap.isOpened())  // if not success, exit program
        cout << "Cannot open the video file" << endl;

    //cap.set(CV_CAP_PROP_POS_MSEC, 300); //start the video at 300ms

    double fps = cap.get(CV_CAP_PROP_FPS); //get the frames per seconds of the video
    cout << "Input video's Frame per seconds : " << fps << endl;

    cap.read(frame);
    LaneDetect detect(frame);

    while (1)
    {
        bool bSuccess = cap.read(frame); // read a new frame from video
        if (!bSuccess)                   //if not success, break loop
        {
            cout << "Cannot read the frame from video file" << endl;
            break;
        }

        cvtColor(frame, frame, CV_BGR2GRAY);

        //start = clock();
        detect.nextFrame(frame);
        //stop =clock();
        // cout<<"fps : "<<1.0/(((double)(stop-start))/ CLOCKS_PER_SEC)<<endl;

        if (waitKey(10) == 27) //wait for 'esc' key press for 10 ms. If 'esc' key is pressed, break loop
        {
            cout << "video paused!, press q to quit, any other key to continue" << endl;
            if (waitKey(0) == 'q')
            {
                cout << "terminated by user" << endl;
                break;
            }
        }
    }
}

int main_test()
{
    makeFromVid("/home/yash/opencv-2.4.10/programs/output.avi");
    // makeFromVid("/home/yash/opencv-2.4.10/programs/road.m4v");
    waitKey(0);
    destroyAllWindows();
    return 0;
}


int DetectLane(Mat& image)
{
    //cap.read(frame);
    LaneDetect detect(image);
    detect.nextFrame(image);

    return 0;
}


//image: gray 
void fillholes(Mat& image)
{
    Mat fillimage = image.clone();
    int rows = image.rows;
    int cols = image.cols;
    for (int j = 1; j < rows - 1; j++)
    {
        uchar* p = image.ptr<uchar>(j);
        uchar* dp = fillimage.ptr<uchar>(j);
        for (int i = 1; i < cols - 1; i++)
        {
            int value = p[i];
            if (value > 0)
                continue;

            double sum = 0;
            int    nsum = 0;
            for (int m = -1; m <= 1; m++)
            {
                uchar* pl = image.ptr<uchar>(j + m);
                for (int n = -1; n <= 1; n++)
                {
                    int g = pl[i + n];
                    if (g > 0)
                    {
                        sum += g;
                        nsum++;
                    }
                }
            }

            if (nsum > 0)
            {
                double ave = sum / double(nsum);
                dp[i] = int(ave + 0.5);
            }
        }
    }

    image.release();
    image = fillimage.clone();
}


int CVDetectLane(string filepath)
{
    int npos = filepath.find_last_of("\\");
    int npos1 = filepath.find_last_of(".");
    string title = filepath.substr(npos + 1, npos1 - npos);
    string outpath = "c:\\temp\\";

    Mat& image = imread(filepath);

    //0.from rgb to gray 
    Mat gray;
    cvtColor(image, gray, CV_RGB2GRAY);
    imwrite(outpath+title+"_ori.jpg", gray);

    //1.fill holes
    fillholes(gray);
    fillholes(gray);
    //equalizeHist(gray, gray);
    imwrite(outpath + title + "_fill.jpg", gray);
    int rows = gray.rows;
    int cols = gray.cols;

    //2. binary image
    int lanewid = 12; //the width of lane 
    printf("%d %d \n", rows, cols);
    //Mat binaryImage(rows, cols, CV_8UC1, 0);
    Mat binaryImage = gray.clone();
    //imwrite("c:\\temp\\lane_fore.jpg", binaryImage);
    for (int j = lanewid; j < rows - lanewid; j++)
    {
        printf(".");
        uchar *p = gray.ptr<uchar>(j);
        uchar *dp = binaryImage.ptr<uchar>(j);
        for (int i = lanewid; i < cols - lanewid; i++)
        {
            dp[i] = 0;

            if (p[i] == 0)
            {
                continue;
            }

            double dx_left = p[i] - p[i - lanewid];
            double dx_right = p[i] - p[i + lanewid];

            if (dx_left > 4 && dx_right > 4)
            {
                dp[i] = 255;
            }
        }
    }
    imwrite(outpath + title + "_binary.jpg", binaryImage);

    //3. remove the noise
    medianBlur(binaryImage, binaryImage, 3);
    int dilate_size = 2;
    Mat dilate_element = getStructuringElement(MORPH_ELLIPSE,
        Size(2 * dilate_size + 1, 2 * dilate_size + 1),
        Point(dilate_size, dilate_size));
    dilate(binaryImage, binaryImage, dilate_element);
    imwrite(outpath + title + "_fore.jpg", binaryImage);

    //4. contours
    Mat drawing = Mat::zeros(gray.size(), CV_8UC3);
    vector<vector<Point>> contours(10000);
    std::vector<Vec4i> hierarchy(10000);
    findContours(binaryImage, contours, hierarchy,
        CV_RETR_EXTERNAL, CV_CHAIN_APPROX_SIMPLE, Point(0, 0));
    printf("contour number: %d \n", contours.size());

    RNG rng(12345);
    for (int i = 0; i < contours.size(); i++)
    {        
        RotatedRect rect = minAreaRect(contours[i]);
        //rectangle(drawing, rect, Scalar(255, 0, 0));
        //printf("size: %f %f \n", rect.size.height, rect.size.width);
        double ht = rect.size.height;
        double wd = rect.size.width;
        double laneLen = max(ht, wd);
        double ratio = max( ht/wd, wd/ht );
        double perimeter = 2 * (ht + wd);
        double len = arcLength(contours[i], true);
        
        if (laneLen < 48)
            continue;
        
        if (len > 1.8 * perimeter)
            continue;

        if(ratio>1.2)
        //if (contourArea(contours[i]) > 8 && arcLength(contours[i], false) > 16)
        //if (contourArea(contours[i]) > 8 && arcLength(contours[i], false) > 16)
        {
            printf("%d ", i);
            Scalar color = Scalar(rng.uniform(0, 255), rng.uniform(0, 255), rng.uniform(0, 255));
            drawContours(drawing, contours, i, color, 2, 8, hierarchy);
            //drawContours(drawing, poly, i, color, 2, 8, vector<Vec4i>(), 0, Point());
        }
    }

    imwrite(outpath + title + "_contours.jpg", drawing);


    return 0;
}
