
#include "xml.hpp"
#include "time.h"

#include <iostream>
#include <string>

//opencv2
#include "opencv2/core/core.hpp"


using namespace cv;
using namespace std;



int TestFileStorageWrite()
{
	string file = "d:\\test.yml";
	FileStorage fs(file, FileStorage::WRITE);
	fs << "frameCount" << 5;
	time_t rawtime; time(&rawtime);
	fs << "calibrationDate" << asctime(localtime(&rawtime));
	Mat cameraMatrix = (Mat_<double>(3,3) << 1000, 0, 320, 0, 1000, 240, 0, 0, 1); //又一种Mat初始化方式
	Mat distCoeffs = (Mat_<double>(5,1) << 0.1, 0.01, -0.001, 0, 0);
	fs << "cameraMatrix" << cameraMatrix << "distCoeffs" << distCoeffs;

	//features为一个大小为3的向量,其中每个元素由随机数x,y和大小为8的uchar数组组成
	fs << "features" << "[";
	for( int i = 0; i < 3; i++ )
	{
		int x = rand() % 640;
		int y = rand() % 480;
		uchar lbp = rand() % 256;
		fs << "{:" << "x" << x << "y" << y << "lbp" << "[:";
		for( int j = 0; j < 8; j++ )
			fs << ((lbp >> j) & 1);
		fs << "]" << "}";
	}
	fs << "]";
	fs.release();

	return 0;
}


int TestFileStorageRead()
{
	FileStorage fs("d:\test.yml", FileStorage::READ);

	//方式一: []操作符
	int frameCount = (int)fs["frameCount"];

	//方式二: FileNode::operator >>()
	string date;
	fs["calibrationDate"] >> date;

	Mat cameraMatrix2, distCoeffs2;
	fs["cameraMatrix"] >> cameraMatrix2;
	fs["distCoeffs"] >> distCoeffs2;

	//注意FileNodeIterator的使用, 似乎只能用一维数组去读取里面所有的数据
	FileNode features   = fs["features"];
	FileNodeIterator it = features.begin(), it_end = features.end();
	int idx = 0;
	std::vector<uchar> lbpval;
	for( ; it != it_end; ++it, idx++ )
	{
		cout << "feature #" << idx << ": ";
		cout << "x=" << (int)(*it)["x"] << ", y=" << (int)(*it)["y"] << ", lbp: (";
		(*it)["lbp"] >> lbpval;  //直接读出一维向量

		for( int i = 0; i < (int)lbpval.size(); i++ )
			cout << " " << (int)lbpval[i];
		cout << ")" << endl;
	}
	fs.release();

	return 0;
}