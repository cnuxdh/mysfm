
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
	Mat cameraMatrix = (Mat_<double>(3,3) << 1000, 0, 320, 0, 1000, 240, 0, 0, 1); //��һ��Mat��ʼ����ʽ
	Mat distCoeffs = (Mat_<double>(5,1) << 0.1, 0.01, -0.001, 0, 0);
	fs << "cameraMatrix" << cameraMatrix << "distCoeffs" << distCoeffs;

	//featuresΪһ����СΪ3������,����ÿ��Ԫ���������x,y�ʹ�СΪ8��uchar�������
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

	//��ʽһ: []������
	int frameCount = (int)fs["frameCount"];

	//��ʽ��: FileNode::operator >>()
	string date;
	fs["calibrationDate"] >> date;

	Mat cameraMatrix2, distCoeffs2;
	fs["cameraMatrix"] >> cameraMatrix2;
	fs["distCoeffs"] >> distCoeffs2;

	//ע��FileNodeIterator��ʹ��, �ƺ�ֻ����һά����ȥ��ȡ�������е�����
	FileNode features   = fs["features"];
	FileNodeIterator it = features.begin(), it_end = features.end();
	int idx = 0;
	std::vector<uchar> lbpval;
	for( ; it != it_end; ++it, idx++ )
	{
		cout << "feature #" << idx << ": ";
		cout << "x=" << (int)(*it)["x"] << ", y=" << (int)(*it)["y"] << ", lbp: (";
		(*it)["lbp"] >> lbpval;  //ֱ�Ӷ���һά����

		for( int i = 0; i < (int)lbpval.size(); i++ )
			cout << " " << (int)lbpval[i];
		cout << ")" << endl;
	}
	fs.release();

	return 0;
}