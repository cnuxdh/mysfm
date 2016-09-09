#ifndef CV_VIDEO_HPP
#define CV_VIDEO_HPP


#include "export.hpp"


class DLL_EXPORT CVideoReadBase
{
public:
	CVideoReadBase(){}
	virtual ~CVideoReadBase(){}

	virtual int Init(){return 0;}
	virtual int Capture(){return 0;}
};




#endif