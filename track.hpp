
#ifndef CV_TRACK_HPP
#define CV_TRACK_HPP

#include "export.hpp"


class DLL_EXPORT CTemplateMatchBase
{
public:
	CTemplateMatchBase(){}
	~CTemplateMatchBase(){}

	virtual int Match(char* srcFile, char* templateFile){return 0;}
};

class DLL_EXPORT CFastTM: public CTemplateMatchBase
{
public:
	CFastTM();
	~CFastTM();

private:
	int Match(char* srcFile, char* templateFile);

};



class DLL_EXPORT CBasicTM: public CTemplateMatchBase
{
public:
	CBasicTM();
	~CBasicTM();

private:
	int Match(char* srcFile, char* templateFile);

};


#endif