

#ifndef CV_EXPORT_H
#define CV_EXPORT_H


#ifdef _WIN32
# define DLL_EXPORT __declspec( dllexport )
#else
# define DLL_EXPORT
#endif

//#define DLL_EXPORT  _declspec(dllexport)





#endif