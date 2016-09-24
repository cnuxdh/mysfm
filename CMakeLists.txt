#generate cvlib


INCLUDE_DIRECTORIES(/usr/local/include  
										../
										../CoreDll
										../Corelib
										../sfm-driver
										../matrix
										../5point
										../imagelib
										../siftlib
										../ann_1.1_char
										../cminpack
										)


ADD_LIBRARY(cvlib
	absOri.cpp
	absOri.hpp
	panorama.cpp
	panorama.hpp
	relativepose.hpp
	relativepose.cpp
	register.cpp
	sift.cpp
	bundlerio.cpp
	feature.cpp
	cali.cpp
	distortion.cpp
	cvInvoke.cpp
	../corelib/matrix.cpp
	../corelib/matrix.h
	../corelib/LatLong-UTMconversion.cpp
	../corelib/LatLong-UTMconversion.h
	../corelib/CommonFuncs.cpp
	../corelib/CommonFuncs.h
	../corelib/commonfile.cpp
	../corelib/commonfile.h
	../coredll/main.cpp
	../coredll/main.h
	../coredll/vim_imgs.cpp
	../coredll/vim_imgs.h
	../coredll/vim_sift.cpp
	../coredll/keys2a.cpp
	../siftlib/siftmatch.cpp
	../siftlib/siftfeature.cpp
	../siftlib/kdtree.cpp
	../siftlib/minpq.cpp
	../siftlib/publicfunction.cpp
	../siftlib/ransac.cpp
	)