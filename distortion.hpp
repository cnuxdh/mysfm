/* 
 *  Copyright (c) 2008-2010  Noah Snavely (snavely (at) cs.cornell.edu)
 *    and the University of Washington
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 */

/* Distortion.h */

#ifndef __distortion_h__
#define __distortion_h__

#include "sfm.h"
#include "vector.h"

void InvertDistortion(int n_in, int n_out, double r0, double r1, 
                      double *k_in, double *k_out);

v2_t UndistortNormalizedPoint(v2_t p, camera_params_t c);


template<typename T>
int Undistort(T sx, T sy, T& dx, T& dy, T k1, T k2)
{
	T radius = sx * sx + sy * sy;

	dx = sx * (1.0 + radius * (k1 + radius * k2));
	dy = sy * (1.0 + radius * (k1 + radius * k2));

	return 0;
}


#endif /* __distortion_h__ */
