//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Rotation estimation from VP
*/

#pragma once
#include <Eigen/Dense>
#include <Eigen/LU>
#include "vp.h"
#include <math.h>

class rotation
{
public:
	Eigen::Matrix3d rot;
	void computeRotationFromVP(vp vp_,Eigen::Matrix3d Kmatrix);
};
