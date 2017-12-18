//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	VP detection and line clustering
*/

#pragma once

#include "geometry.h"

class quadInfo
{
public:
	double x_[4],y_[4];
	int signature[4][2];
	lineParams lineHor,lineVer;
};

class vp
{
public:	
	double vps[6];
	double *lines;
	int noLines;
	double *p;
	quadInfo qi;
	double maxLength;

	double *computePtLineVote(double *lines,int noLines,double *Xpnts,int noPts);
	lnVPDist lineVPdist(double x1,double y1,double x2,double y2,double vpx,double vpy);
	vp getVP(double *lines,int noLines,int w,int h);
	ptsVote removeRedundantPoints(double *Xpnts,int noPts,int noLines,double *Vote,double *VoteArr,int w,int h);
	CalibData orthoVP(double *vp,int w,int h);
	double *getVPLineProbabilities(double *lines,int noLines,double vps[6]);
	double *reOrderVPs(double vps[6],int w,int h);
	int getLineType(double p0, double p1, double p2, double p3);
};


