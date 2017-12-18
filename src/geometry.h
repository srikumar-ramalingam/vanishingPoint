//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Line Processing, and other geometry operations
*/

#pragma once

class ptLnDist
{
public:
	double lambda,dist;
	//double projX,projY;
};

class Lines
{
public:
	double *lines;
	int noLines;
};

class Points
{
public:
	int noPts;
	double *pts;
};

class lnVPDist
{
public:
	double theta;
	int vpinfchk;
};

class ptsVote
{
public:
	double *Xpnts;
	int noPts;
	int noLines;
	double *Vote;
	double *VoteArr;
};

class CalibData
{
public:
	int orthochk;
	double vp1[2],vp2[2],vp3[2];
	int inds_fff,inds_ffi,inds_fii,inds_iii;
	double fsqr,u0,v0;
};

class lineParams
{
public:
	double a,b,c;
};

class geometry
{
public:
	double maxLineLength(double *lines,int noLines);
	ptLnDist ptLineDistance(double x1,double y1,double x2, double y2, double x, double y);
	int *removeNearbyLines(double *lines, int noLines);
	Lines mergeLines(double *lines1,double *lines2,int noLines1,int noLines2);
	Lines discardShortAndCollinearLines(double *inpLines,int noInpLines,double minLen,int colCheck);
	double *crossProduct(double *p1,double *p2,int noLines);
	Points intersectLines(double *lines,int noLines);
	double *normalizeVec(double *l,int dimX,int dimY);
	lineParams getLineParams(double x1,double y1,double x2,double y2);
	Points intersectLinesParams(lineParams l1,lineParams l2);
	Points xRectLine(lineParams l,int w,int h);
};


