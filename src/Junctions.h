//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Using the computed vp and extracted lines, we compute junction scores that are used for finding the best box.
*/

#pragma once

#include "vp.h"
#include "geometry.h"


#define J_SUBSAMP 5
#define SCALE_LEN 100
#define SHORTJUNC_LEN 60
#define NAN_PTHRESH 1000000000000000.0
#define NAN_NTHRESH -1000000000000000.0
#define MIN_NEGVOTE 0
#define EPSILON 0.1


class boxLayout
{
public:
	double boundaryLines[8][4];
	// Polygon ordering 0-left,1-ceiling,2-right,3-floor,4-middle.
	Points poly[5];
	int signature[5][6];
	lineParams polyLines[5][6];
};

class Junctions
{
public:
	double *votesFull,*votesShort,*imgVotes;
	int *topLeft,*topRight,*bottomLeft,*bottomRight;
	int dimVotesX,dimVotesY,noTopCandidates;
	int xbest,ybest;
	
	void computeImgVotes(vp vpL,int w,int h);
	void findTopCandidates(int noCandidates,vp vpL);
	Points findBestBox(int noCandidates,vp vpL,int w,int h,char featureFile[50],int featureWriteMode);
	double computeCornerScore(int x,int y);
	double computeJunctionScore(int x,int y,int juncType[6],double highThreshold,double lowThreshold);
	double computeBoxScore(double box_x[4],double box_y[4],vp vpL,int w,int h);
	boxLayout getBoxLayout(double box_x[4],double box_y[4],vp vpL,int w,int h);
	Points getPoly(double boundaryLines[8][4],int index1,int index2,int axis_type,double threshold1,double threshold2,double x1,double y1,double x2,double y2);
};


