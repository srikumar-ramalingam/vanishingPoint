//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Using the computed vp and extracted lines, we compute votes required for computing the junction scores
*/

#pragma once

#include "vp.h"
#include "geometry.h"

class Votes
{
public:
	double Votes_f,Votes_sf;
	int line_ids[10];	// stores the first 10 lines selected.
	Points epoint;

	quadInfo getQuadrantInfo(double vps[6]);
	int getQuadrant(quadInfo qi,double x,double y);
	Votes getVotes(vp vpL,double x1,double y1,double x2,double y2,int lineType,double Threshold,int fullScale,int shortScale,int markLines);
	Points getEP(double x,double y,quadInfo qi,double vps[6],int w,int h);
};


