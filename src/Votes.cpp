//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Using the computed vp and extracted lines, we compute votes required for computing the junction scores
*/

#include <stdio.h>
#include <stdlib.h>
#include "geometry.h"
#include <math.h>
#include <vector>
#include "IndexSorter.h"
#include "vp.h"
#include "Votes.h"
#include <algorithm>
#include <iostream>

#define NAN_PTHRESH (1/(10^15))
#define NAN_NTHRESH -(1/(10^15))
#define PROB_THRESHOLD 0.3
#define M_PI   3.14159265358979323846 /* pi */


//using namespace std;

// quadInfo stores information that helps us to detect whether a given point is in the first, second, third or the fourth 
// quadrant- we chose four points around the z-vp that lies on the four quadrants. We consider the following two line equations.
// lineHor - line joining z-vp and x-vp, lineVer - line joining z-vp and y-vp
// given a point we check whether it is left or right side of a vertical line and top or bottom of a horizontal line. 
// based on the signature we decide which quadrant a given point lies. In order to do this we chose four representative points 
// in each of these quadrants. 
quadInfo Votes::getQuadrantInfo(double vps[6])
{
	double xh1,yh1,xh2,yh2,xv1,yv1,xv2,yv2,Delta;
	quadInfo qi;
	geometry gm;

	qi.lineHor=gm.getLineParams(vps[4],vps[5],vps[0],vps[1]); // line horizontal joining z-vp and x-vp
	qi.lineVer=gm.getLineParams(vps[4],vps[5],vps[2],vps[3]); // line vertical joining z-vp and y-vp
	
	Delta=5;
	xh1=vps[4]+Delta; // 5 units on the positive x direction.
	yh1=(-qi.lineHor.c-qi.lineHor.a*xh1)/qi.lineHor.b;

	xh2=vps[4]-Delta; // 5 unit on the negative x direction.
	yh2=(-qi.lineHor.c-qi.lineHor.a*xh2)/qi.lineHor.b;

	yv1=vps[5]-Delta; // 5 unit on the negative y direction.
	xv1=(-qi.lineVer.c-qi.lineVer.b*yv1)/qi.lineVer.a;

	yv2=vps[5]+Delta; // 5 unit on the positive y direction.
	xv2=(-qi.lineVer.c-qi.lineVer.b*yv2)/qi.lineVer.a;
	
	//Quadrant 1 - top right
	qi.x_[0]=(xh1+xv1)/2;qi.y_[0]=(yh1+yv1)/2;

	//Quadrant 2 - bottom right
	qi.x_[1]=(xh1+xv2)/2;qi.y_[1]=(yh1+yv2)/2;

	//Quadrant 3 - bottom left
	qi.x_[2]=(xh2+xv2)/2;qi.y_[2]=(yh2+yv2)/2;

	//Quadrant 4 - top left
	qi.x_[3]=(xh2+xv1)/2;qi.y_[3]=(yh2+yv1)/2;

	for(int i=0;i<4;i++)
	{
		//printf("x,y=%lf,%lf\n",qi.x_[i],qi.y_[i]);
		if ((qi.lineHor.a*qi.x_[i]+qi.lineHor.b*qi.y_[i]+qi.lineHor.c)>0)
			qi.signature[i][0]=1;
		else
			qi.signature[i][0]=0;

		if ((qi.lineVer.a*qi.x_[i]+qi.lineVer.b*qi.y_[i]+qi.lineVer.c)>0)
			qi.signature[i][1]=1;
		else
			qi.signature[i][1]=0;
		
		//printf("signature-i=%d-(%d,%d)\n",i,qi.signature[i][0],qi.signature[i][1]);
	}
	return qi;
}

int Votes::getQuadrant(quadInfo qi,double x,double y)
{
	int signature[2],quadrant;

		if ((qi.lineHor.a*x+qi.lineHor.b*y+qi.lineHor.c)>0)
			signature[0]=1;
		else
			signature[0]=0;

		if ((qi.lineVer.a*x+qi.lineVer.b*y+qi.lineVer.c)>0)
			signature[1]=1;
		else
			signature[1]=0;
	
		for(int i=0;i<4;i++)
			if ((signature[0]==qi.signature[i][0])&&(signature[1]==qi.signature[i][1]))
				quadrant=i;
	return quadrant;
}


Points Votes::getEP(double x,double y,quadInfo qi,double vps[6],int w,int h)
{
	geometry gm;
	int quad,quad_;
	lineParams lp;
	Points ep_,ep;

	// Get the quadrant where (x,y) resides.
	quad=getQuadrant(qi,x,y);
	ep.noPts=3;
	ep.pts=new double[3*6]; // First three points corresponds to forward direction and the remaining 3 points correspond to reverse direction.
	for(int i=0;i<3;i++)
	{
		lp=gm.getLineParams(x,y,vps[2*i+0],vps[2*i+1]);
		ep_=gm.xRectLine(lp,w,h);
		quad_=getQuadrant(qi,ep_.pts[3*1+0],ep_.pts[3*1+1]);
		if (i<2) // for the x and y vps, we find the ep that lies in a different quadrant.
		{
			ep.pts[3*i+0]=ep_.pts[(quad_!=quad)*3+0];
			ep.pts[3*i+1]=ep_.pts[(quad_!=quad)*3+1];
			ep.pts[3*i+2]=1;

			ep.pts[3*(3+i)+0]=ep_.pts[(quad_==quad)*3+0];
			ep.pts[3*(3+i)+1]=ep_.pts[(quad_==quad)*3+1];
			ep.pts[3*(3+i)+2]=1;
		}
		else // if we consider the third vp, then we want the ep to lie in the same quadrant.
		{
			ep.pts[3*i+0]=ep_.pts[(quad_==quad)*3+0];
			ep.pts[3*i+1]=ep_.pts[(quad_==quad)*3+1];
			ep.pts[3*i+2]=1;

			ep.pts[3*(3+i)+0]=ep_.pts[(quad_!=quad)*3+0];
			ep.pts[3*(3+i)+1]=ep_.pts[(quad_!=quad)*3+1];
			ep.pts[3*(3+i)+2]=1;
		}
	}
	return ep;
}

Votes Votes::getVotes(vp vpL,double x1,double y1,double x2,double y2,int lineType,double Threshold,int fullScale,int shortRange,int markLines)
{
	Votes vt;
	geometry gm;
	ptLnDist pld1,pld2;

	//(x1,y1) and (x2,y2) is the end points of the tube along which we collect the Votes_f and Votes_sf
	int check1,check2,check3,check4,uppLimit,lowLimit,noMarkedLines=0;
	double dist2ZVP,reverseRatio,lengthToEP;
	int *Scale,shortScale;

	if (markLines==1)
	{
		for(int i=0;i<10;i++)
			vt.line_ids[i]=-1;
	}

	lengthToEP=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	
	Scale=new int[fullScale];
	for(int i=0;i<fullScale;i++)
		Scale[i]=0;
			
	vt.epoint.noPts=1;
	vt.epoint.pts=new double[3*1];

	vt.epoint.pts[3*0+0] = (x2*shortRange + x1*(lengthToEP-shortRange))/(lengthToEP + NAN_PTHRESH);
	vt.epoint.pts[3*0+1] = (y2*shortRange + y1*(lengthToEP-shortRange))/(lengthToEP + NAN_PTHRESH);
	vt.epoint.pts[3*0+2] = 1;
	for(int i=0;i<vpL.noLines;i++)
	{
		if (((vpL.p[4*i+lineType]>=(vpL.p[4*i+0]+PROB_THRESHOLD))+
			 (vpL.p[4*i+lineType]>=(vpL.p[4*i+1]+PROB_THRESHOLD))+
			 (vpL.p[4*i+lineType]>=(vpL.p[4*i+2]+PROB_THRESHOLD)))>=2)
			{
				pld1=gm.ptLineDistance(x1,y1,x2,y2,vpL.lines[7*i+0],vpL.lines[7*i+1]);
				pld2=gm.ptLineDistance(x1,y1,x2,y2,vpL.lines[7*i+2],vpL.lines[7*i+3]);
				check1=(pld1.dist<Threshold);
				check2=(pld2.dist<Threshold);
				if (check1&&check2)
				{
					check3=((pld1.lambda>0)||(pld2.lambda>0)); // One of the end points is on the positive side of (x1,y1).
						if (check3)
						{
							check4=((pld1.lambda<=1)||(pld2.lambda<=1));
							if (check4)
							{
								uppLimit=(int)(fullScale*std::min(1.0,std::max(pld1.lambda,pld2.lambda)));
								lowLimit=(int)(fullScale*std::max(0.0,std::min(pld1.lambda,pld2.lambda)));
					
								for(int k=lowLimit;k<uppLimit;k++)
								Scale[k]=1;

								if ((markLines)&&(noMarkedLines<10))
								{
									vt.line_ids[noMarkedLines]=i;
									noMarkedLines++;
								}
							}
						}
					}
				
			}
	}
	
	
	shortScale=(int)(fullScale*shortRange/lengthToEP);

	vt.Votes_f=0;vt.Votes_sf=0;
	for(int i=0;i<fullScale;i++)
	{
		vt.Votes_f=vt.Votes_f+Scale[i];
		vt.Votes_sf=vt.Votes_sf+(i<shortScale)*Scale[i];
	}
	
	vt.Votes_f=lengthToEP*vt.Votes_f/fullScale;
	vt.Votes_sf=lengthToEP*vt.Votes_sf/fullScale;

	return vt;
}

