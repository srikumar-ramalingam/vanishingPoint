//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Line Processing, and other geometry operations
*/

#include <stdio.h>
#include <stdlib.h>
#include "geometry.h"
#include <math.h>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <Eigen/LU>

#define MAX_LNDIST 3
#define COL_THRESHOLD 50
#define VPINF 10
#define NAN_PTHRESH 0.00000001
#define NAN_NTHRESH -0.00000001
#define M_PI   3.14159265358979323846 /* pi */

using namespace std;

double geometry::maxLineLength(double *lines,int noLines)
{
	double max1,length1;
	for(int i=0;i<noLines;i++)
	{			
		length1=sqrt((lines[7*i+0]-lines[7*i+2])*(lines[7*i+0]-lines[7*i+2])+
					 (lines[7*i+1]-lines[7*i+3])*(lines[7*i+1]-lines[7*i+3]));
		if (max1<length1)
			max1=length1;
	}
	return max1;
}

lineParams geometry::getLineParams(double x1,double y1,double x2,double y2)
{
	lineParams lp;
	double detCheck,normFactor,Threshold;
	Threshold=0.00001;

	detCheck=x1*y2-x2*y1;
	if (abs(detCheck)<NAN_PTHRESH)
	{
		lp.a=(y2-y1);
		lp.b=(-1)*(x2-x1);
		lp.c=0;
	}
	else
	{
		Eigen::Matrix2d A(2,2);
		A(0,0)=x1;A(0,1)=y1;
		A(1,0)=x2;A(1,1)=y2;
		Eigen::Vector2d b;
		b(0)=-1;
		b(1)=-1;
		Eigen::Vector2d x = A.inverse() * b;
		lp.a=x(0);
		lp.b=x(1);
		lp.c=1;
	}	
	normFactor=sqrt(lp.a*lp.a+lp.b*lp.b+lp.c*lp.c);
	lp.a=lp.a/normFactor;lp.b=lp.b/normFactor;lp.c=lp.c/normFactor;
	return lp;
}

Points geometry::xRectLine(lineParams lp,int w,int h)
{
	geometry gm;
	double x[4],y[4],bdSize=1,minPtDist=20;

	lineParams topLine,bottomLine,leftLine,rightLine;
	topLine.a=0.0;		topLine.b=1;		topLine.c=-1.0;
	bottomLine.a=0.0;	bottomLine.b=1;		bottomLine.c=(-1.0)*h;
	leftLine.a=1;		leftLine.b=0;		leftLine.c=-1.0;
	rightLine.a=1;		rightLine.b=0;		rightLine.c=(-1.0)*w;

	Points p;
	p=gm.intersectLinesParams(lp,topLine);
	x[0]=p.pts[0];y[0]=p.pts[1];

	p=gm.intersectLinesParams(lp,rightLine);
	x[1]=p.pts[0];y[1]=p.pts[1];

	p=gm.intersectLinesParams(lp,leftLine);
	x[2]=p.pts[0];y[2]=p.pts[1];

	p=gm.intersectLinesParams(lp,bottomLine);
	x[3]=p.pts[0];y[3]=p.pts[1];

	Points ep;
	ep.noPts=2;
	ep.pts=new double[2*3];
	int epCount=0;
	for(int i=0;i<4;i++)
	{
		if (epCount<2)
		{
			if ((x[i]>=(-1)*bdSize)&&(x[i]<=(w+bdSize))&&(y[i]>=(-1)*bdSize)&&(y[i]<=(h+bdSize)))
			{
				if (epCount==0)
				{
					ep.pts[3*epCount+0]=x[i];
					ep.pts[3*epCount+1]=y[i];
					ep.pts[3*epCount+2]=1;
					epCount=epCount+1;
				}
				if (epCount==1)
					if (((x[i]-ep.pts[3*0+0])*(x[i]-ep.pts[3*0+0])+(y[i]-ep.pts[3*0+1])*(y[i]-ep.pts[3*0+1]))>(minPtDist))
						{
							ep.pts[3*epCount+0]=x[i];
							ep.pts[3*epCount+1]=y[i];
							ep.pts[3*epCount+2]=1;
							epCount=epCount+1;
				}


			}
		}
	}

	return ep;
}
Points geometry::intersectLinesParams(lineParams l1,lineParams l2)
{	
	Points p;
	p.noPts=1;
	p.pts=new double[3];
	double xCheck;

	xCheck=((((l1.a-l2.a)*(l1.a-l2.a)+(l1.b-l2.b)*(l1.b-l2.b))>NAN_PTHRESH)&&
			(((l1.a+l2.a)*(l1.a+l2.a)+(l1.b+l2.b)*(l1.b+l2.b))>NAN_PTHRESH));

	if (xCheck>0)
	{
		Eigen::Matrix2d A(2,2);
		A(0,0)=l1.a;A(0,1)=l1.b;
		A(1,0)=l2.a;A(1,1)=l2.b;
		Eigen::Vector2d b;
		b(0)=-l1.c;
		b(1)=-l2.c;
		Eigen::Vector2d x = A.inverse() * b;
		p.pts[0]=x(0);
		p.pts[1]=x(1);
		p.pts[2]=1;
	}
	else
	{
		p.pts[0]=-10000000;
		p.pts[1]=-10000000;
		p.pts[2]=0;
	}
	return p;
}

ptLnDist geometry::ptLineDistance(double x1,double y1,double x2, double y2, double x, double y)
{
	ptLnDist p;
	p.lambda=((x-x1)*(x2-x1)+(y-y1)*(y2-y1))/((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
	p.dist=sqrt(pow((x1+p.lambda*(x2-x1)-x),2)+pow((y1+p.lambda*(y2-y1)-y),2));
	//p.projX=x1+p.lambda*(x2-x1);
	//p.projY=y1+p.lambda*(y2-y1);
	return p;
}

int *geometry::removeNearbyLines(double *lines, int noLines)
{
	double x1,y1,x2,y2,x1_,y1_,x2_,y2_,length1,length2;
	int *validLines;
	ptLnDist pL1,pL2;

	validLines=new int[noLines];
	for(int i=0;i<noLines;i++)
		validLines[i]=1;
	
	for(int i=0;i<noLines;i++)
		{
			for(int j=i+1;j<noLines;j++)				
			{
				if ((validLines[i]==1)&&(validLines[j]==1))
				{
					x1 = lines[7*i+0]; y1 = lines[7*i+1]; x2 = lines[7*i+2]; y2 = lines[7*i+3];
					x1_= lines[7*j+0]; y1_= lines[7*j+1]; x2_= lines[7*j+2]; y2_= lines[7*j+3];
				
					pL1=ptLineDistance(x1,y1,x2,y2,x1_,y1_);
					if (pL1.dist<(3*MAX_LNDIST))
					{
						pL2=ptLineDistance(x1,y1,x2,y2,x2_,y2_);
						if (pL2.dist<(3*MAX_LNDIST))
						{
							length1=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
							length2=sqrt((x1_-x2_)*(x1_-x2_)+(y1_-y2_)*(y1_-y2_));
							if (length1>length2)
								validLines[j]=0;
							else
								validLines[i]=0;
						}
					}
				}
			}			
	}
	return validLines;
}

Lines 
geometry::discardShortAndCollinearLines(double *inpLines,int noInpLines,double minLen,int colCheck)
{
	int length1;
	Lines outLines;
	geometry gm;
	outLines.lines=new double[7*noInpLines];
	outLines.noLines=0;
	int *validLines;
	if (colCheck==1)
		validLines=gm.removeNearbyLines(inpLines,noInpLines);
	else
	{
		validLines=new int[noInpLines];
		for(int i=0;i<noInpLines;i++)
			validLines[i]=1;
	}

	for(int i=0;i<noInpLines;i++)
		{
			length1=sqrt((inpLines[7*i+0]-inpLines[7*i+2])*(inpLines[7*i+0]-inpLines[7*i+2])+
						(inpLines[7*i+1]-inpLines[7*i+3])*(inpLines[7*i+1]-inpLines[7*i+3]));
			if ((length1>minLen)&&(validLines[i]==1))
			{
				for(int j=0;j<7;j++)
					outLines.lines[outLines.noLines*7+j]=inpLines[7*i+j];
				outLines.noLines++;
			}
		}		
	return outLines;
}

Lines geometry::mergeLines(double *lines1,double *lines2,int noLines1,int noLines2)
{
	double lambda1,lambda2,dist1,dist2;
	double x1,y1,x2,y2,x1_,y1_,x2_,y2_,px1,px2,py1,py2;
	int i,j,counter;
	int *validLines2;
	Lines mL;
	ptLnDist pL;

	validLines2=new int[noLines2];
	for(j=0;j<noLines2;j++)
		validLines2[j]=1;

		for(i=0;i<noLines1;i++)
		{
			for(j=0;j<noLines2;j++)				
			{
				x1 = lines1[7*i+0]; y1 = lines1[7*i+1];	x2 = lines1[7*i+2];	y2 = lines1[7*i+3];
				x1_= lines2[7*j+0]; y1_= lines2[7*j+1]; x2_= lines2[7*j+2]; y2_= lines2[7*j+3];
				
				pL=ptLineDistance(x1,y1,x2,y2,x1_,y1_);
				
				if (pL.dist<MAX_LNDIST)
				{
					lambda1=pL.lambda;
					dist1=pL.dist;
					px1=x1_;
					py1=y1_;
					pL=ptLineDistance(x1,y1,x2,y2,x2_,y2_);
					if (pL.dist<MAX_LNDIST)
					{
						lambda2=pL.lambda;
						dist2=pL.dist;
						px2=x2_;
						py2=y2_;
						// case 1 - no overlap
						//if (((lambda1<0)&&(lambda2<0))||(((lambda1>1)&&(lambda2>1))
						
						// case 2 - if line2 is the superset of line 1 ==> make line1=line2 and kill line2
						if (((lambda1<=0)&&(lambda2>=1))||((lambda1>=1)&&(lambda2<=0)))
						{
							lines1[7*i+0]=x1_;lines1[7*i+1]=y1_;lines1[7*i+2]=x2_;lines1[7*i+3]=y2_;
							lines2[7*j+0]=0;  lines2[7*j+1]=0;	lines2[7*j+2]=0;  lines2[7*j+3]=0;
							validLines2[j]=0;
						}
						
						// case 3 - if line2 is the subset of line 1 ==> keep line1 and kill line2
						if (((lambda1<=1)&&(lambda1>=0))&&((lambda2>=0)&&(lambda2<=1)))
						{
							lines1[7*i+0]=x1;lines1[7*i+1]=y1;lines1[7*i+2]=x2;lines1[7*i+3]=y2;
							lines2[7*j+0]=0;  lines2[7*j+1]=0;	lines2[7*j+2]=0;  lines2[7*j+3]=0;
							validLines2[j]=0;
						}
						// case 4 - partial overlap
						if ((lambda1<0)&&(lambda2>=0)&&(lambda2<=1))
						{
							lines1[7*i+0]=px1;lines1[7*i+1]=py1;lines1[7*i+2]=x2; lines1[7*i+3]=y2;
							lines2[7*j+0]=0;  lines2[7*j+1]=0;	lines2[7*j+2]=0;  lines2[7*j+3]=0;
							validLines2[j]=0;
						}
						// case 5 - partial overlap
						if ((lambda2<0)&&(lambda1>=0)&&(lambda1<=1))
						{
							lines1[7*i+0]=px2;lines1[7*i+1]=py2;lines1[7*i+2]=x2; lines1[7*i+3]=y2;
							lines2[7*j+0]=0;  lines2[7*j+1]=0;	lines2[7*j+2]=0;  lines2[7*j+3]=0;
							validLines2[j]=0;
						}
						// case 6 - partial overlap
						if ((lambda1>=0)&&(lambda1<=1)&&(lambda2>1))
						{
							lines1[7*i+0]=x1; lines1[7*i+1]=y1; lines1[7*i+2]=px2;lines1[7*i+3]=py2;
							lines2[7*j+0]=0;  lines2[7*j+1]=0;	lines2[7*j+2]=0;  lines2[7*j+3]=0;
							validLines2[j]=0;
						}
						// case 7 - partial overlap
						if ((lambda2>=0)&&(lambda2<=1)&&(lambda1>1))
						{
							lines1[7*i+0]=x1; lines1[7*i+1]=y1; lines1[7*i+2]=px1;lines1[7*i+3]=py1;
							lines2[7*j+0]=0;  lines2[7*j+1]=0;	lines2[7*j+2]=0;  lines2[7*j+3]=0;
							validLines2[j]=0;
						}
					}
				}
			}
		}
		mL.noLines=noLines1;
		for(j=0;j<noLines2;j++)
			mL.noLines=mL.noLines+validLines2[j];
		mL.lines=new double[7*mL.noLines];
		
		for(i=0;i<noLines1;i++)
			for(j=0;j<7;j++)
				mL.lines[i*7+j]=lines1[i*7+j];
				
		counter=0;
		for(i=0;i<noLines2;i++)
		{
			if (validLines2[i]==1)
			{
				for(j=0;j<7;j++)
					mL.lines[(noLines1+counter)*7+j]=lines2[i*7+j];
				counter++;
			}
		}
		
		return mL;
}

/*
%%%%%%Computing intersections of all the lines%%%%%%
p1 = [lines(:, [1 3]) ones(size(lines, 1), 1)];
p2 = [lines(:, [2 4]) ones(size(lines, 1), 1)];
% get plane normals for line segments
l = cross(p1, p2);
l = l ./ repmat(sqrt(sum(l.^2,2)), 1, 3);

[XX YY]=meshgrid(1:size(l,1));
ll1=l(XX(:),:);ll2=l(YY(:),:);
Xpnts=cross(ll1,ll2);

%[x1 y1 x2 y2 x3 y3] are colinear if x1(y2-y3)+x2(y3-y1)+x3(y1-y2)=0;
colchck=[lines(XX(:),1) lines(XX(:),3) lines(YY(:),1) lines(YY(:),3) lines(YY(:),2) lines(YY(:),4)];
colchck=colchck(:,1).*(colchck(:,4)-colchck(:,6))+colchck(:,3).*(colchck(:,6)-colchck(:,2))+...
    colchck(:,5).*(colchck(:,2)-colchck(:,4));

keepind=find(abs(colchck)>50);
Xpnts=Xpnts(keepind,:);
Xpnts=[Xpnts(:,1)./Xpnts(:,3) Xpnts(:,2)./Xpnts(:,3)];
*/
Points geometry::intersectLines(double *lines,int noLines)
{
	//p1 = [lines(:, [1 3]) ones(size(lines, 1), 1)];
	//p2 = [lines(:, [2 4]) ones(size(lines, 1), 1)];
	double *p1,*p2;
	p1=new double[noLines*3];
	p2=new double[noLines*3];
	for(int i=0;i<noLines;i++)
	{
		p1[3*i+0]=lines[7*i+0];
		p1[3*i+1]=lines[7*i+1];
		p1[3*i+2]=1;
		p2[3*i+0]=lines[7*i+2];
		p2[3*i+1]=lines[7*i+3];
		p2[3*i+2]=1;
	}
	// get plane normals for line segments
	//l = cross(p1, p2);
	//l = l ./ repmat(sqrt(sum(l.^2,2)), 1, 3);
	double *l;
	l=crossProduct(p1,p2,noLines);
	l=normalizeVec(l,3,noLines);

	//[XX YY]=meshgrid(1:size(l,1));
	//ll1=l(XX(:),:);ll2=l(YY(:),:);
	//Xpnts=cross(ll1,ll2);
	//[x1 y1 x2 y2 x3 y3] are colinear if x1(y2-y3)+x2(y3-y1)+x3(y1-y2)=0;
	//colchck=[lines(XX(:),1) lines(XX(:),3) lines(YY(:),1) lines(YY(:),3) lines(YY(:),2) lines(YY(:),4)];
	//colchck=colchck(:,1).*(colchck(:,4)-colchck(:,6))+colchck(:,3).*(colchck(:,6)-colchck(:,2))+colchck(:,5).*(colchck(:,2)-colchck(:,4));
	//keepind=find(abs(colchck)>50);
	//keepind=find(abs(colchck)>50);
	//Xpnts=Xpnts(keepind,:);
	//Xpnts=[Xpnts(:,1)./Xpnts(:,3) Xpnts(:,2)./Xpnts(:,3)];

	double *ll1,*ll2,x1,y1,x2,y2,x3,y3,colCheck;
	Points Xpnts;
	ll1=new double[3*noLines*noLines];
	ll2=new double[3*noLines*noLines];
	int noValidX=0;
	for(int i=0;i<noLines;i++)
	{
		for(int j=0;j<noLines;j++)
		{
			x1=lines[3*i]; y1=lines[3*i+1]; x2=lines[3*j]; y2=lines[3*j+1]; x3=lines[3*j+2];y3=lines[3*j+3];
			colCheck=x1*(y2-y3)+x2*(y3-y1)+x3*(y1-y2);
			if (colCheck>COL_THRESHOLD)
			{				
				ll1[3*noValidX+0]=l[3*i+0];
				ll1[3*noValidX+1]=l[3*i+1];
				ll1[3*noValidX+2]=l[3*i+2];

				ll2[3*noValidX+0]=l[3*j+0];
				ll2[3*noValidX+1]=l[3*j+1];
				ll2[3*noValidX+2]=l[3*j+2];
				noValidX++;
			}
		}
	}
	Xpnts.pts=crossProduct(ll1,ll2,noValidX);
	int notNanCtr=0;
	for(int i=0;i<noValidX;i++)
	{
		if ((Xpnts.pts[3*i+2]>NAN_PTHRESH)||(Xpnts.pts[3*i+2]<NAN_NTHRESH))
		{
			Xpnts.pts[3*notNanCtr+0]=Xpnts.pts[3*i+0]/Xpnts.pts[3*i+2];
			Xpnts.pts[3*notNanCtr+1]=Xpnts.pts[3*i+1]/Xpnts.pts[3*i+2];
			Xpnts.pts[3*notNanCtr+2]=Xpnts.pts[3*i+2]/Xpnts.pts[3*i+2];
			notNanCtr++;
		}
	//	else
	//		printf("%lf,%lf,%lf\n",Xpnts.pts[3*notNanCtr],Xpnts.pts[3*notNanCtr+1],Xpnts.pts[3*notNanCtr+2]);
	}
	Xpnts.noPts=notNanCtr;
	delete l,ll1,ll2;
	return Xpnts;
}

double *geometry::crossProduct(double *p1,double *p2,int noLines)
{
	double *l,a1,a2,a3,b1,b2,b3,len;
	l=new double[3*noLines];
	for(int i=0;i<noLines;i++)
	{
		a1=p1[3*i+0];a2=p1[3*i+1];a3=p1[3*i+2];
		b1=p2[3*i+0];b2=p2[3*i+1];b3=p2[3*i+2];
		l[3*i+0] = a2*b3 - a3*b2; 
		l[3*i+1] = a3*b1 - a1*b3;
		l[3*i+2] = a1*b2 - a2*b1; 
	}
	return l;
}

double *geometry::normalizeVec(double *l,int dimX,int dimY)
{
	double len;
	for(int i=0;i<dimY;i++)
	{
		len=0;
		for(int j=0;j<dimX;j++)
			len=len+l[dimX*i+j]*l[dimX*i+j];

		len=sqrt(len);
		if (len>0)
		{
			for(int j=0;j<dimX;j++)
				l[dimX*i+j]=l[dimX*i+j]/len;
		}
	}
	return l;
}
