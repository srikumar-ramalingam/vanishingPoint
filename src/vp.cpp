//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	VP detection and line clustering
*/

#include <stdio.h>
#include <stdlib.h>
#include "geometry.h"
#include "vp.h"
#include <math.h>
#include <vector>
#include "IndexSorter.h"
#include <algorithm>

#define MAX_LNDIST 2
#define COL_THRESHOLD 50
#define VPINF 10
#define NAN_PTHRESH (1/(10^15))
#define NAN_NTHRESH -(1/(10^15))
#define M_PI   3.14159265358979323846 /* pi */
#define VP_PROB 0.3
#define MIN_PROB_VP 0.5


using namespace std;


double *vp::computePtLineVote(double *lines,int noLines,double *Xpnts,int noPts)
{
	geometry gm;
	double ta,max1,x1,y1,x2,y2,x,y,length1;
	ta=M_PI/3;	
	max1=gm.maxLineLength(lines,noLines);

	double *VoteArr;
	lnVPDist lv;
	VoteArr=new double[noLines*noPts];
	for(int i=0;i<noLines;i++)
	{		
		x1=lines[7*i+0];y1=lines[7*i+1];x2=lines[7*i+2];y2=lines[7*i+3];
		length1=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
		for(int j=0;j<noPts;j++)
		{
			x=Xpnts[3*j+0];y=Xpnts[3*j+1];
			lv=lineVPdist(x1,y1,x2,y2,x,y);

			lv.theta=lv.theta*(M_PI/180);
			
			//if ((noPts==4)&&(j==2))
			//		ta=ta*2;

			if ((lv.theta<ta)&&(lv.vpinfchk==1))
				VoteArr[i*noPts+j]=1-(1/ta)*lv.theta;
			else
				VoteArr[i*noPts+j]=0;

			//if ((noPts==4)&&(j==2))
			//		ta=ta/2;

			VoteArr[i*noPts+j]=exp(-(1-VoteArr[i*noPts+j])*(1-VoteArr[i*noPts+j])/0.02);
			VoteArr[i*noPts+j]=VoteArr[i*noPts+j]*length1/max1;
			
		}
	}
	return VoteArr;

}

lnVPDist vp::lineVPdist(double x1,double y1,double x2,double y2,double vpx,double vpy)
{
	lnVPDist lv;
	double midx,midy,length1,slope1,slope2,theta,d;
	
	midx=(x1+x2)/2;
	midy=(y1+y2)/2;
	length1=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
	slope1=atan((y2-y1)/(x2-x1));
	
	//line segment s_dash
	slope2=atan((vpy-midy)/(vpx-midx));
	//slope2=(vp(1)-midpntl(1))/(vp(2)-midpntl(2));
	//angle between s and s_dash
	theta=(slope1-slope2);

	
	//plusPi=find(theta<(-pi/2));
	//minusPi=find(theta>(pi/2));
	//theta(plusPi)=theta(plusPi)+pi;
	//theta(minusPi)=theta(minusPi)-pi;
	//theta=theta*180/pi;
	//%theta=atand(abs((ones(size(slope2,1),1)*slope1-slope2(:))./(1+ones(size(slope2,1),1)*slope1.*slope2(:))));

	
	if (theta<(-M_PI/2))
		theta=theta+M_PI;
	
	if (theta>(M_PI/2))
		theta=theta-M_PI;

	theta=theta*180/M_PI;

	//%midpoint and slope of s and s_dash are same 
	//%check if vp lies on rotated s_dash 

	//d=sqrt((vp(:,1)-ones(size(vp,1),1)*midpntl(1)).^2+(vp(:,2)-ones(size(vp,1),1)*midpntl(2)).^2);
	//vpinfchk=zeros(size(vp,1),1);
	//vpinfchk(d > lengthl/2)=1;
	d=sqrt((vpx-midx)*(vpx-midx)+(vpy-midy)*(vpy-midy));
	if (d>(length1/2))
	{
		lv.vpinfchk=1;
	//	printf("x1y1x2y2=%lf,%lf,%lf,%lf\nd=%lf,length1=%lf\n",x1,y1,x2,y2,d,length1);
	}
	if (d<=(length1/2))
		lv.vpinfchk=0;
	
	lv.theta=theta;
	return lv;
}

vp vp::getVP(double *linesInp,int noLines,int w,int h)
{
	vp vpL;
	vpL.p=new double[4*noLines];
	geometry gm;
	double *lines=new double[7*noLines];
	
	for(int i=0;i<noLines;i++)
		for(int j=0;j<7;j++)
			lines[7*i+j]=linesInp[7*i+j];
	
	
	//Xpnts = ComputeIntersectionPoints(lines);
	Points Xpnts;
	Xpnts=gm.intersectLines(lines,noLines);

	//Computing votes for every point from all lines
	//VoteArr = ComputeLinePtVote(lines,Xpnts);
	double *VoteArr;
	VoteArr=new double[noLines*Xpnts.noPts];
	VoteArr=computePtLineVote(lines,noLines,Xpnts.pts,Xpnts.noPts);

	//Vote=sum(VoteArr,1);

	//get the first point & remove the lines of this point
	//[vv ii]=sort(Vote,'descend');
	//vp(1:2)=Xpnts(ii(1),1:2);
	
	double maxVotes,sumVotes;
	int bestPt;
	maxVotes=0;
	for(int i=0;i<Xpnts.noPts;i++)
	{
		sumVotes=0;
		for(int j=0;j<noLines;j++)
			sumVotes=sumVotes+VoteArr[j*Xpnts.noPts+i];
		
		if (sumVotes>=maxVotes)
		{
			maxVotes=sumVotes;
			bestPt=i;
		}
	}
	vpL.vps[0]=Xpnts.pts[3*bestPt+0];
	vpL.vps[1]=Xpnts.pts[3*bestPt+1];

	
	// compute max1.
	double x1,y1,x2,y2,length1,max1;
	max1=gm.maxLineLength(lines,noLines);
	
	//Vote1 = VoteArr(:,ii(1));
	//active_lines = find((Vote1*maxl./All_lines(:,7))<0.8);
	//inactive_lines = find((Vote1*maxl./All_lines(:,7))>=0.8);
	double *Vote1,*Vote1Active,*Vote1Inactive;
	Vote1=new double[noLines];
	Vote1Active=new double[noLines];
	Vote1Inactive=new double[noLines];

	for(int i=0;i<noLines;i++)
		Vote1[i]=VoteArr[i*Xpnts.noPts+bestPt];

	/*
	for(int i=0;i<noLines;i++)
	{
		x1=lines[7*i+0];y1=lines[7*i+1];x2=lines[7*i+2];y2=lines[7*i+3];
		length1=sqrt((x1-y1)*(x1-y1)+(x2-y2)*(x2-y2));
		if ((Vote1[i]*max1/length1)>=0.8)
			printf("(%d,%lf),\n",i,Vote1[i]*max1/length1);
		
	}
	printf("\n");
	*/

	Lines active_lines,inactive_lines;
	active_lines.lines=new double[7*noLines];
	inactive_lines.lines=new double[7*noLines];
	active_lines.noLines=0;
	inactive_lines.noLines=0;
	for(int i=0;i<noLines;i++)
	{
		x1=lines[7*i+0];y1=lines[7*i+1];x2=lines[7*i+2];y2=lines[7*i+3];
		length1=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
		if ((Vote1[i]*max1/length1)<0.8)
		{
			Vote1Active[active_lines.noLines]=Vote1[i];
			for(int j=0;j<7;j++)
				active_lines.lines[7*active_lines.noLines+j]=lines[7*i+j];
			active_lines.noLines++;
		}
		if ((Vote1[i]*max1/length1)>=0.8)
		{						
			Vote1Inactive[inactive_lines.noLines]=Vote1[i];
			for(int j=0;j<7;j++)
				inactive_lines.lines[7*inactive_lines.noLines+j]=lines[7*i+j];
			inactive_lines.noLines++;
		}
	}

	//Vote1 = [Vote1(active_lines);Vote1(inactive_lines)];
	//lines = All_lines(active_lines,:);
	for(int i=0;i<active_lines.noLines;i++)
		Vote1[i]=Vote1Active[i];

	for(int i=0;i<inactive_lines.noLines;i++)
		Vote1[active_lines.noLines+i]=Vote1Inactive[i];

	for(int i=0;i<active_lines.noLines;i++)
		for(int j=0;j<7;j++)
			lines[7*i+j]=active_lines.lines[7*i+j];
	
	noLines=active_lines.noLines;

	//%work with the remaining lines
	//Xpnts = ComputeIntersectionPoints(lines);
	
	//printf("no of intersection points before:%d\n",Xpnts.noPts);
	delete Xpnts.pts;
	Xpnts=gm.intersectLines(lines,noLines);
	//printf("no of intersection points after:%d\n",Xpnts.noPts);

	//VoteArr = ComputeLinePtVote([lines;All_lines(inactive_lines,:)],Xpnts);
	//Vote=sum(VoteArr(1:size(lines,1),:),1);
	VoteArr=new double[(active_lines.noLines+inactive_lines.noLines)*Xpnts.noPts];
	
	delete lines;
	lines=new double[(active_lines.noLines+inactive_lines.noLines)*7];
	for(int i=0;i<active_lines.noLines;i++)
		for(int j=0;j<7;j++)
			lines[i*7+j]=active_lines.lines[i*7+j];

	for(int i=0;i<inactive_lines.noLines;i++)
		for(int j=0;j<7;j++)
			lines[(active_lines.noLines+i)*7+j]=inactive_lines.lines[i*7+j];

	noLines=active_lines.noLines+inactive_lines.noLines;
	
	VoteArr=computePtLineVote(lines,noLines,Xpnts.pts,Xpnts.noPts);
	//printf("no of lines=%d,no of points=%d\n",noLines,Xpnts.noPts);
	double *Vote;
	Vote=new double[Xpnts.noPts];
	for(int i=0;i<Xpnts.noPts;i++)
	{
		Vote[i]=0;
		for(int j=0;j<active_lines.noLines;j++)
			Vote[i]=Vote[i]+VoteArr[j*Xpnts.noPts+i];
	}
			
	//[vv ii]=sort(Vote,'descend');
	//Vote = vv(:);
	//Xpnts=Xpnts(ii,:);
	//VoteArr = VoteArr(:,ii);
	//%Remove some of the points
	//[Xpnts,Vote,VoteArr] = RemoveRedundantPoints2(Xpnts,Vote,VoteArr,w,h);
	int *reOrder=new int[Xpnts.noPts];
	for(int i=0;i<Xpnts.noPts;i++)
		reOrder[i]=i;

	IndexSorter<double>::decSort(Vote, Xpnts.noPts, reOrder);
	
	double *Xpnts_S,*VoteArr_S;
	Xpnts_S=new double[Xpnts.noPts*3];
	VoteArr_S=new double[Xpnts.noPts*noLines];

	for(int i=0;i<Xpnts.noPts;i++)
	{		
		for(int j=0;j<3;j++)
			Xpnts_S[3*i+j]=Xpnts.pts[3*reOrder[i]+j];

		for(int j=0;j<noLines;j++)
			VoteArr_S[j*Xpnts.noPts+i]=VoteArr[j*Xpnts.noPts+reOrder[i]];
	}

	delete Xpnts.pts;
	Xpnts.pts=Xpnts_S;
	delete VoteArr;
	VoteArr=VoteArr_S;
	
	//[Xpnts,Vote,VoteArr] = RemoveRedundantPoints2(Xpnts,Vote,VoteArr,w,h);
	ptsVote pL;	
	pL=removeRedundantPoints(Xpnts.pts,Xpnts.noPts,noLines,Vote,VoteArr,w,h);
	Xpnts.pts=pL.Xpnts;
	Xpnts.noPts=pL.noPts;
	Vote=pL.Vote;
	VoteArr=pL.VoteArr;	

	int noCandidates=0;
	CalibData cd,bestcd;
	double vps[6];
	
	double bestVote,curVote;
	bestVote=0;
	int bestIndex=0;
	vps[0]=vpL.vps[0];vps[1]=vpL.vps[1];			
	for(int i=1;i<Xpnts.noPts;i=i+1)
		for(int j=0;j<i;j=j+1)
		{
			vps[2]=Xpnts.pts[3*i+0];vps[3]=Xpnts.pts[3*i+1];
			vps[4]=Xpnts.pts[3*j+0];vps[5]=Xpnts.pts[3*j+1];
			cd=orthoVP(vps,w,h);
			if (cd.orthochk)
			{
				curVote=0;
				for(int k=0;k<noLines;k++)
					curVote=curVote+max(Vote1[k],max(VoteArr[k*Xpnts.noPts+i],VoteArr[k*Xpnts.noPts+j]));

				if (curVote>bestVote)
				{
					bestVote=curVote;
					vpL.vps[2]=Xpnts.pts[3*i+0];
					vpL.vps[3]=Xpnts.pts[3*i+1];
					vpL.vps[4]=Xpnts.pts[3*j+0];
					vpL.vps[5]=Xpnts.pts[3*j+1];
					bestcd=cd;
				}
			}
		}
	//printf("detected vp:(%lf,%lf),(%lf,%lf),(%lf,%lf)\n",vpL.vps[0],vpL.vps[1],vpL.vps[2],vpL.vps[3],vpL.vps[4],vpL.vps[5]);
	
	double *vpNew=reOrderVPs(vpL.vps,w,h);
	
	vpL.p=getVPLineProbabilities(lines,noLines,vpNew);

	for(int i=0;i<6;i++)
		vpL.vps[i]=vpNew[i];

	vpL.lines=lines;
	vpL.noLines=noLines;
	return vpL;

}

double *vp::reOrderVPs(double vps[6],int w,int h)
{
		double mx,my,distVal[3],x1,y1,x2,y2,x3,y3,slope2,slope3;
		int zIndex,xIndex,yIndex,idx[3];
		double *vpNew;

		mx=w/2;
		my=h/2;

		distVal[0]=sqrt((vps[0]-mx)*(vps[0]-mx)+(vps[1]-my)*(vps[1]-my));
		distVal[1]=sqrt((vps[2]-mx)*(vps[2]-mx)+(vps[3]-my)*(vps[3]-my));
		distVal[2]=sqrt((vps[4]-mx)*(vps[4]-mx)+(vps[5]-my)*(vps[5]-my));

		idx[0]=0;idx[1]=1;idx[2]=2;
		IndexSorter<double>::Sort(distVal, 3, idx);
		
		x1=vps[2*idx[0]];y1=vps[2*idx[0]+1];
		x2=vps[2*idx[1]];y2=vps[2*idx[1]+1];
		x3=vps[2*idx[2]];y3=vps[2*idx[2]+1];

        slope2=abs((y2-y1)/(x2-x1));
        slope3=abs((y3-y1)/(x3-x1));
    
		zIndex=idx[0];
        if (slope2<slope3)
		{
            xIndex=idx[1];
            yIndex=idx[2];
		}
        else
		{
            xIndex=idx[2];
            yIndex=idx[1];
		}
		
		vpNew=new double[6];
		vpNew[0]=vps[2*xIndex];vpNew[1]=vps[2*xIndex+1];
		vpNew[2]=vps[2*yIndex];vpNew[3]=vps[2*yIndex+1];
		vpNew[4]=vps[2*zIndex];vpNew[5]=vps[2*zIndex+1];

		return vpNew;
}

double *vp::getVPLineProbabilities(double *lines,int noLines,double vps[6])
{
	double *probVP,*p;
	double temp,x1,y1,x2,y2,length1,max1;
	geometry gm;
	probVP=new double[3*4];

	probVP[3*0+0]=vps[0];probVP[3*0+1]=vps[1];
	probVP[3*1+0]=vps[2];probVP[3*1+1]=vps[3];
	probVP[3*2+0]=vps[4];probVP[3*2+1]=vps[5];
	probVP[3*3+0]=vps[0];probVP[3*3+1]=vps[1];

	p=computePtLineVote(lines,noLines,probVP,4);
	max1=gm.maxLineLength(lines,noLines);
	for(int i=0;i<noLines;i++)
	{
		x1=lines[7*i+0];y1=lines[7*i+1];x2=lines[7*i+2];y2=lines[7*i+3];
		length1=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));

		// normalize.
		for(int j=0;j<3;j++)
			p[4*i+j]=p[4*i+j]*max1/length1;

		p[4*i+3]=0;

		if (max(p[4*i+0],max(p[4*i+1],p[4*i+2]))<0.5)
			p[4*i+3]=1;

		temp=0;
		for(int j=0;j<4;j++)
			temp=temp+p[4*i+j];

		for(int j=0;j<4;j++)
			p[4*i+j]=p[4*i+j]/temp;

		//printf("line no:%d,%lf,%lf,%lf,%lf\n",i,p[4*i+0],p[4*i+1],p[4*i+2],p[4*i+3]);
	}
	return p;
}

/*
function [Xpnts2,Vote2,VoteArr2] = RemoveRedundantPoints(Xpnts,Vote,VoteArr,w,h)
return;
*/

ptsVote vp::removeRedundantPoints(double *Xpnts,int noPts,int noLines,double *Vote,double *VoteArr,int w,int h)
{
	ptsVote pV;
	
	double *dist,dist0,thres,x0,y0;
	dist=new double[noPts];

	int nRemainingPoints;
	nRemainingPoints = noPts;

	int i = 0;
	while (i<nRemainingPoints)	
	{
		dist0=((Xpnts[3*i]-w)*(Xpnts[3*i]-w)+(Xpnts[3*i+1]-h)*(Xpnts[3*i+1]-h))/((0.5*w)*(0.5*w)+(0.5*h)*(0.5*h));
		if (dist0 <  1)
			thres=10*10;
		if (dist0 >= 1)
			thres=20*20*dist0;
		
		int j = i+1;
		while (j<nRemainingPoints)		
		{
			double diff = (Xpnts[i*3] - Xpnts[j*3])*(Xpnts[i*3] - Xpnts[j*3]) + 
						  (Xpnts[i*3+1] - Xpnts[j*3+1])*(Xpnts[i*3+1] - Xpnts[j*3+1]);

			if ((diff*diff) < thres)
			{
				Xpnts[j*3] = Xpnts[(nRemainingPoints-1)*3];
				Xpnts[j*3+1] = Xpnts[(nRemainingPoints-1)*3+1];
				Xpnts[j*3+2] = Xpnts[(nRemainingPoints-1)*3+2];
				Vote[j]=Vote[nRemainingPoints-1];
				
				for(int k=0;k<noLines;k++)
					VoteArr[k*noPts+j]=VoteArr[k*noPts+(nRemainingPoints-1)];

				nRemainingPoints--;
			}						
			else
			{
				j++;
			}
		}
		i++;
	}

	pV.Xpnts=Xpnts;
	pV.noPts=nRemainingPoints;
	pV.noLines=noLines;
	pV.Vote=Vote;
	pV.VoteArr=new double[pV.noPts*pV.noLines];
	
	for(int i=0;i<pV.noLines;i++)
		for(int j=0;j<pV.noPts;j++)
			pV.VoteArr[i*pV.noPts+j]=VoteArr[i*noPts+j];

	return pV;
}

CalibData vp::orthoVP(double *vp,int w,int h)
{
	double vp1[2],vp2[2],vp3[2],vptemp[2];
	
	//orthochk=zeros(size(vp1s,1),1);
	int inf1,inf2,inf3,inds_fff,inds_ffi,inds_fii,inds_iii,orthochk=0;

	vp1[0]=vp[0];vp1[1]=vp[1];
	vp2[0]=vp[2];vp2[1]=vp[3];
	vp3[0]=vp[4];vp3[1]=vp[5];
	
	//inf conditions
	inf1 = (abs(vp1[0])>VPINF*w) || (abs(vp1[1])>VPINF*h);
	inf2 = (abs(vp2[0])>VPINF*w) || (abs(vp2[1])>VPINF*h);
	inf3 = (abs(vp3[0])>VPINF*w) || (abs(vp3[1])>VPINF*h);

	inds_fff = ((!inf1) && (!inf2) && (!inf3));
	
	// if there is only one infinity term then it goes to the third column.
	/*inds = find(inf1 & ~inf2 & ~inf3);
	temp = vp1s(inds,:);
	vp1s(inds,:) = vp3s(inds,:);
	vp3s(inds,:) = temp;*/
	if (inf1 && (!inf2) && (!inf3))
	{
		vptemp[0]=vp1[0];vptemp[1]=vp1[1];
		vp1[0]=vp3[0];vp1[1]=vp3[1];
		vp3[0]=vptemp[0];vp3[1]=vptemp[1];
	}

	/*inds = find(~inf1 & inf2 & ~inf3);
	temp = vp2s(inds,:);
	vp2s(inds,:) = vp3s(inds,:);
	vp3s(inds,:) = temp;*/
	if ((!inf1) && inf2 && (!inf3))
	{
		vptemp[0]=vp2[0];vptemp[1]=vp2[1];
		vp2[0]=vp3[0];vp2[1]=vp3[1];
		vp3[0]=vptemp[0];vp3[1]=vptemp[1];
	}

	inds_ffi = ((inf1+inf2+inf3)==1);

	/*
	inds = find(inf1 & ~inf2 & inf3);
	temp = vp2s(inds,:);
	vp2s(inds,:) = vp1s(inds,:);
	vp1s(inds,:) = temp;*/
	if (inf1 && (!inf2) && inf3)
	{
		vptemp[0]=vp2[0];vptemp[1]=vp2[1];
		vp2[0]=vp1[0];vp2[1]=vp1[1];
		vp1[0]=vptemp[0];vp1[1]=vptemp[1];
	}

	/*inds = find(inf1 & inf2 & ~inf3);
	temp = vp3s(inds,:);
	vp3s(inds,:) = vp1s(inds,:);
	vp1s(inds,:) = temp;*/
	if (inf1 & inf2 && (!inf3))
	{
		vptemp[0]=vp3[0];vptemp[1]=vp3[1];
		vp3[0]=vp1[0];vp3[1]=vp1[1];
		vp1[0]=vptemp[0];vp1[1]=vptemp[1];
	}
	//printf("%lf,%lf,%lf,%lf,%lf,%lf\n",vp[0],vp[1],vp[2],vp[3],vp[4],vp[5]);
	//printf("%lf,%lf,%lf,%lf,%lf,%lf\n",vp1[0],vp1[1],vp2[0],vp2[1],vp3[0],vp3[1]);

	inds_fii = ((inf1+inf2+inf3)==2);
	
	inds_iii = (inf1 && inf2 && inf3);

	//printf("fff:%d,ffi:%d,fii:%d,iii:%d\n",inds_fff,inds_ffi,inds_fii,inds_iii);
	
	// case when all three vps are finite.
	double Mats[3][3],A[2][2],b[2],fsqr,u0,v0,detA;
	
	if (inds_fff)
	{	    
		Mats[0][0] = vp1[0]+vp2[0];
		Mats[0][1] = vp1[1]+vp2[1];
		Mats[0][2] = vp1[0]*vp2[0]+vp1[1]*vp2[1];
		Mats[1][0] = vp1[0]+vp3[0];
		Mats[1][1] = vp1[1]+vp3[1];
		Mats[1][2] = vp1[0]*vp3[0]+vp1[1]*vp3[1];
		Mats[2][0] = vp3[0]+vp2[0];
		Mats[2][1] = vp3[1]+vp2[1];
		Mats[2][2] = vp3[0]*vp2[0]+vp3[1]*vp2[1];

		A[0][0] = Mats[0][0]-Mats[1][0]; A[0][1] = Mats[0][1]-Mats[1][1];
		A[1][0] = Mats[0][0]-Mats[2][0]; A[1][1] = Mats[0][1]-Mats[2][1];
		b[0] = Mats[0][2]-Mats[1][2]; b[1] = Mats[0][2]-Mats[2][2];
		detA = A[0][0]*A[1][1]-A[0][1]*A[1][0];
		u0 = (A[1][1]*b[0]-A[0][1]*b[1])/detA;
		v0 = (A[0][0]*b[1]-A[1][0]*b[0])/detA;
		fsqr = Mats[0][0]*u0+Mats[0][1]*v0-Mats[0][2]-u0*u0-v0*v0;

		orthochk=(u0 <= 0.7*w) && (u0 >= 0.3*w) && 
				 (v0 <= 0.7*h) && (v0 >= 0.3*h) && 
				 (fsqr > 0) && (fsqr<=(5000*5000));
	}

	double r,vec1[2],vec2[2],dot12,norm1,norm2;
	if (inds_ffi)
	{
    
		r=((w/2-vp1[0])*(vp2[0]-vp1[0])+(h/2-vp1[1])*(vp2[1]-vp1[1]))/
			((vp2[0]-vp1[0])*(vp2[0]-vp1[0])+(vp2[1]-vp1[1])*(vp2[1]-vp1[1]));

	    u0= vp1[0] + r*(vp2[0]-vp1[0]);
		v0= vp1[1] + r*(vp2[1]-vp1[1]);

	    fsqr=u0*(vp1[0]+vp2[0])+v0*(vp2[1]+vp1[1])-(vp1[0]*vp2[0]+vp2[1]*vp1[1]+u0*u0+v0*v0);
	    
	    vec1[0]=vp2[0]-vp1[0];vec1[1]=vp2[1]-vp1[1];
		vec2[0]=vp3[0];		  vec2[1]=vp3[1];
		
		dot12 = vec1[0]*vec2[0]+vec1[1]*vec2[1];//sum(vec1*vec2,2); 
		norm1 = sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]);//(sum(vec1*vec1,2).^.5);
		norm2 = sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]);

		orthochk=(r > 0) && (r < 1) && 
				 (u0 <= 0.7*w) && (u0 >= 0.3*w) && 
				 (v0 <= 0.7*h) && (v0 >= 0.3*h) && 
				 (fsqr > 0) && (fsqr<=(5000*5000)) && 
				 (abs(dot12/(norm1*norm2)) < 0.1);
	}

	
	if (inds_fii)
	{
    
    vec1[0]=vp2[0];		  vec1[1]=vp2[1];
	vec2[0]=vp3[0];		  vec2[1]=vp3[1];
    
	dot12 = vec1[0]*vec2[0]+vec1[1]*vec2[1];//sum(vec1*vec2,2); 
		norm1 = sqrt(vec1[0]*vec1[0]+vec1[1]*vec1[1]);//(sum(vec1*vec1,2).^.5);
		norm2 = sqrt(vec2[0]*vec2[0]+vec2[1]*vec2[1]);
    u0=vp1[0];
	v0=vp1[1];
	orthochk=(r > 0) && (r < 1) && 
				 (u0 <= 0.7*w) && (u0 >= 0.3*w) && 
				 (v0 <= 0.7*h) && (v0 >= 0.3*w) &&  
				 (abs(dot12/(norm1*norm2)) < 0.1);
	}
if (inds_iii)
    orthochk = 0;

/*
if ((orthochk==1)&&(inds_fii==1))
{
		printf("fff:%d,ffi:%d,fii:%d,iii:%d\n",inds_fff,inds_ffi,inds_fii,inds_iii);
		printf("%lf,%lf,%lf,%lf,%lf,%lf\n",vp1[0],vp1[1],vp2[0],vp2[1],vp3[0],vp3[1]);
		printf("u0=%lf,v0=%lf,w/2=%d,h/2=%d\n\n",u0,v0,w/2,h/2);
}*/

CalibData cd;
cd.orthochk=orthochk;
cd.inds_fff=inds_fff;
cd.inds_ffi=inds_ffi;
cd.inds_fii=inds_fii;
cd.inds_iii=inds_iii;
cd.u0=u0;
cd.v0=v0;
cd.vp1[0]=vp1[0];cd.vp1[1]=vp1[1];
cd.vp2[0]=vp2[0];cd.vp2[1]=vp2[1];
cd.vp3[0]=vp3[0];cd.vp3[1]=vp3[1];
cd.fsqr=fsqr;

return cd;
}


int vp::getLineType(double p0, double p1, double p2, double p3)
{
	if ((p0 >= (p1 + MIN_PROB_VP)) && (p0 >= (p2 + MIN_PROB_VP)) && (p0 >= (p3 + MIN_PROB_VP)))
		return 0;

	if ((p1 >= (p0 + MIN_PROB_VP)) && (p1 >= (p2 + MIN_PROB_VP)) && (p1 >= (p3 + MIN_PROB_VP)))
		return 1;

	if ((p2 >= (p0 + MIN_PROB_VP)) && (p2 >= (p1 + MIN_PROB_VP)) && (p2 >= (p3 + MIN_PROB_VP)))
		return 2;

	if ((p3 >= (p0 + MIN_PROB_VP)) && (p3 >= (p1 + MIN_PROB_VP)) && (p3 >= (p2 + MIN_PROB_VP)))
		return 3;
}
