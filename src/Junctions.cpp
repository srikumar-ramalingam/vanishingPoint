//
// (c) MERL 2012 - 2013
//
/* 
	Created by: Srikumar Ramalingam	
	Using the computed vp and extracted lines, we compute junction scores that are used for finding the best box.
*/

#include <stdio.h>
#include <stdlib.h>
#include "geometry.h"
#include <math.h>
#include <vector>
#include "IndexSorter.h"
#include "vp.h"
#include "Votes.h"
#include "Junctions.h"
#include <algorithm>
#include <iostream>
#include <highgui.h>
#include <cv.h>
#include <cvaux.h>

void Junctions::computeImgVotes(vp vpL,int w,int h)
{
		Votes vt;
		int quad;
		int junctionYTypes[6];
		junctionYTypes[0]=1;junctionYTypes[1]=1;junctionYTypes[2]=1;junctionYTypes[3]=0;junctionYTypes[4]=0;junctionYTypes[5]=0;	
		dimVotesX=(int)(w/J_SUBSAMP);
		dimVotesY=(int)(h/J_SUBSAMP);
		votesFull=new double[6*dimVotesX*dimVotesY];
		votesShort=new double[6*dimVotesX*dimVotesY];
		imgVotes=new double[dimVotesX*dimVotesY];
		
		Points ep;	
		double Threshold=5,maxVotes_f[3],maxVotes_b[3],maxVoteValue=0;
		double votes0_f,votes0_sf,votes1_f,votes1_sf,votes2_f,votes2_sf;
		maxVotes_f[0]=0;maxVotes_f[1]=0;maxVotes_f[2]=0;
		maxVotes_b[0]=0;maxVotes_b[1]=0;maxVotes_b[2]=0;
				
		for(int i=0;i<dimVotesX;i++)
			for(int j=0;j<dimVotesY;j++)
			{
				ep=vt.getEP((double)(i*J_SUBSAMP),(double)(j*J_SUBSAMP),vpL.qi,vpL.vps,w,h);
				for(int k=0;k<3;k++)
				{
					vt=vt.getVotes(vpL,i*J_SUBSAMP,j*J_SUBSAMP,ep.pts[3*k+0],ep.pts[3*k+1],k,Threshold,SCALE_LEN,SHORTJUNC_LEN,0);
					
					votesFull[6*(dimVotesX*j+i)+k]=vt.Votes_f;

					votesShort[6*(dimVotesX*j+i)+k]=vt.Votes_sf;

					if (k!=2)
						vt=vt.getVotes(vpL,i*J_SUBSAMP,j*J_SUBSAMP,ep.pts[3*(k+3)+0],ep.pts[3*(k+3)+1],k,Threshold,SCALE_LEN,SHORTJUNC_LEN,0);
					else
						vt=vt.getVotes(vpL,i*J_SUBSAMP,j*J_SUBSAMP,vpL.vps[4],vpL.vps[5],k,Threshold,SCALE_LEN,SHORTJUNC_LEN,0);

					votesFull[6*(dimVotesX*j+i)+3+k]=vt.Votes_f;
					votesShort[6*(dimVotesX*j+i)+3+k]=std::max(0.0,vt.Votes_sf-MIN_NEGVOTE);
					
					if (votesFull[6*(dimVotesX*j+i)+k]>maxVotes_f[k])
						maxVotes_f[k]=votesFull[6*(dimVotesX*j+i)+k];

					if (votesFull[6*(dimVotesX*j+i)+3+k]>maxVotes_b[k])
						maxVotes_b[k]=votesFull[6*(dimVotesX*j+i)+3+k];
				}
			}
		double Votes_f1,Votes_f2,Votes_f3,Votes_b1,Votes_b2,Votes_b3,Votes_sf1,Votes_sf2,Votes_sf3;
		for(int i=0;i<dimVotesX;i++)
			for(int j=0;j<dimVotesY;j++)
			{	
				quad=vt.getQuadrant(vpL.qi,(double)(i*J_SUBSAMP),(double)(j*J_SUBSAMP));
				if ((quad==0)&&(quad==3))
				{
					Votes_f1=std::min(double(w/3),votesFull[6*(dimVotesX*j+i)+0]);
					Votes_f2=std::min(double(h/3),votesFull[6*(dimVotesX*j+i)+1]);
				}
				else
				{
					Votes_f1=std::min(double(w/3),votesFull[6*(dimVotesX*j+i)+0]);
					Votes_f2=std::min(double(h/3),votesFull[6*(dimVotesX*j+i)+1]);				
				}

				Votes_sf1=votesShort[6*(dimVotesX*j+i)+0];
				Votes_sf2=votesShort[6*(dimVotesX*j+i)+1];
				Votes_sf3=votesShort[6*(dimVotesX*j+i)+2];


				Votes_f3=votesFull[6*(dimVotesX*j+i)+2];
				Votes_b1=votesFull[6*(dimVotesX*j+i)+3+0];
				Votes_b2=votesFull[6*(dimVotesX*j+i)+3+1];
				Votes_b3=votesFull[6*(dimVotesX*j+i)+3+2];
				
				imgVotes[dimVotesX*j+i]=15*Votes_f1*Votes_f2*Votes_f3+
										15*Votes_sf1*Votes_sf2*Votes_sf3+
										 5*Votes_f1*Votes_f2+
										 5*Votes_f2*Votes_f3+
										 5*Votes_f1*Votes_f3+
										 5*Votes_sf1*Votes_sf2+
										 5*Votes_sf2*Votes_sf3+
										 5*Votes_sf1*Votes_sf3-
										 10*Votes_b1*Votes_b1*Votes_b1*Votes_b1-
										 10*Votes_b2*Votes_b2*Votes_b2*Votes_b2-
										 10*Votes_b3*Votes_b3*Votes_b3*Votes_b3;
				if (imgVotes[dimVotesX*j+i]>maxVoteValue)
				{
					xbest=i*J_SUBSAMP;
					ybest=j*J_SUBSAMP;
					maxVoteValue=imgVotes[dimVotesX*j+i];
				}
			}
			printf("Best Score:%lf\n",maxVoteValue);
}

void Junctions::findTopCandidates(int noCandidates,vp vpL)
{
	int *reOrder=new int[dimVotesX*dimVotesY];
	int x,y;
	int quad,counter,cornersFound[4];
	Votes vt;
	double *Votes=new double[dimVotesX*dimVotesY];
	
	for(int i=0;i<(dimVotesX*dimVotesY);i++)
		reOrder[i]=i;
	
	// Top right best candidates.

	for(int i=0;i<dimVotesX;i++)
			for(int j=0;j<dimVotesY;j++)
	{
			Votes[j*dimVotesX+i]=imgVotes[j*dimVotesX+i];
	}

	// We divide the whole image to have 100 blocks and we allow only one non-zero value in each block.
	int hSize=(int)(dimVotesX/41);
	int vSize=(int)(dimVotesY/41);
	int chosenX,chosenY;
	double maxValue=NAN_NTHRESH;
	
	for(int i=0;i<dimVotesX;i=i+hSize)
		for(int j=0;j<dimVotesY;j=j+vSize)
		{
			maxValue=NAN_NTHRESH;
			for(int k=i;k<std::min(dimVotesX,(i+hSize));k++)
				for(int l=j;l<std::min(dimVotesY,(j+vSize));l++)
				{
					Votes[dimVotesX*l+k]=std::min(0.0,imgVotes[dimVotesX*l+k]);
					if (imgVotes[dimVotesX*l+k]>maxValue)
					{
						maxValue=imgVotes[dimVotesX*l+k];
						chosenX=k;
						chosenY=l;
					}
				}
			Votes[chosenY*dimVotesX+chosenX]=imgVotes[chosenY*dimVotesX+chosenX];
		}
		

	IndexSorter<double>::decSort(Votes, dimVotesX*dimVotesY, reOrder);
	topRight=new int[noCandidates];
	topLeft=new int[noCandidates];
	bottomLeft=new int[noCandidates];
	bottomRight=new int[noCandidates];

	counter=0;
	for(int i=0;i<4;i++)
		cornersFound[i]=0;

	while((counter<(dimVotesX*dimVotesY))&&((cornersFound[0]<noCandidates)||(cornersFound[1]<noCandidates)||(cornersFound[2]<noCandidates)||(cornersFound[3]<noCandidates)))
	{
		y=(int)(reOrder[counter]/dimVotesX);
		x=reOrder[counter]-y*dimVotesX;

		quad=vt.getQuadrant(vpL.qi,(double)(x*J_SUBSAMP),(double)(y*J_SUBSAMP));
		if (cornersFound[quad]<noCandidates)
		{
			if (quad==0)
				topRight[cornersFound[quad]]=reOrder[counter];

			if (quad==1)
				bottomRight[cornersFound[quad]]=reOrder[counter];

			if (quad==2)
				bottomLeft[cornersFound[quad]]=reOrder[counter];

			if (quad==3)
				topLeft[cornersFound[quad]]=reOrder[counter];

			cornersFound[quad]=cornersFound[quad]+1;
		}
		counter=counter+1;
	}
	delete Votes,reOrder;
}

Points Junctions::findBestBox(int noCandidates,vp vpL,int w,int h, char featureFile[50],int featureWriteMode)
{
	double bestScore=NAN_NTHRESH,curScore,maxWinScore,bestCornerScores[4],curCornerScores[4];
	double box_x[4],box_y[4];

	lineParams lineHor_i,lineHor_j,lineVer_i,lineVer_j;
	geometry gm;
	Points p,pbest;
	int x,y,chosenPair;
	pbest.noPts=4;
	pbest.pts=new double[3*4];

	// Selecting the pair from top right corner and bottom left corner
	for(int i=0;i<noCandidates;i++)
		for(int j=0;j<noCandidates;j++)
		{
			
			y=(int)(topLeft[i]/dimVotesX);
			x=topLeft[i]-y*dimVotesX;			
			box_x[3]=(double)(x*J_SUBSAMP);
			box_y[3]=(double)(y*J_SUBSAMP);			
		
			lineHor_i=gm.getLineParams(box_x[3],box_y[3],vpL.vps[0],vpL.vps[1]);
			lineVer_i=gm.getLineParams(box_x[3],box_y[3],vpL.vps[2],vpL.vps[3]);
			
			y=(int)(bottomRight[j]/dimVotesX);
			x=bottomRight[j]-y*dimVotesX;			
			box_x[1]=(double)(x*J_SUBSAMP);
			box_y[1]=(double)(y*J_SUBSAMP);
			
			lineHor_j=gm.getLineParams(box_x[1],box_y[1],vpL.vps[0],vpL.vps[1]);
			lineVer_j=gm.getLineParams(box_x[1],box_y[1],vpL.vps[2],vpL.vps[3]);
			
			p=gm.intersectLinesParams(lineVer_i,lineHor_j); // Intersection of a vert line from topleft and hor line from bottom right.
			box_x[2]=p.pts[0];
			box_y[2]=p.pts[1];
			
			p=gm.intersectLinesParams(lineHor_i,lineVer_j);
			box_x[0]=p.pts[0];
			box_y[0]=p.pts[1];
			
			curScore=computeBoxScore(box_x,box_y,vpL,w,h);
			if (curScore>bestScore)
			{
				bestScore=curScore;
				for(int k=0;k<4;k++)
				{
					pbest.pts[3*k+0]=box_x[k];
					pbest.pts[3*k+1]=box_y[k];
					chosenPair=1;
				}
			}
		}
		

	// Selecting the pair from top right corner and bottom left corner
	for(int i=0;i<noCandidates;i++)
		for(int j=0;j<noCandidates;j++)
		{
			y=(int)(topRight[i]/dimVotesX);
			x=topRight[i]-y*dimVotesX;
			box_x[0]=(double)(x*J_SUBSAMP);
			box_y[0]=(double)(y*J_SUBSAMP);
			
			lineHor_i=gm.getLineParams(box_x[0],box_y[0],vpL.vps[0],vpL.vps[1]);
			lineVer_i=gm.getLineParams(box_x[0],box_y[0],vpL.vps[2],vpL.vps[3]);
			
			y=(int)(bottomLeft[j]/dimVotesX);
			x=bottomLeft[j]-y*dimVotesX;			
			box_x[2]=(double)(x*J_SUBSAMP);
			box_y[2]=(double)(y*J_SUBSAMP);
			
			lineHor_j=gm.getLineParams(box_x[2],box_y[2],vpL.vps[0],vpL.vps[1]);
			lineVer_j=gm.getLineParams(box_x[2],box_y[2],vpL.vps[2],vpL.vps[3]);
			p=gm.intersectLinesParams(lineVer_i,lineHor_j);
			box_x[1]=p.pts[0];
			box_y[1]=p.pts[1];
	
			p=gm.intersectLinesParams(lineHor_i,lineVer_j);
			box_x[3]=p.pts[0];
			box_y[3]=p.pts[1];

			
			curScore=computeBoxScore(box_x,box_y,vpL,w,h);
			if (curScore>bestScore)
			{
				bestScore=curScore;
				for(int k=0;k<4;k++)
				{
					pbest.pts[3*k+0]=box_x[k];
					pbest.pts[3*k+1]=box_y[k];
					chosenPair=0;
				}
			}
		}
		if (chosenPair==0)
			printf("ChosenPair: top right and bottom left\n");
		if (chosenPair==1)
			printf("ChosenPair: top left and bottom right\n");
		
		return pbest;
}

double 
Junctions::computeCornerScore(int x,int y)
{
	
	double maxWinScore=NAN_NTHRESH,curScore;
	int winSize=2;
	if ((x>0)&&(x<dimVotesX)&&(y>0)&&(y<dimVotesY))
		maxWinScore=imgVotes[y*dimVotesX+x];
	
	for(int Delta_x=-(winSize-1);Delta_x<(winSize-1);Delta_x++)
		for(int Delta_y=-(winSize-1);Delta_y<(winSize-1);Delta_y++)
			if (((x+Delta_x)>0)&&((x+Delta_x)<dimVotesX)&&((y+Delta_y)>0)&&((y+Delta_y)<dimVotesY))
			{
				curScore=imgVotes[dimVotesX*(y+Delta_y)+(x+Delta_x)];
				if (maxWinScore<curScore)
					maxWinScore=curScore;
			}
	return maxWinScore;
	
}

double 
Junctions::computeBoxScore(double box_x[4],double box_y[4],vp vpL,int w,int h)
{
	double boxScore=0,juncYThreshold=SHORTJUNC_LEN*SHORTJUNC_LEN;
	double boundaryLines[8][4],x,y,Delta_x,Delta_y,lineLength,juncYScore;
	int junctionYTypes[6];
	junctionYTypes[0]=1;junctionYTypes[1]=1;junctionYTypes[2]=1;junctionYTypes[3]=0;junctionYTypes[4]=0;junctionYTypes[5]=0;
	Votes votes_;
	Points ep;
	
	for(int i=0;i<4;i++)
	{
		boxScore=boxScore+computeCornerScore((int)(box_x[i]/J_SUBSAMP),(int)(box_y[i]/J_SUBSAMP));				
	}
	
	/*
	for(int i=0;i<4;i++)
	{
		boundaryLines[i][0]=box_x[i];
		boundaryLines[i][1]=box_y[i];
		boundaryLines[i][2]=box_x[(i+1)%4];
		boundaryLines[i][3]=box_y[(i+1)%4];
	}
	for(int i=4;i<8;i++)
	{
		boundaryLines[i][1]=box_y[i-4];
		ep=votes_.getEP(box_x[i-4],box_y[i-4],vpL.qi,vpL.vps,w,h);

		boundaryLines[i][2]=ep.pts[3*2+0];
		boundaryLines[i][3]=ep.pts[3*2+1];
	}

	for(int i=0;i<8;i++)
	{
		Delta_x=boundaryLines[i][2]-boundaryLines[i][0];
		Delta_y=boundaryLines[i][3]-boundaryLines[i][1];
		lineLength=sqrt(Delta_x*Delta_x+Delta_y*Delta_y);
		Delta_x=2*(J_SUBSAMP*Delta_x)/lineLength;
		Delta_y=2*(J_SUBSAMP*Delta_y)/lineLength;
		for(int j=0;j<(lineLength/(2*J_SUBSAMP));j++)
		{
			x=boundaryLines[i][0]+j*Delta_x;
			y=boundaryLines[i][1]+j*Delta_y;
			junctionYTypes[0]=1;junctionYTypes[1]=1;junctionYTypes[2]=1;junctionYTypes[3]=0;junctionYTypes[4]=0;junctionYTypes[5]=0;
			juncYScore=computeJunctionScore((int)(x/J_SUBSAMP),(int)(y/J_SUBSAMP),junctionYTypes,SHORTJUNC_LEN/3,SHORTJUNC_LEN/6);
			if (juncYScore>juncYThreshold)
				boxScore=boxScore+juncYScore;

			
			junctionYTypes[0]=1;junctionYTypes[1]=1;junctionYTypes[2]=0;junctionYTypes[3]=0;junctionYTypes[4]=0;junctionYTypes[5]=1;
			juncYScore=computeJunctionScore((int)(x/J_SUBSAMP),(int)(y/J_SUBSAMP),junctionYTypes,SHORTJUNC_LEN/3,SHORTJUNC_LEN/6);
			if (juncYScore>juncYThreshold)
				boxScore=boxScore+juncYScore;
				
		}
	}*/

	return boxScore;
}

double 
Junctions::computeJunctionScore(int x,int y,int juncType[6],double highThreshold,double lowThreshold)
{
	double votes[6],junctionScore;
	junctionScore=1;
	for(int i=0;i<6;i++)
		if (junctionScore>0)
		{
			if (juncType[i]==1)
			{
				if (votesFull[6*(dimVotesX*y+x)+i]>highThreshold)
					junctionScore=junctionScore*votesFull[6*(dimVotesX*y+x)+i];
				else
					junctionScore=0;
			}
			else
			{
				if (votesShort[6*(dimVotesX*y+x)+i]>lowThreshold)					
					junctionScore=0;
			}
		}
	return junctionScore;
}

boxLayout 
Junctions::getBoxLayout(double box_x[4],double box_y[4],vp vpL,int w,int h)
{
	boxLayout boxLayout_;
	Votes votes_;
	Points ep;
	geometry geometry_;
	int counter;
	for(int i=0;i<4;i++)
	{
		boxLayout_.boundaryLines[i][0]=box_x[i];
		boxLayout_.boundaryLines[i][1]=box_y[i];
		boxLayout_.boundaryLines[i][2]=box_x[(i+1)%4];
		boxLayout_.boundaryLines[i][3]=box_y[(i+1)%4];
	}
	for(int i=4;i<8;i++)
	{
		boxLayout_.boundaryLines[i][0]=box_x[i-4];
		boxLayout_.boundaryLines[i][1]=box_y[i-4];
		ep=votes_.getEP(box_x[i-4],box_y[i-4],vpL.qi,vpL.vps,w,h);
		boxLayout_.boundaryLines[i][2]=ep.pts[3*2+0];
		boxLayout_.boundaryLines[i][3]=ep.pts[3*2+1];
		delete ep.pts;
	}
	
	// 0-left,1-ceiling,2-right,3-floor,4-middle.
	counter=0;
	int index1=2,index2=3,axis_type=0;
	double x1=0,y1=h,x2=0,y2=0;	
	boxLayout_.poly[0]	=	getPoly(boxLayout_.boundaryLines,2,3,0,1-EPSILON,1+EPSILON,0,h,0,0);
	boxLayout_.poly[1]	=	getPoly(boxLayout_.boundaryLines,3,0,1,1-EPSILON,1+EPSILON,0,0,w,0);
	boxLayout_.poly[2]	=	getPoly(boxLayout_.boundaryLines,0,1,0,w-EPSILON,w+EPSILON,w,0,w,h);
	boxLayout_.poly[3]	=	getPoly(boxLayout_.boundaryLines,1,2,1,h-EPSILON,h+EPSILON,w,h,0,h);
	
	boxLayout_.poly[4].pts=new double[3*4];
	boxLayout_.poly[4].noPts=4;
	for(int i=0;i<4;i++)
	{
		boxLayout_.poly[4].pts[3*i+0]=boxLayout_.boundaryLines[i][0];
		boxLayout_.poly[4].pts[3*i+1]=boxLayout_.boundaryLines[i][1];
	}

	// left wall signature.
	double insideX,insideY;

	for(int i=0;i<5;i++)
	{
		insideX=0;insideY=0;
		for(int j=0;j<boxLayout_.poly[i].noPts;j++)
		{
			insideX=insideX+boxLayout_.poly[i].pts[3*i+0];
			insideY=insideX+boxLayout_.poly[i].pts[3*i+1];
		}
		insideX=insideX/boxLayout_.poly[i].noPts;
		insideY=insideY/boxLayout_.poly[i].noPts;
		for(int j=0;j<boxLayout_.poly[i].noPts;j++)
		{
			boxLayout_.polyLines[i][j]=geometry_.getLineParams(boxLayout_.poly[i].pts[3*j+0],
															   boxLayout_.poly[i].pts[3*j+1],
															   boxLayout_.poly[i].pts[3*((j+1)%boxLayout_.poly[i].noPts)+0],
															   boxLayout_.poly[i].pts[3*((j+1)%boxLayout_.poly[i].noPts)+1]);
			if ((boxLayout_.polyLines[i][j].a*insideX+boxLayout_.polyLines[i][j].b*insideY+boxLayout_.polyLines[i][j].c)>0)
				 boxLayout_.signature[i][j]=0;
			else
				 boxLayout_.signature[i][j]=1;
		}
	}


	return boxLayout_;
}

Points 
Junctions::getPoly(double boundaryLines[8][4],int index1,int index2,int axis_type,double threshold1,double threshold2,double x1,double y1,double x2,double y2)
{
	Points ep;
	ep.pts=new double[3*6];
	int counter=0;
	printf("boundary lines:%lf,%lf\n",boundaryLines[4+index1][axis_type+2],boundaryLines[4+index2][axis_type+2]);

	if ((boundaryLines[4+index1][axis_type+2]>threshold1)&&(boundaryLines[4+index1][axis_type+2]<threshold2))
	{
		ep.pts[3*counter+0]=boundaryLines[4+index1][2];
		ep.pts[3*counter+1]=boundaryLines[4+index1][3];
		counter=counter+1;		
	}
	else
	{
		ep.pts[3*counter+0]=boundaryLines[4+index1][2];
		ep.pts[3*counter+1]=boundaryLines[4+index1][3];
		counter=counter+1;
		ep.pts[3*counter+0]=x1;
		ep.pts[3*counter+1]=y1;
		counter=counter+1;
	}
	
	if ((boundaryLines[4+index2][axis_type+2]>threshold1)&&(boundaryLines[4+index2][axis_type+2]<threshold2))
	{
		ep.pts[3*counter+0]=boundaryLines[4+index2][2];
		ep.pts[3*counter+1]=boundaryLines[4+index2][3];
		counter=counter+1;
	}
	else
	{
		ep.pts[3*counter+0]=x2;
		ep.pts[3*counter+1]=y2;
		counter=counter+1;
		ep.pts[3*counter+0]=boundaryLines[4+index2][2];
		ep.pts[3*counter+1]=boundaryLines[4+index2][3];
		counter=counter+1;		
	}		
	ep.pts[3*counter+0]=boundaryLines[index2][0];
	ep.pts[3*counter+1]=boundaryLines[index2][1];
	counter=counter+1;
	ep.pts[3*counter+0]=boundaryLines[index1][0];
	ep.pts[3*counter+1]=boundaryLines[index1][1];
	counter=counter+1;
	ep.noPts=counter;
	return ep;
}
